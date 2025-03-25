# Core imports
from __future__ import annotations
import argparse
import traceback
from packaging.version import Version
import logging
import time
from timeit import default_timer as timer
from pathlib import Path
from typing import Any

# Third party imports
from readfish.read_until.read_cache import AccumulatingCache
from readfish.read_until import ReadUntilClient
from minknow_api import protocol_service

# Library
from readfish._cli_args import DEVICE_BASE_ARGS, Chemistry
from readfish._read_until_client import RUClient
from readfish._config import Action, Conf, make_decision, _Condition
from readfish._statistics import ReadfishStatistics
from readfish.__about__ import __version__
from readfish._compatibility import (
    _get_minknow_version,
    check_compatibility,
    MINKNOW_COMPATIBILITY_RANGE,
    DIRECTION,
)
from readfish._utils import (
    get_device,
    send_message,
    ChunkTracker,
    Severity,
)
from readfish.plugins.abc import AlignerABC, CallerABC
from readfish.plugins.utils import (
    Decision,
    PreviouslySentActionTracker,
    Result,
    DuplexTracker,
    Strand,
)

from readfish.plugins.utils import Targets, Action, Decision, Result

_help = "Run targeted sequencing"
_cli = DEVICE_BASE_ARGS + (
    (
        "--toml",
        dict(
            metavar="TOML",
            required=True,
            help="TOML file specifying experimental parameters",
        ),
    ),
    (
        "--no-debug-log",
        dict(
            help="Disable debug output of information about chunks seen into a .tsv formatted log. Default enabled.",
            action="store_false",
            dest="debug_log",
        ),
    ),
    (
        "--padding",
        dict(
            help="Number of bases to pad the target sequences with",
            default=0,
            type=int,
            metavar="PADDING",
        ),
    ),
)
# When sequencing in duplex mode, overriding a decided `Action` on a currently sequenced molecule
# is not allowed if the previous molecules decision was one of these.
DISALLOWED_DUPLEX_DECISIONS = {Decision.first_read_override, Decision.duplex_override}


class Analysis:
    """
    Analysis class where the read until magic happens. Comprises of one run
    function that is run threaded in the run function at the base of this file.
    Arguments listed in the __init__ docs.

    :param client: An instance of the ReadUntilClient object.
    :param conf: An instance of the Conf object.
    :param logger: The command level logger for this module.
    :param debug_log: Whether to output the Debug Log. log Name is generated.
    :param throttle: The time interval (seconds) between requests to the ReadUntilClient.
    :param unblock_duration: Time, in seconds, to apply unblock voltage.
    :param dry_run: If True unblocks are replaced with `stop_receiving` commands.
    :param toml: The path to the toml file containing experiment conf. Used as the path for checking if the TOML needs reloading.
    :param chemistry: Instance of Chemistry Enum, representing the chemistry of the run (Simplex/Duplex). Used for
        decision making on strands that may be part of a duplex pair.
    """

    def __init__(
        self,
        client: ReadUntilClient,
        conf: Conf,
        logger: logging.Logger,
        debug_log: bool,
        throttle: float,
        unblock_duration: float,
        dry_run: bool,
        toml: str,
        chemistry: Chemistry,
    ):
        self.client = client
        self.conf = conf
        self.logger = logger
        self.debug_log = debug_log
        self.throttle = throttle
        self.unblock_duration = unblock_duration
        self.dry_run = dry_run
        self.live_toml = Path(f"{toml}_live").resolve()
        self.run_information = self.client.connection.protocol.get_run_info()
        self.chemistry = chemistry
        # Generate a run specific read log
        read_log_name = (
            f"{self.run_information.run_id}_readfish.tsv" if debug_log else None
        )

        self.logger.info("Fetching Run Configuration")
        self.break_reads_after_seconds = (
            self.client.connection.analysis_configuration.get_analysis_configuration().read_detection.break_reads_after_seconds.value
        )
        self.sample_rate = self.client.connection.device.get_sample_rate().sample_rate
        self.logger.info("Run Configuration Received")
        self.logger.info(f"run_id={self.run_information.run_id}")
        self.logger.info(f"break_reads_after_seconds={self.break_reads_after_seconds}")
        self.logger.info(f"sample_rate={self.sample_rate}")
        # Create our statistics tracker
        self.loop_statistics = ReadfishStatistics(
            read_log_name, self.break_reads_after_seconds
        )
        logger.info("Initialising Caller")
        self.caller: CallerABC = self.conf.caller_settings.load_object(
            "Caller", run_information=self.run_information, sample_rate=self.sample_rate
        )
        logger.info("Caller initialised")
        caller_description = self.caller.describe()
        self.logger.info(caller_description)
        send_message(self.client.connection, caller_description, Severity.INFO)
        logger.info("Initialising Aligner")
        self.mapper: AlignerABC = self.conf.mapper_settings.load_object("Aligner")
        self.logger.info("Aligner initialised")
        # count how often a read is seen
        self.chunk_tracker = ChunkTracker(self.client.channel_count)

        # This is an object to keep track of the last action sent to the client for each channel
        self.previous_action_tracker = PreviouslySentActionTracker()
        # Keep track of previous alignments
        self.duplex_tracker = DuplexTracker()

        # We assume that sequencing is already running.
        # If the run is not in sequencing phase when the read until loop starts will
        # be set to false and the first read seen can be unblocked.
        self.readfish_started_during_sequencing = True

        # This is a flag to prevent repeated logging of the same message
        self.log_once_in_loop = True

    @property
    def wait_for_sequencing(self) -> bool:
        """
        Wait for minKNOW to report PHASE_SEQUENCING before starting readfish tight loop.
        The check occurs in out RUClient wrapper.

        :return: True if we are waiting for PHASE_SEQUENCING, False otherwise

        """
        if self.client.wait_for_sequencing_to_start:
            if self.log_once_in_loop:
                self.logger.info(
                    f"MinKNOW is reporting {protocol_service.ProtocolPhase.Name(self.client.current_protocol_phase)}, waiting for PHASE_SEQUENCING to begin."
                )
                self.log_once_in_loop = not self.log_once_in_loop
            self.readfish_started_during_sequencing = False  # We are not in sequencing phase, so we can unblock the first read we see as we will be sequencing it from the start
            return True
        return False

    def reload_toml(self, last_toml_mtime: float) -> float:
        """
        Reload the toml to refresh the conf with any updates.
        Reloading is determined by checking the modified time of the toml file.
        If it is more recent, reload the conf.

        :param last_live_mtime: The last modified time for the toml file.

        :return: The last modified time for the toml file, updated if changed.
        """
        if (
            self.live_toml.is_file()
            and self.live_toml.stat().st_mtime > last_toml_mtime
        ):
            try:
                self.conf = Conf.from_file(
                    self.live_toml, self.client.channel_count, self.logger
                )
            # FIXME: Broad exception
            except Exception as e:
                if hasattr(e, "exceptions"):
                    self.logger.error(getattr(e, "exceptions"))
                    self.logger.error(traceback.format_exc())
            last_toml_mtime = self.live_toml.stat().st_mtime
        return last_toml_mtime

    def check_override_action(
        self,
        control: bool,
        action: Action,
        result: Result,
        seen_count: int,
        condition: _Condition,
        stop_receiving_action_list: list[tuple[int, int]],
        unblock_batch_action_list: list[tuple[int, int]],
    ) -> tuple[Action, bool, str | None]:
        """
        Check the chosen Action and amend it based on conditional checks.
        The action lists are appended to in place, so no return is required.

        Checks include:
            1. If the read is in a control region, the action is always stop_receiving.
            1. If the read is below the minimum chunks, use value in toml or default to proceed
            1. If the read is above the maximum chunks, use value in toml or default unblock--throttle
            1. First read seen for channel and readfish started during sequencing, override to stop_receiving
            1. If action is unblock and we are dry-running, override to stop_receiving
            1. If we are running in duplex chemistry, check the previous reads final decision and Action, and potentially sequence
                the current read, instead of unblocking it.

        :param control: Indicates read from a channel in a control region
        :param action: What action was decided for this read before any meddling
        :param result: Information about the current read.
        :param seen_count: Number of times other chunks from the read have been observed.
        :param condition: The set of conditions for deciding the action.
        :param stop_receiving_action_list: List to append channels and read numbers for which 'stop receiving' action is decided.
        :param unblock_batch_action_list: List to append channels, read numbers, and read IDs for which 'unblock' action is decided.

        :return: A tuple containing the previous action taken for this read,
          boolean indicating if the action was overridden, and the name of the action overridden too.

        """

        # Easy dub
        if control:
            action = Action.stop_receiving
        else:
            # TODO: Document the less than logic here
            below_min_chunks = seen_count < condition.min_chunks
            above_max_chunks = seen_count > condition.max_chunks

            # TODO: This will also factor into the precedence and documentation
            # If we have seen this read more than the max chunks and want to
            #   evaluate it again (Action.proceed) then we will overrule that
            #   action using the above_max_chunks_action, unblock by default
            if above_max_chunks and action is Action.proceed:
                action = condition.above_max_chunks
                result.decision = Decision.above_max_chunks

            # If we are below min chunks and we get an action that is not PROCEED
            #   then we will overrule that action using the below_min_chunks_action
            #   which by default is proceed.
            if below_min_chunks and action is not Action.proceed:
                action = condition.below_min_chunks
                result.decision = Decision.below_min_chunks

        # previous_action will be None if the read has not been seen before.
        previous_action = self.previous_action_tracker.get_action(result.channel)
        action_overridden = False
        # If --duplex flag override decisions made based on the strand and contig alignment of the previous read.
        # Unfinished bruv
        if (
            self.chemistry is Chemistry.DUPLEX
            # Easy checks first, so wdon't do more complex processing unless we have to
            and action == Action.unblock
            and previous_action is Action.stop_receiving
        ):
            # Check if we think this read is possibly duplex
            possible_duplex = any(
                self.duplex_tracker.possible_duplex(
                    result.channel, result.read_id, al.ctg, al.strand
                )
                for al in result.alignment_data
            )
            # Check the previous decision for this channel was not already an override
            previous_decision_allowed = (
                self.duplex_tracker.get_previous_decision(result.channel)
                not in DISALLOWED_DUPLEX_DECISIONS
            )
            if possible_duplex and previous_decision_allowed:
                self.logger.debug(
                    f"Overriding read {result.read_id} as it is possibly second half of a duplex"
                    f"- previous read action {previous_action}, current_action: {action},"
                    f" previous_decision: {self.duplex_tracker.get_previous_decision(result.channel)}"
                )
                action_overridden = True
                result.decision = Decision.duplex_override
                action = Action.stop_receiving
        # Duplex
        elif (
            self.chemistry is Chemistry.DUPLEX_SIMPLE
            and previous_action is Action.stop_receiving
            and action is Action.unblock
        ):
            previous_decision_allowed = (
                self.duplex_tracker.get_previous_decision(result.channel)
                not in DISALLOWED_DUPLEX_DECISIONS
            )
            if previous_decision_allowed:
                self.logger.debug(
                    f"Overriding to duplex - previous read action {previous_action}, current_action: {action},"
                    f" previous_decision: {self.duplex_tracker.get_previous_decision(result.channel)}"
                )
                action = Action.stop_receiving
                action_overridden = True
                result.decision = Decision.duplex_override

        # Override to stop receiving if this is the first read ona channel and we started mid sequencing
        if previous_action is None and self.readfish_started_during_sequencing:
            self.logger.debug(
                f"This is the first suitable read chunk from channel {result.channel}. Translocated read length unknown, sequencing."
            )
            action_overridden = True
            result.decision = Decision.first_read_override
            action = Action.stop_receiving

        if action is Action.stop_receiving:
            stop_receiving_action_list.append((result.channel, result.read_id))

        elif action is Action.unblock:
            if self.dry_run:
                # Log an 'unblock' action to previous action, but send a 'stop receiving' to prevent further read processing.
                action_overridden = True
                stop_receiving_action_list.append((result.channel, result.read_id))
            else:
                unblock_batch_action_list.append((result.channel, result.read_id))

        # If we have made a final decision for this read and we shouldn't see it again!
        if action is Action.unblock or action is Action.stop_receiving:
            # Add decided Action
            self.previous_action_tracker.add_action(result.channel, action)
            # Add duplex based tracking if we are in duplex mode
            if self.chemistry is Chemistry.DUPLEX_SIMPLE:
                self.duplex_tracker.set_decision(result.channel, result.decision)
            elif self.chemistry is Chemistry.DUPLEX:
                self.duplex_tracker.set_decision(result.channel, result.decision)
                self.duplex_tracker.set_alignments(
                    result.channel,
                    [(al.ctg, Strand(al.strand)) for al in result.alignment_data],
                )

        return (
            action,
            previous_action,
            action_overridden,
            action.name if action_overridden else None,
        )

    def run(self):
        """Run the read until loop, in one continuous while loop."""

        # TODO: Swap this for a CSV record later
        self.conf.write_channels_toml(self.client.mk_run_dir)

        # TODO: This could still be passed through to the basecaller to prevent
        #       rebasecalling data that is already being unblocked or sequenced
        loop_counter = 0

        last_live_toml_mtime = 0
        self.logger.info("Starting main loop")
        self.logger.info("Generating aligner description, if possible...")
        mapper_description = self.mapper.describe(self.conf.regions, self.conf.barcodes)
        self.logger.info(mapper_description)
        send_message(self.client.connection, mapper_description, Severity.INFO)

        while self.client.is_sequencing:
            t0 = timer()
            # Check if we have started readfish before PHASE_SEQUENCING,
            if self.wait_for_sequencing:
                time.sleep(self.throttle)
                continue
            # Set back to true for when we re-enter a non sequencing phase
            self.log_once_in_loop = True

            if self.readfish_started_during_sequencing and loop_counter == 0:
                self.logger.info(
                    "readfish started in PHASE_SEQUENCING. Fully sequencing first read from each channel."
                )
            if not self.mapper.initialised:
                self.logger.warning(
                    "readfish main loop started but mapper is not initialised. Please check your aligners plugin documentation."
                    "If you are using mappy or mappy-rs this is definitely an error. Please open an issue here - "
                    "https://github.com/LooseLab/readfish/issues"
                )
                time.sleep(self.throttle)
                continue

            last_live_toml_mtime = self.reload_toml(last_live_toml_mtime)
            ########### Main Loop ###########
            loop_counter += 1
            number_reads = 0
            unblock_batch_action_list = []
            stop_receiving_action_list = []

            chunks = self.client.get_read_chunks(self.client.channel_count, last=True)
            
            for channel, read in chunks:
                number_reads += 1

                read_id = read_id
                raw_data = np.frombuffer(read.raw_data, self.client.signal_dtype)

                ## processing
                # channel, read_id, raw_data

                ## unblock_batch_action_list
                # List to append channels, read numbers, and read IDs for which 'unblock' action is decided.
                # example:
                # unblock_batch_action_list = [(channel, read_id), (channel, read_id), (channel, read_id)]

                unblock_batch_action_list.append((channel, read_id))
            
            #######################################################################
            # Compile actions to be sent
            self.client.unblock_read_batch(
                unblock_batch_action_list, duration=self.unblock_duration
            )
            self.client.stop_receiving_batch(stop_receiving_action_list)

            t1 = timer()
            if number_reads > 0:
                self.loop_statistics.add_batch_performance(
                    number_of_reads=number_reads, batch_time=t1 - t0
                )
                self.logger.info(self.loop_statistics.get_batch_performance())

            # limit the rate at which we make requests
            if t0 + self.throttle > t1:
                time.sleep(self.throttle + t0 - t1)
        else:
            send_message(
                self.client.connection,
                "Readfish client stopped.",
                Severity.WARN,
            )
            self.caller.disconnect()
            self.mapper.disconnect()
            self.logger.info("Finished analysis of reads as client stopped.")


def run(
    parser: argparse.ArgumentParser, args: argparse.ArgumentParser, extras: list[Any]
) -> int:
    """Run function for targets.py

    Imported in `_cli_base.py`.
    Sets up the read until client and starts the analysis thread above.

    :param parser: Argparse onject - unused but must be taken due as may be needed
    :param args: The arguments passed to ArgParse
    :param extras: Extra stuff, I guess

    :returns: An exit code integer, 0 for success
    """
    # Setup logger used in this entry point, this one should be passed through
    logger = logging.getLogger(f"readfish.{args.command}")

    # Check MinKNOW version

    minknow_version = _get_minknow_version(host=args.host, port=args.port)
    if (
        action := check_compatibility(minknow_version, MINKNOW_COMPATIBILITY_RANGE)
    ) in (
        DIRECTION.UPGRADE,
        DIRECTION.DOWNGRADE,
    ):
        lower_bound, upper_bound = MINKNOW_COMPATIBILITY_RANGE
        logger.warning(
            f"""This readfish version ({__version__}) is tested for compatibility with MinKNOW v{lower_bound} to v{upper_bound}.
This version of minknow is {minknow_version}.
If readfish fails please try to {action.value} readfish.
If there isn't a newer version of readfish and readfish is failing, please open an issue:
    https://github.com/LooseLab/readfish/issues"""
        )

    if minknow_version < Version("6.0.0"):
        logger.critical(
            f"This version of readfish ({__version__}) is not compatible with less than MinKNOW 6.X.X, you downgrade to at least readfish 2024.2.0"
            f"This won't work, exiting..."
        )
        raise SystemExit(1)

    # Fetch sequencing device
    position = get_device(args.device, host=args.host, port=args.port)

    # Create a read until client
    read_until_client = RUClient(
        mk_host=position.host,
        mk_port=position.description.rpc_ports.secure,
        filter_strands=True,
        cache_type=AccumulatingCache,
        timeout=args.wait_for_ready,
        prefilter_classes={
            "strand",
            "strand2",
            "short_strand",
            "adapter",
            "unknown_positive",
        },
    )

    # Load TOML configuration
    conf = Conf.from_file(args.toml, read_until_client.channel_count, logger=logger)

    # Set the padding if it is specified.
    if padding := getattr(args, "padding", None):
        for region in conf.regions:
            region.targets.padding = padding
        for barcode in conf.barcodes:
            conf.barcodes[barcode].targets.padding = padding
    logger.info(conf.describe_experiment())

    send_message(
        read_until_client.connection,
        f"'readfish {args.command}' connected to this device.",
        Severity.WARN,
    )

    # start the client running
    read_until_client.run(
        first_channel=1,
        last_channel=read_until_client.channel_count,
        max_unblock_read_length_seconds=args.max_unblock_read_length_seconds,
        accepted_first_chunk_classifications=[
            "strand",
            "strand2",
            "short_strand",
            "adapter",
            "unknown_positive",
        ],
    )

    worker = Analysis(
        read_until_client,
        conf=conf,
        logger=logger,
        debug_log=args.debug_log,
        unblock_duration=args.unblock_duration,
        throttle=args.throttle,
        dry_run=args.dry_run,
        toml=args.toml,
        chemistry=Chemistry(args.chemistry),
    )

    # begin readfish function
    try:
        worker.run()
    except KeyboardInterrupt:
        logger.info("Keyboard interrupt received, stopping readfish.")
        pass
    finally:
        read_until_client.reset()

    send_message(
        read_until_client.connection,
        "Readfish disconnected from this device. Sequencing will proceed normally.",
        Severity.WARN,
    )
    return 0
