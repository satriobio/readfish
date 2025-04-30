"""Mapping interface for readfish, using Minimap2 mappy, or mappy-rs as dictated by the experiment TOML
`mapper_settings.<PLUGIN>` section. See {ref}`plugin configuration <plugins-config>` section.
"""

from __future__ import annotations
from enum import Enum
from itertools import chain, repeat
from pathlib import Path
from typing import Optional, Iterable

from readfish._loggers import setup_logger
from readfish._config import Barcode, Region
from readfish.plugins.abc import AlignerABC
from readfish.plugins.utils import (
    Result,
    Strand,
    count_dict_elements,
    get_contig_lengths,
    _summary_percent_reference_covered,
)
from readfish._utils import nice_join

from collections import defaultdict
from datetime import datetime
import random

class Alignment():
    def __init__(self, ctg):
        self.ctg = ctg
        self.strand = 1
        self.r_en = 0
        self.r_st = 0

class Aligner(AlignerABC):
    """
    A wrapper class for mappy.Aligner providing an interface to either `mappy` or `mappy-rs` aligner implementations.

    :param Aligners leveling_impl: Specifies which mappy implementation to use, represented by an `Aligners` enum instance.
    :param Optional[str] debug_log: Specifies the file to log debug information to. If `None`, no logging is performed.
    :param kwargs: Additional keyword arguments to be passed to the mappy aligner constructor.

    :raises ValueError: When an invalid `leveling_impl` value is provided.

    Example:
        .. code-block:: python

            aligner = _Aligner(leveling_impl=Aligners.C_MAPPY, debug_log="debug.log", preset="map-ont")

    .. note::
        The actual aligner instance is created during initialization based on the provided `leveling_impl`, and is accessible via the `aligner` attribute of the instance.

    .. seealso::
        {ref}`plugin configuration <plugins-config>` for details on how `leveling_impl` is determined from the experiment TOML configuration.
    """

    def __init__(self, debug_log: Optional[str] = None, **kwargs):
        self.leveling_params = kwargs
        self.barcode_db = defaultdict(int)
        self.barcode_db_unblock = defaultdict(int)

    def validate(self) -> None:
        """
        Check that this aligner can be initialised without any issues. Catches any problems and raises helpful errors.
        Currently checks:
         1. that the Reference (fn_idx_in) exists, IF one is provided
         2. That the reference is an .mmi file, or a FASTA or FASTQ, either uncompressed or Gzipped, IF a fn_idx_in is provided.

        :return: None, if the Aligner is setup with valid paths and permissions
        """
        pass

    def disconnect(self) -> None:
        return
    
    def leveling_enrichment(self, barcode_db, min_reads=1000, upper_threshold=1.5, lower_threshold=0.5):
        """Identify usable barcodes based on enrichment/depletion analysis.
        
        Args:
            barcode_db: Dictionary {barcode: count} of read counts
            min_reads: Minimum absolute read count to consider a barcode
            upper_threshold: Fold-enrichment above which barcode is always unblocked
            lower_threshold: Fold-enrichment below which barcode is considered depleted
            
        Returns:
            tuple: (unblocked_barcodes, is_depleted)
                unblocked_barcodes: List of barcodes passing filters
                is_depleted: Boolean indicating if depletion was detected
        """
        if not barcode_db:
            return [], False

        total_counts = sum(barcode_db.values())
        barcodes_num = len(barcode_db)

        if barcodes_num == 0 or total_counts == 0:
            return [], False

        expected_count = total_counts / barcodes_num
        
        # Calculate enrichment and filter by min_reads
        enrichment_scores = {
            barcode: count / expected_count
            for barcode, count in barcode_db.items()
            if count >= min_reads
        }
        
        # Check for depletion (any barcode below lower threshold)
        is_depleted = any(
            count < (lower_threshold * expected_count)
            for count in barcode_db.values()
        )
        
        # Determine unblocked barcodes
        if is_depleted:
            # When depletion exists, use barcodes at or above expected level
            unblocked_barcodes = [
                bc for bc, score in enrichment_scores.items()
                if score >= 1.0
            ]
        else:
            # Otherwise use barcodes above upper threshold (if any)
            unblocked_barcodes = [
                bc for bc, score in enrichment_scores.items()
                if score >= upper_threshold
            ]
            # If none meet upper threshold, use all passing min_reads
            if not unblocked_barcodes:
                unblocked_barcodes = list(enrichment_scores.keys())
        
        return unblocked_barcodes
    
    def make_decision(self, reads):
        """Process barcode calls and make alignment decisions, logging cumulative barcode stats.
        
        Args:
            reads: List of (channel, read) tuples
            calls: List of (barcode, mapq, baseq) tuples
            
        Yields:
            Tuple of (channel, read, barcode, mapq, baseq, aligner)
        """
        timestamp = datetime.now().isoformat()

        reads = [x for x in reads if x is not None]
        
        batch_stats = {
            'batch_size': len(reads),
            'qc_passed': 0,
            'target_reads': 0,
            'qc_failed': 0
        }

        # Setup logging files
        with open("live_barcode_stats.tsv", "a") as stats_file, \
             open("live_barcode_batch.tsv", "a") as batch_file, \
             open("live_unblock_stats.tsv", "a") as unblock_file, \
             open("live_read_decisions.tsv", "a") as decision_file:

            # Initialize headers if files are empty
            if stats_file.tell() == 0:
                stats_file.write("timestamp\ttotal_reads\tbarcode\tcount\tproportion\n")
            if batch_file.tell() == 0:
                batch_file.write("timestamp\ttotal_reads\tqc_passed\tqc_failed\ttarget\n")
            if decision_file.tell() == 0:
                decision_file.write("timestamp\tchannel\tread_id\tbarcode\tmapq\tbaseq\tclassification\n")

            total_reads = sum(self.barcode_db.values())

            # Filter valid target barcodes once
            barcode_db_target = {
                bc: self.barcode_db[bc]
                for bc in self.barcode_db.keys()
                if 'target_barcodes' not in self.leveling_params or bc in self.leveling_params['target_barcodes']
            }
            total_reads_valid = sum(barcode_db_target.values())

            for barcode, count in self.barcode_db.items():
                if 'target_barcodes' in self.leveling_params and barcode in self.leveling_params['target_barcodes']:
                    proportion = count / total_reads_valid if total_reads_valid > 0 else 0
                else:
                    proportion = "not a valid barcode"

                stats_file.write(
                    f"{timestamp}\t{total_reads}\t{barcode}\t{count}\t"
                    f"{proportion if isinstance(proportion, str) else f'{proportion:.4f}'}\n"
                )

            # Find barcode with highest frequency
            unblocked_barcodes = self.leveling_enrichment(barcode_db_target, upper_threshold=self.leveling_params['upper_threshold'], lower_threshold=self.leveling_params['lower_threshold'])
            
            for read in reads:
                channel = read.channel
                read_id = read.read_id
                barcode = str(read.basecall_data[0])
                mapq = read.basecall_data[1]
                baseq = read.basecall_data[2]

                # Filter calls based on baseq
                if baseq > self.leveling_params['min_baseq']:
                    batch_stats['qc_passed'] += 1

                    # Label significantly frequent barcode as target to unblock
                    # classification = "TARGET" if barcode in unblocked_barcodes else "NONTARGET"
                    classification = "TARGET" if (barcode in unblocked_barcodes and random.random() < self.leveling_params['likelihood']) else "NONTARGET"

                    if classification == "TARGET":
                        batch_stats['target_reads'] += 1
                        self.barcode_db_unblock[barcode] += 1
                    else:
                        self.barcode_db[barcode] += 1

                    aln = Alignment(ctg=classification)

                    # Log decision
                    decision_file.write(
                        f"{timestamp}\t{channel}\t{read_id}\t{barcode}\t"
                        f"{mapq}\t{baseq}\t{classification}\n"
                    )

                    yield (channel, read, barcode, mapq, baseq, aln)
                else:
                    batch_stats['qc_failed'] += 1
                #     decision_file.write(
                #         f"{timestamp}\t{channel}\t{read.id}\t{barcode}\t"
                #         f"{mapq}\t{baseq}\tNONE\tqc_failed\n"
                #     )

            # Log batch summary
            batch_file.write(
                f"{timestamp}\t"
                f"{batch_stats['batch_size']}\t"
                f"{batch_stats['qc_passed']}\t"
                f"{batch_stats['qc_failed']}\t"
                f"{batch_stats['target_reads']}\n"
            )

            # for barcode, prob in barcode_db_prob.items():
            #     unblock_file.write(
            #         f"{timestamp}\t{barcode}\t{prob}\n"
            #     )
    
    @property
    def initialised(self) -> bool:
        return True

    def map_reads(self, basecall_results: Iterable[Result]) -> Iterable[Result]:
        """
        Make decision for an iterable of base-called data using either the `` or `` implementation,
        based on the `leveling_impl` provided during the instantiation of the class.

        :param Iterable[Result] basecall_results: An iterable of basecalled read results to be mapped.
        :return: An iterable of mapped read results.

        :raises NotImplementedError: If the aligner is not configured (i.e., if `leveling_impl` is neither `Strategies.` nor `Strategies.`).

        Example:
            .. code-block:: python

                mapped_results = aligner.map_reads(basecall_results)

        """
        results = self.make_decision(basecall_results)

        for channel, read, barcode, mapq, baseq, aln in results:
            yield Result(
                channel=channel,
                read_id=read.read_id,
                seq=[],
                barcode=None,
                basecall_data=(barcode, mapq, baseq),
                alignment_data=[aln],
            )

    def describe(self, regions: list[Region], barcodes: list[Barcode]) -> str:
        """
        Describe the leveling plugin

        :return: Description of parameters passed to this leveling plugin
        """
        description = ["Utilising the Leveling plugin:"]
        return "\n".join(description)