"""Seqtagger plugin module

"""

from __future__ import annotations
import logging
import os
import time
from collections import namedtuple
from pathlib import Path
from typing import Iterable, TYPE_CHECKING
from packaging.version import parse as parse_version

import numpy as np
import numpy.typing as npt
from minknow_api.protocol_pb2 import ProtocolRunInfo
from minknow_api.device_pb2 import GetSampleRateResponse

from readfish._loggers import setup_logger
from readfish.plugins.abc import CallerABC
from readfish.plugins.utils import Result
from readfish._utils import nice_join

if TYPE_CHECKING:
    import minknow_api

__all__ = ["Caller"]

logger = logging.getLogger("RU_basecaller")
CALIBRATION = namedtuple("calibration", "scaling offset")

from datetime import datetime

class DefaultDAQValues:
    """Provides default calibration values

    Mimics the read_until_api calibration dict value from
    https://github.com/nanoporetech/read_until_api/blob/2319bbe/read_until/base.py#L34
    all keys return scaling=1.0 and offset=0.0
    """

    calibration = CALIBRATION(1.0, 0.0)

    def __getitem__(self, _):
        return self.calibration


_DefaultDAQValues = DefaultDAQValues()

from collections import Counter
from demux import load_demux_model, chunks2barcodes, get_med_mad

from collections import defaultdict
from collections import Counter
import os

class Aligner():
    def __init__(self, ctg):
        self.ctg = ctg
        self.strand = 1
        self.r_en = 0
        self.r_st = 0

class Caller():
    def __init__(
        self,
        run_information: ProtocolRunInfo = None,
        sample_rate: GetSampleRateResponse = None,
        debug_log=None,
        **kwargs,
    ):
        self.logger = setup_logger("readfish_seqtagger_logger", log_file=debug_log)
        self.supported_barcode_kits = None
        self.supported_basecall_models = None
        self.run_information = run_information

        self.seqtagger_params = kwargs
        # GET batch size

        if sample_rate:
            self.sample_rate = float(sample_rate)
        else:
            self.sample_rate = float(5000)

        self.validate()
        
        self.demux_model = load_demux_model(self.seqtagger_params['model'], self.seqtagger_params['batchsize'], pwd=self.seqtagger_params['model_secret'])
        self.demux_chunksize = self.demux_model.config['basecaller']['chunksize']

        self.barcode_db = defaultdict(int)
        self.barcode_db_unblock = defaultdict(int)

    def validate(self) -> None:
        """Validate the parameters passed to SeqTagger

        Currently checks:
            1. 

        :return: None, if the parameters pass all the checks

        """
        model = self.seqtagger_params.get("model", "").strip()
        model_path = Path(f"{model}.zip")

        if not model_path.is_file():
            raise FileNotFoundError(f"Model file '{model_path}' does not exist.")

        if model_path.suffix.lower() != ".zip":
            raise RuntimeError(f"Provided model file '{model_path}' is of an incorrect type; expected a .zip file.")
    
    def is_sig(self, barcode, min_reads=10, min_fold_enrichment=2.0):
        """Determine if a barcode count is significantly higher than others
        
        Args:
            barcode: Barcode to test
            min_reads: Minimum absolute read count to consider
            min_fold_enrichment: Minimum fold-enrichment over expected uniform distribution
            
        Returns:
            bool: Whether the barcode is significantly enriched
        """
        if not self.barcode_db or barcode not in self.barcode_db:
            return False
            
        top_count = self.barcode_db[barcode]
        total_counts = sum(self.barcode_db.values())
        
        # Filter valid barcodes (if you have such a list)
        valid_barcodes = [bc for bc in self.barcode_db.keys() if bc in self.valid_barcodes] \
                        if hasattr(self, 'valid_barcodes') else list(self.barcode_db.keys())
        
        if not valid_barcodes:
            return False
            
        barcodes_num = len(valid_barcodes)
        expected_uniform_count = total_counts / barcodes_num
        
        return (top_count >= min_reads and 
                top_count >= min_fold_enrichment * expected_uniform_count)

    def make_decision(self, reads, calls):
        """Process barcode calls and make alignment decisions, logging cumulative barcode stats.
        
        Args:
            reads: List of (channel, read) tuples
            calls: List of (barcode, mapq, baseq) tuples
            
        Yields:
            Tuple of (channel, read, barcode, mapq, baseq, aligner)
        """
        timestamp = datetime.now().isoformat()
        
        batch_stats = {
            'batch_size': len(reads),
            'qc_passed': 0,
            'target_reads': 0,
            'qc_failed': 0
        }

        # Setup logging files
        with open("live_barcode_stats.tsv", "a") as stats_file, \
             open("live_barcode_batch.tsv", "a") as batch_file, \
             open("live_read_decisions.tsv", "a") as decision_file:

            # Initialize headers if files are empty
            if stats_file.tell() == 0:
                stats_file.write("timestamp\ttotal_reads\tbarcode\tcount\tproportion\n")
            if batch_file.tell() == 0:
                batch_file.write("timestamp\ttotal_reads\tqc_passed\tqc_failed\ttarget\n")
            if decision_file.tell() == 0:
                decision_file.write("timestamp\tchannel\tread_id\tbarcode\tmapq\tbaseq\tclassification\tqc\n")

            # Log current proportion stats before processing this batch
            for barcode in self.barcode_db.keys():
                stats_file.write(
                    f"{timestamp}\t{0}\t{barcode}\t{self.barcode_db[barcode]}\t"
                    f"{(self.barcode_db[barcode]/sum(self.barcode_db.values())):.4f}\n"
                )

            # Find barcode with high frequency
            top_barcode = max(self.barcode_db, key=self.barcode_db.get) if self.barcode_db else None
            top_is_sig = self.is_sig(top_barcode) if top_barcode else False           

            for (channel, read), (barcode, mapq, baseq) in zip(reads, calls):
                barcode = str(barcode)

                # Filter calls based on baseq
                if baseq > self.seqtagger_params['min_baseq']:
                    batch_stats['qc_passed'] += 1

                    # Label significantly frequent barcode as target to unblock
                    classification = "TARGET" if (barcode == top_barcode and top_is_sig) else "NONTARGET"
                    if classification == "TARGET":
                        batch_stats['target_reads'] += 1
                        self.barcode_db_unblock[barcode] += 1
                    else:
                        self.barcode_db[barcode] += 1

                    aln = Aligner(ctg=classification)

                    # Log decision
                    decision_file.write(
                        f"{timestamp}\t{channel}\t{read.id}\t{barcode}\t"
                        f"{mapq}\t{baseq}\t{classification}\tqc_passed\n"
                    )

                    yield (channel, read, barcode, mapq, baseq, aln)
                # else:
                #     batch_stats['qc_failed'] += 1
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

    def preprocess_read(self, read: np.ndarray, chunksize: int, offset: int) -> np.ndarray:
        """
        Normalize the first `chunksize` values of `read` using median and MAD.
        
        Parameters:
            read (np.ndarray): The raw signal data.
            chunksize (int): The number of data points to process.
        
        Returns:
            np.ndarray: The normalized signal.
        """
        median, mad = get_med_mad(read[offset:chunksize+offset])
        normalized_signal = (read[:chunksize] - median) / mad
        return normalized_signal

    def basecall(
        self,
        reads: Iterable[tuple[int, minknow_api.data_pb2.GetLiveReadsResponse.ReadData]],
        signal_dtype: npt.DTypeLike,
        daq_values: dict[int, namedtuple] = None,
    ):
        """Call live data from minknow RPC

        :param reads: List or generator of tuples containing (channel, MinKNOW.rpc.Read)
        :param signal_dtype: Numpy dtype of the raw data
        :param daq_values: Dictionary mapping channel to offset and scaling values.
                           If not provided default values of 1.0 and 0.0 are used.
        :yield:
        :rtype: readfish.plugins.utils.Result
        """

        # reads is a list of (channel, read) tuples
        batch_size = self.seqtagger_params['batchsize']
        signal_array = np.zeros((batch_size, self.demux_chunksize), dtype=np.float16)

        # Process all reads in one go
        for channel, read in reads:
            signal = np.frombuffer(read.raw_data, dtype=signal_dtype)
            signal_pA = (signal + daq_values[channel].offset) * daq_values[channel].scaling
            signal_array[channel] = self.preprocess_read(signal_pA, self.demux_chunksize, offset=0)

        # Initial barcode processing
        calls = chunks2barcodes(self.demux_model, signal_array)
        calls = calls[0:len(reads)]

        # # Check which reads need reprocessing with different offset
        # resubmit_indices = [i for i, (_, _, baseq) in enumerate(calls) if baseq < 50]

        # # Only process channels that need retrying
        # if resubmit_indices:  # Only if we have reads to retry
        #     signal_array_retry = np.zeros_like(signal_array)
        #     for idx, (channel, read) in enumerate(reads):
        #         if idx in resubmit_indices:
        #             signal = np.frombuffer(read.raw_data, dtype=signal_dtype)
        #             signal_pA = (signal + daq_values[channel].offset) * daq_values[channel].scaling
        #             signal_array_retry[channel] = self.preprocess_read(
        #                 signal_pA, 
        #                 self.demux_chunksize, 
        #                 offset=self.demux_chunksize//2
        #             )
            
        #     # Second barcode processing for retried signals
        #     retry_calls = chunks2barcodes(self.demux_model, signal_array_retry)
            
        #     # Update only the calls that were retried
        #     for idx in resubmit_indices:
        #         calls[idx] = retry_calls[idx]

        # Combine results and make final decision
        results = self.make_decision(reads, calls)

        for channel, read, barcode, mapq, baseq, aln in results:  # Fixed syntax and iteration logic
            yield Result(
                channel=channel,
                read_id=read.id,
                seq=[],
                barcode=None,
                basecall_data=(barcode, mapq, baseq),
                alignment_data=[aln],
            )
        
    def describe(self) -> str:
        """
        Describe the SeqTAgger Caller

        :return: Description of parameters passed to this SeqTagger plugin
        """
        description = ["Utilising the SeqTagger plugin:"]
        for param in self.seqtagger_params.keys():
            description.append(f"\t- {param}: {self.seqtagger_params[param]}")
        return "\n".join(description)
