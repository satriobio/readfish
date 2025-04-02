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
        self.demux_chunksize_max = self.demux_chunksize*2 

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
    
    def is_sig(self, barcode, barcode_db, min_reads=10, min_fold_enrichment=2.0):
        """Determine if a barcode count is significantly higher than others
        
        Args:
            barcode: Barcode to test
            min_reads: Minimum absolute read count to consider
            min_fold_enrichment: Minimum fold-enrichment over expected uniform distribution
            
        Returns:
            bool: Whether the barcode is significantly enriched
        """
        top_count = barcode_db.get(barcode, 0)  # Ensure barcode exists in filtered list
        total_counts = sum(barcode_db.values())
        barcodes_num = len(barcode_db)
        
        expected_uniform_count = total_counts / barcodes_num

        return (top_count >= min_reads and 
                top_count >= min_fold_enrichment * expected_uniform_count)

    def barcode_significance(self, barcode_db, min_reads=10):
        """Compute fold enrichment and depletion for each barcode.
        
        Args:
            barcode_db: Dictionary of barcode counts
            min_reads: Minimum absolute read count
            
        Returns:
            dict: Barcode enrichment scores, dict: Barcode depletion scores
        """
        if not barcode_db:
            return {}, {}

        total_counts = sum(barcode_db.values())
        barcodes_num = len(barcode_db)

        if barcodes_num == 0 or total_counts == 0:
            return {}, {}

        expected_uniform_count = total_counts / barcodes_num

        enrichment_scores = {
            barcode: count / expected_uniform_count
            for barcode, count in barcode_db.items()
        }

        # Identify depleted barcodes
        depleted_barcodes = {
            barcode: count for barcode, count in barcode_db.items()
            if count < (0.5 * expected_uniform_count)  # If <50% of expected
        }

        return enrichment_scores, depleted_barcodes


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
                decision_file.write("timestamp\tchannel\tread_id\tbarcode\tmapq\tbaseq\tclassification\n")

            total_reads = sum(self.barcode_db.values())

            # Filter valid target barcodes once
            barcode_db_target = {
                bc: self.barcode_db[bc]
                for bc in self.barcode_db.keys()
                if 'target_barcodes' not in self.seqtagger_params or bc in self.seqtagger_params['target_barcodes']
            }
            total_reads_valid = sum(barcode_db_target.values())

            for barcode, count in self.barcode_db.items():
                if 'target_barcodes' in self.seqtagger_params and barcode in self.seqtagger_params['target_barcodes']:
                    proportion = count / total_reads_valid if total_reads_valid > 0 else 0
                else:
                    proportion = "not a valid barcode"

                stats_file.write(
                    f"{timestamp}\t{total_reads}\t{barcode}\t{count}\t"
                    f"{proportion if isinstance(proportion, str) else f'{proportion:.4f}'}\n"
                )

            # Find barcode with highest frequency
            # top_barcode = max(self.barcode_db, key=self.barcode_db.get) if self.barcode_db else None
            # top_is_sig = self.is_sig(top_barcode, barcode_db_target) if top_barcode else False

            # Compute enrichment & find depleted barcodes
            enrichment_scores, depleted_barcodes = self.barcode_significance(barcode_db_target)

            # Find the most enriched barcode
            top_barcode = max(enrichment_scores, key=enrichment_scores.get) if enrichment_scores else None

            # Threshold for significant enrichment
            threshold = 1.5  
            top_is_sig = enrichment_scores.get(top_barcode, 0) >= threshold if top_barcode else False

            if depleted_barcodes:
                unblocked_barcodes = [bc for bc in enrichment_scores if enrichment_scores[bc] >= 1.0]
            else:
                unblocked_barcodes = [top_barcode] if top_is_sig else []

            for idx, (channel, read) in enumerate(reads):
                barcode = str(calls[idx][0])
                mapq = calls[idx][1]
                baseq = calls[idx][2]

                # Filter calls based on baseq
                if baseq > self.seqtagger_params['min_baseq']:
                    batch_stats['qc_passed'] += 1

                    # Label significantly frequent barcode as target to unblock
                    classification = "TARGET" if barcode in unblocked_barcodes else "NONTARGET"
                    if classification == "TARGET":
                        batch_stats['target_reads'] += 1
                        self.barcode_db_unblock[barcode] += 1
                    else:
                        self.barcode_db[barcode] += 1

                    aln = Aligner(ctg=classification)

                    # Log decision
                    decision_file.write(
                        f"{timestamp}\t{channel}\t{read.id}\t{barcode}\t"
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

    def preprocess_read(self, read: np.ndarray, offset: int) -> np.ndarray:
        """
        Normalize the first `chunksize` values of `read` using median and MAD.
        
        Parameters:
            read (np.ndarray): The raw signal data.
            chunksize (int): The number of data points to process.
        
        Returns:
            np.ndarray: The normalized signal.
        """
        median, mad = get_med_mad(read[:self.demux_chunksize])
        normalized_signal = (read[:self.demux_chunksize_max] - median) / mad
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

        signal_array = np.zeros((self.seqtagger_params['batchsize'], self.demux_chunksize_max), dtype=np.float16)

        # Process all reads in one go
        for idx, (channel, read) in enumerate(reads):
            signal = np.frombuffer(read.raw_data, dtype=signal_dtype)
            signal_pA = (signal + daq_values[channel].offset) * daq_values[channel].scaling
            signal_array[idx] = self.preprocess_read(signal_pA, offset=0)

        # Initial barcode processing
        calls = chunks2barcodes(self.demux_model, signal_array[:, 0:self.demux_chunksize])
        
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
