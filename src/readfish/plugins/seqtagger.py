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
from readfish.plugins.demux import load_demux_model, chunks2barcodes, get_med_mad

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

        self.seqtagger_params['model'] = "/mnt/869990e7-a61f-469f-99fe-a48d24ac44ca/git/crg/b100_RNA004" # b100_RNA004.zip
        self.seqtagger_params['model_secret'] = "iY8Hfro9"
        self.seqtagger_params['batchsize'] = 512
        self.seqtagger_params['min_baseq'] = 50
        self.validate()
        
        self.demux_model = load_demux_model(self.seqtagger_params['model'], self.seqtagger_params['batchsize'], pwd=self.seqtagger_params['model_secret'])
        self.demux_chunksize = self.demux_model.config['basecaller']['chunksize']

        self.barcode_db = defaultdict(dict)

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
        
    def preprocess_read(self, read: np.ndarray, chunksize: int) -> np.ndarray:
        """
        Normalize the first `chunksize` values of `read` using median and MAD.
        
        Parameters:
            read (np.ndarray): The raw signal data.
            chunksize (int): The number of data points to process.
        
        Returns:
            np.ndarray: The normalized signal.
        """
        median, mad = get_med_mad(read[:chunksize])
        normalized_signal = (read[:chunksize] - median) / mad
        return normalized_signal

    def make_decision(self, reads, calls):
        # Count barcode occurrences from the database
        barcode_counts = Counter(self.barcode_db.values())
        most_common = barcode_counts.most_common()

        # Get the top barcode if any exist
        top_barcode = most_common[0][0] if most_common else None

        with open("stat.tsv", "a") as log_file:
            log_file.write(f"{barcode_counts}\n")
            log_file.write(f"{top_barcode}\n")

        # Process each read and call
        with open("debug_log.tsv", "a") as log_file:  # 'a' mode for appending
            for (channel, read), (barcode, mapq, baseq) in zip(reads, calls):
                # Update barcode counts in database for high-quality bases
                if baseq > 50:
                    self.barcode_db[barcode] = self.barcode_db.get(barcode, 0) + 1

                    # Determine alignment type (default is "P")
                    ctg = "U" if (top_barcode is not None and barcode == top_barcode) else "P"
                    aln = Aligner(ctg=ctg)
                    
                    # Prepare data for logging
                    log_entry = f"{channel}\t{read.id}\t{barcode}\t{mapq}\t{baseq}\t{ctg}\n"
                    log_file.write(log_entry)
                    
                    yield (channel, read, barcode, mapq, baseq, aln)

    def call(
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

        batch_size = self.seqtagger_params['batchsize']
        signal_array = np.zeros((batch_size, self.demux_chunksize), dtype=np.float16)

        for channel, read in reads:  # Fixed syntax error (extra closing parenthesis)
            signal_array[channel] = self.preprocess_read(
                np.frombuffer(read.raw_data, dtype=signal_dtype), self.demux_chunksize
            )

        # Process barcodes in one batch  
        calls = chunks2barcodes(self.demux_model, signal_array)

        results = self.make_decision(reads, calls)

        for channel, read, barcode, mapq, baseq, aln in results:  # Fixed syntax and iteration logic
            yield Result(
                channel=channel,
                read_id=read.id,
                seq=None,
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
