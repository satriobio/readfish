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

from demux import load_demux_model, chunks2barcodes, get_med_mad

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

    def disconnect(self) -> None:
        pass

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

        for idx, (channel, read) in enumerate(reads):
            signal = np.frombuffer(read.raw_data, dtype=signal_dtype)
            signal_pA = (signal + daq_values[channel].offset) * daq_values[channel].scaling
            signal_array[idx] = self.preprocess_read(signal_pA, offset=0)

        tmp_calls = []
        final_call = []

        for s in [0, 1000, 2000]:
            raw_calls = chunks2barcodes(self.demux_model, signal_array[:, s:s+self.demux_chunksize])
            tmp_calls.append(raw_calls[:len(reads)])

        calls_by_read = list(zip(*tmp_calls))

        # Select the best call for each read (highest baseq)
        for read_calls in calls_by_read:
            best_call = max(read_calls, key=lambda x: x[2])
            final_call.append(best_call)
        
        # results = self.make_decision(reads, final_call)

        # for channel, read, barcode, mapq, baseq, aln in results:
        #     yield Result(
        #         channel=channel,
        #         read_id=read.id,
        #         seq=[],
        #         barcode=None,
        #         basecall_data=(barcode, mapq, baseq),
        #         alignment_data=[aln],
        #     )

        # for idx, (channel, read, barcode, mapq, baseq, aln in enumerate(reads):
        for idx, (channel, read) in enumerate(reads):
            yield Result(
                channel=channel,
                read_id=read.id,
                seq=[],
                barcode=None,
                basecall_data=final_call[idx],
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
