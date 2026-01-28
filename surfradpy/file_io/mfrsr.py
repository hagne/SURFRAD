from __future__ import annotations


import subprocess
import pandas as _pd
import xarray as _xr
import numpy as np
import io
import atmPy.radiation.retrievals.spectral_irradiance as atmradobs
import warnings
import pathlib as pl

import math
from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Optional



class FileCorruptError(Exception):
    """Custom exception raised when a file is detected as corrupt."""
    def __init__(self, filename, message="File is corrupt or unreadable."):
        self.filename = filename
        self.message = f"{message} ({filename})"
        super().__init__(self.message)
# -----------------------------
# Constants (from rsrlib.h / tu.h / setup_pp defaults)
# -----------------------------
TYPE1_HEADER_SIZE = 27
TYPE2_HEADER_SIZE = 35

Unknown_Instrument = 0
Multi_Filter_16 = 1
Multi_Filter_32 = 2
Type2 = 3
Single_Channel = 4


# -----------------------------
# Time conversion (from rsrlibc.c)
# -----------------------------
DAYS_1900_1970 = 25568
SECONDS_PER_DAY = 86400


def rsr_j2unix(jdays: float) -> int:
    if jdays < DAYS_1900_1970:
        return -1
    jdays -= float(DAYS_1900_1970)
    # +0.5 for rounding like C
    secs = int(jdays * 86400.0 + 0.5)
    return secs if secs >= 0 else -1


# -----------------------------
# RSR structures / parsing (ported from rsrlibc.c/h)
# -----------------------------
@dataclass
class RSRHeader:
    soft_rev: int = 0
    logger_id: str = ""
    head_id: str = ""  # type2
    longitude: float = 0.0
    latitude: float = 0.0
    flags: int = 0
    immed_out: int = 0
    low_power: int = 0
    band_on: int = 0
    met_on: int = 0
    dummy: int = 0
    volt_dog: int = 0
    oao: int = 0
    halt: int = 0

    raw_extra: int = 0
    raw_bipolar: int = 0
    extra: List[int] = field(default_factory=lambda: [0] * 16)
    bipolar: List[int] = field(default_factory=lambda: [0] * 16)

    sample_rate: int = 0
    avg_period: int = 0
    secs_today: int = 0
    days_1900: int = 0

    rsr_gain: int = 0
    rsr_offset: int = 0  # single-channel only
    diodes: int = 0

    # type2 channel masks
    daytime: List[int] = field(default_factory=lambda: [0] * 32)
    all_the_time: List[int] = field(default_factory=lambda: [0] * 32)
    counter: List[int] = field(default_factory=lambda: [0] * 6)
    raw_daytime: int = 0
    raw_all_the_time: int = 0
    raw_counter: int = 0
    err: int = 0


@dataclass
class RSRRec:
    obs_time: int = 0
    is_gap: int = 0
    n_data: int = -2
    data: List[int] = field(default_factory=list)


@dataclass
class RSRFile:
    filename: str
    inst_type: int
    header_length: int
    head: RSRHeader
    data: bytes
    data_length: int

    start_time: int
    curr_rec_no: int = 0
    cdp_pos: int = 0  # index into data (points at length byte)
    rec_n_day: int = 0
    rec_n_nite: int = 0
    rec_len_day: int = 0
    rec_len_nite: int = 0
    polar_day: Optional[List[int]] = None
    polar_nite: Optional[List[int]] = None
    record: RSRRec = field(default_factory=RSRRec)

    _gap_skip: int = -1  # mirrors static skip in rsr_next_record

    def __post_init__(self):
        self._dataset = None

    @property
    def dataset(self) -> _xr.Dataset:
        if isinstance(self._dataset, type(None)):
            # loop records to get the data
            i = 0
            data = []
            # while i < 100:
            while True:
                i += 1
                # r = rsr_next_record(inp)
                r = self.next_record
                if r.n_data < 0:
                    break

                if r.is_gap:
                    vec = None
                else:
                    vec = [float(v) for v in r.data[:r.n_data]]
                # print('doing it')
                emit_data_out = self.pad_data(vec, 
                                        # sep=DEFAULT_SEP, 
                                        # nighttime_placeholder=DEFAULT_NIGHTTIME, 
                                        # gap_placeholder=DEFAULT_GAP,
                                        #   out = out
                                        )

                row = [self.record.obs_time] +  emit_data_out['parts_f']
                data.append(row)

            data = np.array(data)
            self.tp_data = data
            if data.shape[0] == 0:
                raise FileCorruptError(self.filename, message="No valid records found in file")
            
            df = _pd.DataFrame(data, index = _pd.to_datetime(data[:,0], unit='s'))
            df.index.name = 'datetime'

            if self.head.band_on == 0:
                instrument = 'mfr'
            elif self.head.band_on == 1:
                instrument = 'mfrsr'
            else:
                assert(False), f'nop, not possible! '
        
            di = 7
            si = 2
            alltime = df.iloc[:,si: si+di]
            
            # out['df'] = df.copy()
            if instrument == 'mfr':
                alltime.columns = range(7)
                alltime.columns.name = 'channel'
                
                ds = _xr.Dataset()
                ds['alltime'] = alltime
                # obs = atmradobs.GlobalHorizontalIrradiation(ds)
            
            elif instrument == 'mfrsr':
                si = 11
                global_horizontal = df.iloc[:,si: si+di]
                si =  18
                diffuse_horizontal = df.iloc[:,si: si+di]
                si =  25
                direct = df.iloc[:,si: si+di]
                
                for dft in [alltime, global_horizontal, diffuse_horizontal, direct]:
                    dft.columns = range(7)
                    dft.columns.name = 'channel'
                
                ds = _xr.Dataset()
                # ds['sun_position'] = sun_position
                ds['alltime'] = alltime
                ds['global_horizontal'] = global_horizontal
                ds['diffuse_horizontal'] = diffuse_horizontal
                ds['direct_horizontal'] = direct

            # add metadata
            ds.attrs['instrument'] = instrument
            ds.attrs['start_time'] = _pd.to_datetime(self.start_time, unit='s').__str__()
            ds.attrs['avg_period'] = self.head.avg_period
            ds.attrs['sample_rate'] = self.head.sample_rate
            ds.attrs['logger_id'] = self.head.logger_id
            ds.attrs['head_id'] = self.head.head_id if self.inst_type == Type2 else 'N/A'
            ds.attrs['soft_rev'] = self.head.soft_rev
            ds.attrs['latitude'] = self.head.latitude
            ds.attrs['longitude'] = self.head.longitude
            ds.attrs['band_on'] = self.head.band_on
            ds.attrs['diodes'] = self.head.diodes
            ds.attrs['path2file'] = pl.Path(self.filename).as_posix()
            self._dataset = ds
            
        return self._dataset

    @property
    def next_record(self) -> RSRRec:
        # End-of-buffer
        if self.cdp_pos >= self.data_length:
            self.record.n_data = -1
            return self.record
        self.record.is_gap = 0
        reclen = self.data[self.cdp_pos]

        if reclen == self.rec_len_day or reclen == self.rec_len_nite:
            # Ensure record payload fits within data buffer
            end = self.cdp_pos + reclen
            if end >= self.data_length:
                warnings.warn(
                    f"Truncated record: reclen={reclen} at data index {self.cdp_pos}, "
                    f"data_length={self.data_length}"
                )
            rec_n = self.rec_n_day if reclen == self.rec_len_day else self.rec_n_nite
            if rec_n > 0:
                if self.inst_type == Type2:
                    vals = self._unpack_2(rec_n)
                else:
                    vals = self._unpack_1(rec_n)
                # store into record.data (day-sized buffer in C)
                self.record.data[:rec_n] = vals
            self.record.n_data = rec_n
            self.cdp_pos += reclen + 1

        elif reclen == 0xFE:
            # gap record indicator; mirrors static skip in C
            if self._gap_skip == -1:
                if self.cdp_pos + 4 >= self.data_length:
                    raise ValueError("Truncated gap record")
                self._gap_skip = (
                    (self.data[self.cdp_pos + 1] << 24)
                    + (self.data[self.cdp_pos + 2] << 16)
                    + (self.data[self.cdp_pos + 3] << 8)
                    + self.data[self.cdp_pos + 4]
                )
                if self._gap_skip <= 0:
                    raise ValueError("Gap record skip count must be > 0")
                for i in range(self.rec_n_day):
                    self.record.data[i] = -9999

            if self._gap_skip == 1:
                self.cdp_pos += 5
                self._gap_skip = -1
            else:
                self._gap_skip -= 1

            self.record.is_gap = 1
            self.record.n_data = self.rec_n_day

        elif reclen == 0xFF:
            self.record.n_data = -1
            return self.record

        else:
            warnings.warn(
                f"Bad record length indicator {reclen} at byte {self.cdp_pos + self.header_length}; "
                "treating as end-of-data"
            )
            self.record.n_data = -1
            self.cdp_pos = self.data_length
            return self.record

        self.curr_rec_no += 1
        # set current time
        self.record.obs_time = self.start_time + self.head.avg_period * (self.curr_rec_no - 1)
        return self.record

    def pad_data(self: RSRFile,
                    vec: Optional[List[float]],
                    # sep: str = DEFAULT_SEP,
                    # nighttime_placeholder: Optional[int] = DEFAULT_NIGHTTIME,
                    # gap_placeholder: int = DEFAULT_GAP,
                    out = None,
                ) -> dict:
        if isinstance(out, type(None)):
            out = {}

        partsf: List[str] = []

        # pre-type2: if nighttime record, prepend placeholders for missing daytime-only channels
        if self.inst_type != Type2:# and nighttime_placeholder is not None:
            if self.record.n_data == self.rec_n_nite:
                missing = self.rec_n_day - self.record.n_data
                for _ in range(missing):
                    partsf.append(np.nan)

        # actual values
        n = self.record.n_data
        for i in range(n):
            if self.record.is_gap:
                partsf.append(np.nan)
            else:
                assert vec is not None
                partsf.append(vec[i])

        # type2: pad remaining with nighttime placeholders (if requested)
        if self.inst_type == Type2:
            for _ in range(self.rec_n_day - self.record.n_data):
                partsf.append(np.nan)
        out['parts_f'] = partsf
        return out

    def _unpack_2(self: RSRFile, n: int) -> List[int]:
        # p points to the first byte that contains 12-bit groups
        p = self.cdp_pos + ((n + 11) // 8)
        nib = ((n - 1) % 8) < 4  # matches C: nib = (n - 1)%8 < 4;
        out = [0] * n
        n_read = 0
        for i in range(n):
            try:
                if nib:
                    out[i] = ((self.data[p] & 0x0F) << 8) | self.data[p + 1]
                    p += 2
                else:
                    out[i] = ((self.data[p + 1] >> 4) | (self.data[p] << 4)) & 0xFFFF
                    p += 1
                nib = not nib
                n_read += 1
            except IndexError:
                warnings.warn(
                    f"Truncated record while unpacking 12-bit values at data index {p}, data_length={self.data_length}")
                break # TODO: could more data eventually be salvaged?

        # apply signs from packed sign bits (start at cdp+1)
        p = self.cdp_pos + 1
        bm = 0x80
        for i in range(n_read): # only loop over number of n that were actually read
            if self.data[p] & bm:
                out[i] = -out[i]
            bm >>= 1
            if bm == 0:
                p += 1
                bm = 0x80
        return out

    def _unpack_1(self, n: int) -> List[int]:
        p = self.cdp_pos + 1
        out = [0] * n
        i = 0
        while i < n:
            b0 = self.data[p]
            b1 = self.data[p + 1]
            b2 = self.data[p + 2]
            out[i] = ((b1 >> 4) | (b0 << 4)) & 0xFFFF
            if i + 1 < n:
                out[i + 1] = (((b1 & 0x0F) << 8) | b2) & 0xFFFF
            i += 2
            p += 3

        # apply bipolar conversion
        poles = self.polar_day if n == self.rec_n_day else self.polar_nite
        if poles is None:
            return out

        for i in range(n):
            if poles[i]:
                v = out[i]
                if v >= 2048:
                    out[i] = 2 * (v - 4096)
                else:
                    out[i] = 2 * v
        return out
    
    def _rec_lengths(self: RSRFile) -> None:
        h = self.head
        if self.inst_type == Type2:
            c = sum(1 for x in h.counter if x)
            d = sum(1 for x in h.daytime if x)
            n = sum(1 for x in h.all_the_time if x)
            if h.band_on:
                d += 3 * h.diodes

            t = c + n + d
            self.rec_n_day = t
            self.rec_len_day = ((t + 3) // 4 + t * 3 + 1) // 2
            t = c + n
            self.rec_n_nite = t
            self.rec_len_nite = 0 if t == 0 else ((t + 3) // 4 + t * 3 + 1) // 2

            self.polar_day = None
            self.polar_nite = None
        else:
            # nighttime channels
            if h.met_on:
                self.rec_n_nite = 10 if self.inst_type == Multi_Filter_32 else 3
            else:
                self.rec_n_nite = 0
            for i in range(16):
                if h.extra[i]:
                    self.rec_n_nite += 1

            self.rec_n_day = self.rec_n_nite
            if h.band_on:
                if self.inst_type == Single_Channel:
                    self.rec_n_day += 4
                else:
                    self.rec_n_day += 3 * h.diodes

            self.rec_len_day = int(math.ceil(self.rec_n_day * 1.5))
            self.rec_len_nite = int(math.ceil(self.rec_n_nite * 1.5))

            self.polar_day = [0] * self.rec_n_day
            self.polar_nite = [0] * self.rec_n_nite if self.rec_n_nite > 0 else None
            _find_poles(self)
        self.record.data = [0] * self.rec_n_day
        self.record.n_data = -2
        self.record.obs_time = self.start_time
        return

    def _find_poles(self: RSRFile) -> None:
        h = self.head
        assert self.polar_day is not None

        # shadowband channels are always unipolar (0)
        d = 0
        if h.band_on:
            if self.inst_type == Single_Channel:
                for _ in range(4):
                    self.polar_day[d] = 0
                    d += 1
            else:
                for _ in range(3 * h.diodes):
                    self.polar_day[d] = 0
                    d += 1

        n = 0
        # met channels are unipolar; multi_filter_32 has some bipolar met channels set per C
        if h.met_on:
            if self.inst_type != Multi_Filter_32:
                for _ in range(3):
                    self.polar_day[d] = 0
                    d += 1
                    if self.polar_nite is not None:
                        self.polar_nite[n] = 0
                    n += 1
            else:
                for _ in range(10):
                    self.polar_day[d] = 0
                    d += 1
                    if self.polar_nite is not None:
                        self.polar_nite[n] = 0
                    n += 1
                # C sets two of the met channels bipolar (indices relative to end)
                self.polar_day[d - 5] = 1
                self.polar_day[d - 8] = 1
                if self.polar_nite is not None:
                    self.polar_nite[n - 5] = 1
                    self.polar_nite[n - 8] = 1

        # extras inherit bipolar bit
        for i in range(16):
            if h.extra[i]:
                self.polar_day[d] = h.bipolar[i]
                d += 1
                if self.polar_nite is not None:
                    self.polar_nite[n] = h.bipolar[i]
                n += 1


def _unhead_1(inst_type: int, buf: bytes) -> tuple[RSRHeader, int]:
    h = RSRHeader()
    h.soft_rev = buf[0]
    h.logger_id = f'{(buf[1] << 8) + buf[2]:X}'
    h.longitude = 360.0 * ((buf[3] << 8) + buf[4]) / 65536.0
    h.latitude = 360.0 * ((buf[5] << 8) + buf[6]) / 65536.0

    h.flags = buf[7]
    h.immed_out = 1 if (h.flags & 0x80) else 0
    h.low_power = 1 if (h.flags & 0x40) else 0
    h.band_on = 1 if (h.flags & 0x20) else 0
    h.met_on = 1 if (h.flags & 0x10) else 0
    h.dummy = 1 if (h.flags & 0x08) else 0
    h.volt_dog = 1 if (h.flags & 0x04) else 0
    h.oao = 1 if (h.flags & 0x02) else 0
    h.halt = 1 if (h.flags & 0x01) else 0

    h.raw_extra = (buf[8] << 8) + buf[9]
    h.raw_bipolar = (buf[10] << 8) + buf[11]

    # extra bits
    i = 0
    bp = 0x80
    while bp:
        h.extra[i] = 1 if (buf[8] & bp) else 0
        i += 1
        bp >>= 1
    bp = 0x80
    while bp:
        h.extra[i] = 1 if (buf[9] & bp) else 0
        i += 1
        bp >>= 1

    # bipolar bits
    i = 0
    bp = 0x80
    while bp:
        h.bipolar[i] = 1 if (buf[10] & bp) else 0
        i += 1
        bp >>= 1
    bp = 0x80
    while bp:
        h.bipolar[i] = 1 if (buf[11] & bp) else 0
        i += 1
        bp >>= 1

    h.sample_rate = (buf[12] << 8) + buf[13]
    h.avg_period = (buf[14] << 8) + buf[15]
    h.secs_today = (buf[16] << 16) + (buf[17] << 8) + buf[18]
    h.days_1900 = (buf[19] << 8) + buf[20]
    h.rsr_gain = (buf[21] << 8) + buf[22]

    if inst_type == Single_Channel:
        h.rsr_offset = (buf[23] << 8) + buf[24]
        h.diodes = 1
    else:
        h.diodes = buf[23]
        h.rsr_offset = 0

    data_length = (buf[25] << 8) + buf[26]
    return h, data_length


def _unhead_2(buf: bytes) -> tuple[RSRHeader, int]:
    h = RSRHeader()
    h.soft_rev = buf[0]
    h.logger_id = f'{(buf[1] << 8) + buf[2]:X}'
    h.head_id = f'{(buf[3] << 8) + buf[4]:X}'
    h.longitude = 360.0 * ((buf[5] << 8) + buf[6]) / 65536.0
    h.latitude = 360.0 * ((buf[7] << 8) + buf[8]) / 65536.0

    h.flags = buf[9]
    h.immed_out = 1 if (h.flags & 0x80) else 0
    h.low_power = 1 if (h.flags & 0x40) else 0
    h.band_on = 1 if (h.flags & 0x20) else 0
    h.met_on = 1 if (h.flags & 0x10) else 0
    h.dummy = 1 if (h.flags & 0x08) else 0
    h.volt_dog = 1 if (h.flags & 0x04) else 0
    h.oao = 1 if (h.flags & 0x02) else 0
    h.halt = 1 if (h.flags & 0x01) else 0

    h.sample_rate = (buf[10] << 8) + buf[11]
    h.avg_period = (buf[12] << 8) + buf[13]
    h.days_1900 = (buf[14] << 8) + buf[15]
    h.secs_today = (buf[16] << 16) + (buf[17] << 8) + buf[18]
    h.diodes = buf[19]

    tul = (buf[20] << 24) + (buf[21] << 16) + (buf[22] << 8) + buf[23]
    h.raw_daytime = tul
    bp = 0x80000000
    i = 31
    while bp:
        h.daytime[i] = 1 if (tul & bp) else 0
        i -= 1
        bp >>= 1

    tul = (buf[24] << 24) + (buf[25] << 16) + (buf[26] << 8) + buf[27]
    h.raw_all_the_time = tul
    bp = 0x80000000
    i = 31
    while bp:
        h.all_the_time[i] = 1 if (tul & bp) else 0
        i -= 1
        bp >>= 1

    tul = (buf[28] << 16) + (buf[29] << 8) + buf[30]
    h.raw_counter = tul
    bp = 0x00F00000
    i = 5
    while bp:
        h.counter[i] = 1 if (tul & bp) else 0
        i -= 1
        bp >>= 4

    data_length = (buf[31] << 16) + (buf[32] << 8) + buf[33]
    h.err = buf[34]
    return h, data_length


def open_rsr(path: str) -> RSRFile:
    # print(f"Opening RSR file: {path}")
    with open(path, "rb") as f:
        raw = f.read()

    if len(raw) < 1:
        raise ValueError("File is empty")

    c = raw[0]
    if c in (10, 11):
        inst_type = Multi_Filter_16
        header_length = TYPE1_HEADER_SIZE
    elif c == 12:
        inst_type = Multi_Filter_32
        header_length = TYPE1_HEADER_SIZE
    elif c == 13:
        inst_type = Type2
        header_length = TYPE2_HEADER_SIZE
    elif c in (77, 78):
        inst_type = Single_Channel
        header_length = TYPE1_HEADER_SIZE
    else:
        raise ValueError(f"Not an RSR file: byte 1 is {c}")

    if len(raw) < header_length:
        raise ValueError("File too short to contain a full header")

    hb = raw[:header_length]
    if inst_type == Type2:
        head, data_length = _unhead_2(hb)
    else:
        head, data_length = _unhead_1(inst_type, hb)

    # print(head)

    # start_time: rsr_j2unix(days_1900 + secs_today/86400.0)
    start_time = rsr_j2unix(head.days_1900 + head.secs_today / 86400.0)
    if start_time < 0:
        raise ValueError("Suspicious header time; cannot compute start_time")

    data = raw[header_length:]
    if data_length == 0:
        data_length = len(data)
    else:
        data_length = min(data_length, len(data))

    rf = RSRFile(
        filename=path,
        inst_type=inst_type,
        header_length=header_length,
        head=head,
        data=data[:data_length],
        data_length=data_length,
        start_time=start_time,
    )

    rf._rec_lengths()
    rf.cdp_pos = 0
    rf.curr_rec_no = 0
    rf._gap_skip = -1
    return rf


