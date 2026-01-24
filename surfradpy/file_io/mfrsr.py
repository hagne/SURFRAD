from __future__ import annotations


import subprocess
import pandas as _pd
import xarray as _xr
import numpy as np
import io
import atmPy.radiation.retrievals.spectral_irradiance as atmradobs


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

REPORT_END = 0  # default in setup_pp.c
DATE_JOE = "joe"  # only mode we implement

DEFAULT_GAP = -9998
DEFAULT_NIGHTTIME = -9999
DEFAULT_SEP = " "  # setup_pp.c default
DEFAULT_TZ_SECONDS = 0  # setup_pp.c default "timezone=0"


# -----------------------------
# Minimal fmt4/sstoa equivalents
# -----------------------------
def sstoa(v: int) -> str:
    return str(int(v))


def fmt4(x: float) -> str:
    """
    Best-effort stand-in for the C fmt4(): prints a float in a compact way.
    This is only used for non-gap values when separator != NULL (default is " ").
    """
    # Preserve exact integers nicely
    if math.isfinite(x) and abs(x - round(x)) < 1e-12:
        return str(int(round(x)))

    ax = abs(x)
    if ax == 0:
        return "0"

    # Keep ~4 significant digits-ish without being too noisy.
    if ax >= 10000 or ax < 0.001:
        s = f"{x:.4e}"
    elif ax >= 1000:
        s = f"{x:.1f}"
    elif ax >= 100:
        s = f"{x:.2f}"
    elif ax >= 10:
        s = f"{x:.3f}"
    else:
        s = f"{x:.4f}"

    # Trim trailing zeros/dot for readability
    if "e" not in s and "E" not in s:
        s = s.rstrip("0").rstrip(".")
    return s


# -----------------------------
# Time conversion (from rsrlibc.c)
# -----------------------------
DAYS_1900_1970 = 25568
SECONDS_PER_DAY = 86400


# def rsr_unix2j(secs: int) -> float:
#     if secs < 0:
#         return -1.0
#     return float(DAYS_1900_1970) + (float(secs) / SECONDS_PER_DAY)


def rsr_j2unix(jdays: float) -> int:
    if jdays < DAYS_1900_1970:
        return -1
    jdays -= float(DAYS_1900_1970)
    # +0.5 for rounding like C
    secs = int(jdays * 86400.0 + 0.5)
    return secs if secs >= 0 else -1


# -----------------------------
# Solar geometry (ported from sunae.c/h)
# -----------------------------
# @dataclass
# class AEP:
#     az: float = 0.0
#     el: float = 0.0
#     ra: float = 0.0
#     dec: float = 0.0
#     ha: float = 0.0
#     eqt: float = 0.0
#     tst: float = 0.0


# def _mod2pi(x: float) -> float:
#     twopi = 2.0 * math.pi
#     y = x % twopi
#     return y


# def sunae(year: int, doy: int, hour: float, lat_deg: float, lon_deg: float, aep: AEP) -> float:
#     """
#     Port of sunae() from sunae.c (sufficient for correct_direct()).
#     Returns tst (true solar time hours), and fills aep with angles (radians).
#     """
#     rlat = math.radians(lat_deg)
#     rlon = math.radians(lon_deg)

#     # Leap year handling like C
#     leap = (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0))
#     ydays = 366 if leap else 365

#     gamma = (2.0 * math.pi / ydays) * (doy - 1 + (hour - 12.0) / 24.0)

#     eqtime = 229.18 * (
#         0.000075
#         + 0.001868 * math.cos(gamma)
#         - 0.032077 * math.sin(gamma)
#         - 0.014615 * math.cos(2 * gamma)
#         - 0.040849 * math.sin(2 * gamma)
#     )

#     decl = (
#         0.006918
#         - 0.399912 * math.cos(gamma)
#         + 0.070257 * math.sin(gamma)
#         - 0.006758 * math.cos(2 * gamma)
#         + 0.000907 * math.sin(2 * gamma)
#         - 0.002697 * math.cos(3 * gamma)
#         + 0.00148 * math.sin(3 * gamma)
#     )

#     time_offset = eqtime + 4.0 * lon_deg  # minutes (timezone handled by caller; tu uses UTC default)
#     tst = hour * 60.0 + time_offset  # minutes
#     ha = math.radians((tst / 4.0) - 180.0)

#     cos_zen = math.sin(rlat) * math.sin(decl) + math.cos(rlat) * math.cos(decl) * math.cos(ha)
#     cos_zen = max(-1.0, min(1.0, cos_zen))
#     zen = math.acos(cos_zen)
#     el = (math.pi / 2.0) - zen

#     # Azimuth
#     sin_az = -(math.sin(ha) * math.cos(decl)) / max(1e-12, math.cos(el))
#     cos_az = (math.sin(decl) - math.sin(rlat) * math.sin(el)) / max(1e-12, (math.cos(rlat) * math.cos(el)))
#     az = math.atan2(sin_az, cos_az)
#     az = _mod2pi(az)

#     # Fill outputs
#     aep.az = az
#     aep.el = el
#     aep.dec = decl
#     aep.ha = ha
#     aep.eqt = eqtime
#     aep.tst = tst / 60.0  # hours

#     return aep.tst


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


def _rec_lengths(rf: RSRFile) -> None:
    h = rf.head
    if rf.inst_type == Type2:
        c = sum(1 for x in h.counter if x)
        d = sum(1 for x in h.daytime if x)
        n = sum(1 for x in h.all_the_time if x)
        if h.band_on:
            d += 3 * h.diodes

        t = c + n + d
        rf.rec_n_day = t
        rf.rec_len_day = ((t + 3) // 4 + t * 3 + 1) // 2

        t = c + n
        rf.rec_n_nite = t
        rf.rec_len_nite = 0 if t == 0 else ((t + 3) // 4 + t * 3 + 1) // 2

        rf.polar_day = None
        rf.polar_nite = None
    else:
        # nighttime channels
        if h.met_on:
            rf.rec_n_nite = 10 if rf.inst_type == Multi_Filter_32 else 3
        else:
            rf.rec_n_nite = 0
        for i in range(16):
            if h.extra[i]:
                rf.rec_n_nite += 1

        rf.rec_n_day = rf.rec_n_nite
        if h.band_on:
            if rf.inst_type == Single_Channel:
                rf.rec_n_day += 4
            else:
                rf.rec_n_day += 3 * h.diodes

        rf.rec_len_day = int(math.ceil(rf.rec_n_day * 1.5))
        rf.rec_len_nite = int(math.ceil(rf.rec_n_nite * 1.5))

        rf.polar_day = [0] * rf.rec_n_day
        rf.polar_nite = [0] * rf.rec_n_nite if rf.rec_n_nite > 0 else None
        _find_poles(rf)

    rf.record.data = [0] * rf.rec_n_day
    rf.record.n_data = -2
    rf.record.obs_time = rf.start_time


def _find_poles(rf: RSRFile) -> None:
    h = rf.head
    assert rf.polar_day is not None

    # shadowband channels are always unipolar (0)
    d = 0
    if h.band_on:
        if rf.inst_type == Single_Channel:
            for _ in range(4):
                rf.polar_day[d] = 0
                d += 1
        else:
            for _ in range(3 * h.diodes):
                rf.polar_day[d] = 0
                d += 1

    n = 0
    # met channels are unipolar; multi_filter_32 has some bipolar met channels set per C
    if h.met_on:
        if rf.inst_type != Multi_Filter_32:
            for _ in range(3):
                rf.polar_day[d] = 0
                d += 1
                if rf.polar_nite is not None:
                    rf.polar_nite[n] = 0
                n += 1
        else:
            for _ in range(10):
                rf.polar_day[d] = 0
                d += 1
                if rf.polar_nite is not None:
                    rf.polar_nite[n] = 0
                n += 1
            # C sets two of the met channels bipolar (indices relative to end)
            rf.polar_day[d - 5] = 1
            rf.polar_day[d - 8] = 1
            if rf.polar_nite is not None:
                rf.polar_nite[n - 5] = 1
                rf.polar_nite[n - 8] = 1

    # extras inherit bipolar bit
    for i in range(16):
        if h.extra[i]:
            rf.polar_day[d] = h.bipolar[i]
            d += 1
            if rf.polar_nite is not None:
                rf.polar_nite[n] = h.bipolar[i]
            n += 1


def rsr_open_file(path: str) -> RSRFile:
    print(f"Opening RSR file: {path}")
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

    _rec_lengths(rf)
    rf.cdp_pos = 0
    rf.curr_rec_no = 0
    rf._gap_skip = -1
    return rf


def _set_current_time(rf: RSRFile) -> None:
    rf.record.obs_time = rf.start_time + rf.head.avg_period * (rf.curr_rec_no - 1)


def _unpack_1(rf: RSRFile, n: int) -> List[int]:
    p = rf.cdp_pos + 1
    out = [0] * n
    i = 0
    while i < n:
        b0 = rf.data[p]
        b1 = rf.data[p + 1]
        b2 = rf.data[p + 2]
        out[i] = ((b1 >> 4) | (b0 << 4)) & 0xFFFF
        if i + 1 < n:
            out[i + 1] = (((b1 & 0x0F) << 8) | b2) & 0xFFFF
        i += 2
        p += 3

    # apply bipolar conversion
    poles = rf.polar_day if n == rf.rec_n_day else rf.polar_nite
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


def _unpack_2(rf: RSRFile, n: int) -> List[int]:
    # p points to the first byte that contains 12-bit groups
    p = rf.cdp_pos + ((n + 11) // 8)
    nib = ((n - 1) % 8) < 4  # matches C: nib = (n - 1)%8 < 4;
    out = [0] * n
    for i in range(n):
        if nib:
            out[i] = ((rf.data[p] & 0x0F) << 8) | rf.data[p + 1]
            p += 2
        else:
            out[i] = ((rf.data[p + 1] >> 4) | (rf.data[p] << 4)) & 0xFFFF
            p += 1
        nib = not nib

    # apply signs from packed sign bits (start at cdp+1)
    p = rf.cdp_pos + 1
    bm = 0x80
    for i in range(n):
        if rf.data[p] & bm:
            out[i] = -out[i]
        bm >>= 1
        if bm == 0:
            p += 1
            bm = 0x80
    return out


def rsr_next_record(rf: RSRFile) -> RSRRec:
    # End-of-buffer
    if rf.cdp_pos >= rf.data_length:
        rf.record.n_data = -1
        return rf.record

    rf.record.is_gap = 0
    reclen = rf.data[rf.cdp_pos]

    if reclen == rf.rec_len_day or reclen == rf.rec_len_nite:
        rec_n = rf.rec_n_day if reclen == rf.rec_len_day else rf.rec_n_nite
        if rec_n > 0:
            if rf.inst_type == Type2:
                vals = _unpack_2(rf, rec_n)
            else:
                vals = _unpack_1(rf, rec_n)
            # store into record.data (day-sized buffer in C)
            rf.record.data[:rec_n] = vals
        rf.record.n_data = rec_n
        rf.cdp_pos += reclen + 1

    elif reclen == 0xFE:
        # gap record indicator; mirrors static skip in C
        if rf._gap_skip == -1:
            if rf.cdp_pos + 4 >= rf.data_length:
                raise ValueError("Truncated gap record")
            rf._gap_skip = (
                (rf.data[rf.cdp_pos + 1] << 24)
                + (rf.data[rf.cdp_pos + 2] << 16)
                + (rf.data[rf.cdp_pos + 3] << 8)
                + rf.data[rf.cdp_pos + 4]
            )
            if rf._gap_skip <= 0:
                raise ValueError("Gap record skip count must be > 0")
            for i in range(rf.rec_n_day):
                rf.record.data[i] = -9999

        if rf._gap_skip == 1:
            rf.cdp_pos += 5
            rf._gap_skip = -1
        else:
            rf._gap_skip -= 1

        rf.record.is_gap = 1
        rf.record.n_data = rf.rec_n_day

    elif reclen == 0xFF:
        rf.record.n_data = -1
        return rf.record

    else:
        raise ValueError(f"Bad record length indicator {reclen} at byte {rf.cdp_pos + rf.header_length}")

    rf.curr_rec_no += 1
    _set_current_time(rf)
    return rf.record


# -----------------------------
# unpack.c logic needed for your invocation
# -----------------------------
# def _call_sunae(inp: RSRFile, report_time: int, midpt: bool) -> tuple[AEP, float]:
#     aep = AEP()
#     t = inp.record.obs_time + 5  # C: obs_time + 5 seconds
#     if midpt and inp.head.sample_rate != inp.head.avg_period:
#         t -= inp.head.sample_rate / 2.0
#     if report_time != REPORT_END:
#         # only END is used in your invocation; included for completeness
#         if report_time == 1:  # Mid
#             t -= inp.head.avg_period / 2.0
#         elif report_time == 2:  # Start
#             t -= float(inp.head.avg_period)

#     # gmtime
#     dt = datetime.utcfromtimestamp(int(t))
#     year = dt.year
#     doy = int(dt.strftime("%j"))
#     hour = dt.hour + dt.minute / 60.0 + dt.second / 3600.0

#     # C uses lon = -longitude
#     lon = -1.0 * inp.head.longitude
#     lat = inp.head.latitude

#     tst = sunae(year, doy, hour, lat, lon, aep)
#     return aep, tst


# def correct_direct(inp: RSRFile, report_time: int, vec: List[float]) -> None:
#     """
#     Port of correct_direct() from unpack.c, acting in-place on vec.
#     """
#     # only when sample_rate == avg_period, non-gap, not single-channel, band_on, diodes != 0
#     if inp.head.sample_rate != inp.head.avg_period:
#         return
#     if inp.record.n_data == inp.rec_n_nite:
#         return
#     if inp.inst_type == Single_Channel:
#         return
#     if not inp.head.band_on:
#         return
#     if inp.head.diodes == 0:
#         return

#     aep, _tst = _call_sunae(inp, report_time=report_time, midpt=False)
#     if aep.el < (10.0 * math.pi / 180.0):
#         di_idx = inp.head.diodes * 2
#         cos_el = math.cos((math.pi / 2.0) - aep.el)
#         if cos_el > 0:
#             vec[di_idx] = vec[di_idx] / cos_el


def emit_header(inp: RSRFile) -> str:
    h = inp.head
    if inp.inst_type == Type2:
        daytime = "".join(str(int(x)) for x in h.daytime)
        allthetime = "".join(str(int(x)) for x in h.all_the_time)
        counter = "".join(str(int(x)) for x in h.counter)
        return (
            f"FORMAT VERSION={h.soft_rev} UNIT_ID={h.unit_id} HEAD_ID={h.head_id} "
            f"LONGITUDE={h.longitude:.4f} LATITUDE={h.latitude:.4f} FLAGS=0x{h.flags:x} "
            f"AVERAGING PERIOD={h.avg_period} SAMPLE RATE={h.sample_rate} "
            f"DIODES={h.diodes} DAYTIME=0x{h.raw_daytime:x} ALLTHETIME=0x{h.raw_all_the_time:x} "
            f"COUNTERS=0x{h.raw_counter:x} DAYTIME_FLAGS={daytime} "
            f"ALLTHETIME_FLAGS={allthetime} COUNTER_FLAGS={counter} BYTES={inp.data_length}"
        )
    else:
        return (
            f"FORMAT VERSION={h.soft_rev} UNIT_ID={h.unit_id} LONGITUDE={h.longitude:.4f} "
            f"LATITUDE={h.latitude:.4f} FLAGS=0x{h.flags:x} AVERAGING PERIOD={h.avg_period} "
            f"SAMPLE RATE={h.sample_rate} GAIN={h.rsr_gain} OFFSET={h.rsr_offset} DIODES={h.diodes} "
            f"EXTRA=0x{h.raw_extra:x} BIPOLAR=0x{h.raw_bipolar:x} BYTES={inp.data_length}"
        )


# def emit_date_time_joe(inp: RSRFile, tz_seconds: int = 0) -> str:
#     # C: t = obs_time; adjust by report_time (END => none), then add timezone
#     t = inp.record.obs_time + tz_seconds
#     return f"{rsr_unix2j(int(t)):.5f}"


def emit_data(
                inp: RSRFile,
                vec: Optional[List[float]],
                sep: str = DEFAULT_SEP,
                nighttime_placeholder: Optional[int] = DEFAULT_NIGHTTIME,
                gap_placeholder: int = DEFAULT_GAP,
                out = None,
            ) -> str:
    if isinstance(out, type(None)):
        out = {}

    h = inp.head
    parts: List[str] = []
    partsf: List[str] = []
    if inp.record.is_gap:
        cvt = sstoa(gap_placeholder)
    else:
        cvt = ""

    # pre-type2: if nighttime record, prepend placeholders for missing daytime-only channels
    if inp.inst_type != Type2 and nighttime_placeholder is not None:
        if inp.record.n_data == inp.rec_n_nite:
            missing = inp.rec_n_day - inp.record.n_data
            for _ in range(missing):
                parts.append(sstoa(nighttime_placeholder))

    # actual values
    n = inp.record.n_data
    for i in range(n):
        if inp.record.is_gap:
            parts.append(cvt)
            partsf.append(np.nan)
        else:
            assert vec is not None
            parts.append(fmt4(vec[i]))
            out['veci'] = vec[i]
            partsf.append(vec[i])

    # type2: pad remaining with nighttime placeholders (if requested)
    if inp.inst_type == Type2 and nighttime_placeholder is not None:
        for _ in range(inp.rec_n_day - inp.record.n_data):
            parts.append(sstoa(nighttime_placeholder))
            partsf.append(np.nan)
    out['parts_str'] = sep + sep.join(parts) if parts else ""
    out['parts_f'] = partsf
    return out


def tu_like_unpack(path: str, out = None) -> None:

    if isinstance(out, type(None)):
        out = {}

    inp = rsr_open_file(path)
    out['inp'] = inp

    # loop records
    i = 0
    data = []
    # while i < 100:
    while True:
        i += 1
        r = rsr_next_record(inp)
        if r.n_data < 0:
            break

        if r.is_gap:
            vec = None
        else:
            vec = [float(v) for v in r.data[:r.n_data]]
            # correct_direct(inp, report_time=REPORT_END, vec=vec)
        emit_data_out = emit_data(inp, vec, 
                                  sep=DEFAULT_SEP, 
                                  nighttime_placeholder=DEFAULT_NIGHTTIME, 
                                  gap_placeholder=DEFAULT_GAP,
                                #   out = out
                                  )
        # line += emit_data_out['parts_str']
        # lines.append(line)
        row = [inp.record.obs_time] +  emit_data_out['parts_f']
        data.append(row)
    # out['lines'] = lines
    out['data'] = data
    # out['data_str'] = '\n'.join(lines)
    return out



# if __name__ == "__main__":
#     main()
# python tu_py.py -d joe -H /path/to/file > out.txt


def read_raw(path2file, out = None, use_str = False):
    """
    Opens a MFRSR or MFR raw file. Those are the files that have an extension 
    like .mtm, .xmd, .rsr. Surfrad files are stored here:
    /nfs/grad/Inst/MFR/SURFRAD/{site}/mfrsr/raw/
    This requires the "tu" program to be installed on the system. You can get 
    it at: https://github.com/HagenTelg/tu_mfrsr_reader

    Parameters
    ----------
    path2file : TYPE
        DESCRIPTION.
    read_header_only : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    if isinstance(out, type(None)):
        out = {}


    out = tu_like_unpack(path2file, 
                        #   True, 'joe', 
                         out = out 
                         )
    
    if use_str:
        print('=====================')
        print('++++ Don"t use this path anymore, use use_str = False ++++')
        outstr = out['data_str']
        # out['rawstring'] = outstr
        # print(f'error: {err}')


        #### the header
        header = outstr.split('\n')[0]    
        # if read_header_only:
        #     return header
        
        bla = header.split()[5]
        if bla == '0':
            instrument = 'mfr'
        else:
            instrument = 'mfrsr'
        
        
        #### format the data
        df = _pd.read_csv(io.StringIO(outstr), sep = r'\s+', header = None, skiprows=1)
        df.index = df.apply(lambda row: _pd.to_datetime('19000101') + _pd.to_timedelta(row[0] - 1, 'd'), axis = 1)
        df.index.name = 'datetime'
        df = df.where(df!=-9999)
        
        # photodiode (channel) collumns
        di = 7
        si = 2
        alltime = df.iloc[:,si: si+di]
        out['df'] = df.copy()
        if instrument == 'mfr':
            alltime.columns = range(7)
            alltime.columns.name = 'channel'
            
            ds = _xr.Dataset()
            # ds['sun_position'] = sun_position
            ds['alltime'] = alltime
            ds.attrs = attrs
            # obs = atmradobs.GlobalHorizontalIrradiation(ds)
        
        elif instrument == 'mfrsr':
            si = 11
            global_horizontal = df.iloc[:,si: si+di]
            si =  18
            diffuse_horizontal = df.iloc[:,si: si+di]
            si =  25
            direct = df.iloc[:,si: si+di]
            
            for dft in [alltime, global_horizontal, diffuse_horizontal, direct]:
                out['dft'] = dft
                dft.columns = range(7)
                dft.columns.name = 'channel'
            
            ds = _xr.Dataset()
            # ds['sun_position'] = sun_position
            ds['alltime'] = alltime
            ds['global_horizontal'] = global_horizontal
            ds['diffuse_horizontal'] = diffuse_horizontal
            ds['direct_horizontal'] = direct
            # ds.attrs = attrs
            
            # obs = atmradobs.CombinedGlobalDiffuseDirect(ds, 
            #                                             # raise_error_if_no_channel_wavelength = False
            #               
        else:
            assert(False), 'nop, not possible!' 
    else:
        data = np.array(out['data'])
        inp = out['inp']
        df = _pd.DataFrame(data, index = _pd.to_datetime(data[:,0], unit='s'))
        df.index.name = 'datetime'

        # out['rawstring'] = outstr
        # print(f'error: {err}')


        #### the header
        # header = outstr.split('\n')[0]    
        # if read_header_only:
        #     return header
        # return
        # bla = header.split()[5]
        if inp.head.band_on == 0:
            instrument = 'mfr'
        elif inp.head.band_on == 1:
            instrument = 'mfrsr'
        else:
            assert(False), f'nop, not possible! '
        
        
        # #### format the data
        # df = _pd.read_csv(io.StringIO(outstr), sep = r'\s+', header = None, skiprows=1)
        # df.index = df.apply(lambda row: _pd.to_datetime('19000101') + _pd.to_timedelta(row[0] - 1, 'd'), axis = 1)
        # df.index.name = 'datetime'
        # df = df.where(df!=-9999)


        # if test == 1:
        #     return df, header
        #### assign collumns
        # attrs = dict(site_longitude = - float(header.split()[4]),
        #              site_latitude = float(header.split()[3]),
        #              site_elevation = 0,
        #              site = 'TMP',
        #              site_name = 'unknown',
        #              calibrated_irradiance = 'False',
        #              calibrated_cosine = 'False',
        #              info = 'This is the raw file',
        #              file_type = int(header.split()[0]),
        #              serial_no = header.split()[1],
        #              path2file = path2file,
        #              measurement_sequenc = ', '.join(header.split()[5:8]),
        #              instrument_type = instrument,
        #             )
        # non photodiod collums
        # 0: time
        # 1: no idea
        # 9: no idea
        # 10: no idea
        
        # photodiode (channel) collumns
        di = 7
        si = 2
        alltime = df.iloc[:,si: si+di]
        
        out['df'] = df.copy()
        if instrument == 'mfr':
            alltime.columns = range(7)
            alltime.columns.name = 'channel'
            
            ds = _xr.Dataset()
            # ds['sun_position'] = sun_position
            ds['alltime'] = alltime
            # ds.attrs = attrs
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
            # ds.attrs = attrs
            
            # obs = atmradobs.CombinedGlobalDiffuseDirect(ds, 
            #                                             # raise_error_if_no_channel_wavelength = False
            #                                             )
    

    out['ds'] = ds
    return out