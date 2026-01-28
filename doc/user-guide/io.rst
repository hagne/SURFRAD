IO
==

This section describes the file input/output helpers provided by the
``surfradpy.file_io`` package. At the moment, the only concrete reader is for
MFRSR/RSR binary files.

Module overview
---------------

``surfradpy.file_io.mfrsr``
   Reader and parser for RSR binary files (MFRSR/MFR). It unpacks the RSR
   header and records, handles daytime/nighttime channel layouts, and exposes
   a dataset-oriented interface for downstream analysis.

MFRSR RSR reader
----------------
This module mimics the most basic functions of the tu programm, which is reading 
in header and data. It will not apply any correction, this is thesedays done using atmPy

The RSR reader provides a small API to parse raw RSR files and expose the
results as an ``xarray.Dataset`` with time-indexed data and metadata.

Key types and functions
^^^^^^^^^^^^^^^^^^^^^^^

``open_rsr(path)``
   Open an RSR binary file and return an ``RSRFile`` instance. This validates
   the header, determines the instrument type, and prepares the record layout.

``RSRFile``
   Parsed representation of an RSR file. The ``dataset`` property iterates
   through all records, pads missing values with NaNs, and returns an
   ``xarray.Dataset`` with time as the index.

``RSRHeader`` / ``RSRRec``
   Internal dataclasses representing the file header and each data record.

Data outputs
^^^^^^^^^^^^

The ``RSRFile.dataset`` property builds a dataset with:

- ``datetime`` index derived from the record timestamps.
- Variables for irradiance channels, depending on the instrument:
  - MFR: ``alltime``
  - MFRSR: ``alltime``, ``global_horizontal``, ``diffuse_horizontal``,
    ``direct_horizontal``
- Attributes with metadata (instrument, location, logger/head IDs, sampling
  rate, averaging period, and file path).

Errors and data quality
^^^^^^^^^^^^^^^^^^^^^^

The reader raises ``FileCorruptError`` if no valid records can be parsed, and
uses NaNs to represent gaps or missing channels in the record stream. It also
emits warnings for truncated records or malformed record length indicators.

Minimal example
^^^^^^^^^^^^^^^

.. code-block:: python

   from surfradpy.file_io.mfrsr import open_rsr

   rsr = open_rsr("path/to/file.rsr")
   ds = rsr.dataset
   print(ds)
