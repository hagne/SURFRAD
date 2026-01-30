Command-line tools
==================

This section lists the console entry points installed with ``surfradpy``.
Each command maps to a function in ``surfradpy.scripts``.

Available commands
------------------

- ``produce_surfrad_mfrsr2nc`` → ``surfradpy.scripts.surfrad_mfrsr_raw2netcdf:run``
- ``produce_aod2netcdf`` → ``surfradpy.scripts.aod2netcdf:run``
- ``qcrad2ncei`` → ``surfradpy.scripts.qcrad2ncei:main``
- ``aodinv1`` → ``surfradpy.scripts.aodinversion:produce_aodinversion1_01``
- ``aodinv1cu`` → ``surfradpy.scripts.aodinversion:produce_aodinversion1_01_catchup``
- ``srf_aodrealtime`` → ``surfradpy.scripts.produce_realtime_aod:routine``
- ``produce_aod_realtime_catchup`` → ``surfradpy.scripts.produce_realtime_aod:catchup``
- ``radflux2netcdf`` → ``surfradpy.scripts.radflux2netcdf:main``
- ``radiation2netcdf`` → ``surfradpy.scripts.radiation2netcdf:main``
- ``srf_langleys`` → ``surfradpy.scripts.createlangleys:main``

Notes
-----

- These commands are registered in ``pyproject.toml`` under ``[project.scripts]``.
- See the corresponding modules in ``surfradpy/scripts`` for usage details.

produce_surfrad_mfrsr2nc
------------------------

Convert SURFRAD MFRSR raw ``.xmd`` files to daily NetCDF outputs. The script
loops over standard SURFRAD sites and processes the last 60 days by default.

Defaults (from ``surfradpy.scripts.surfrad_mfrsr_raw2netcdf``):

- ``prefix``: ``/nfs/grad/``
- ``db_path``: ``/home/grad/htelg/prog/SURFRAD/notebooks/databases/surfrad_database.db``
- ``log_folder``: ``/home/grad/htelg/.processlogs/``
- ``test``: ``False`` (if ``True``, processes a single site without saving)

This entry point does not parse command-line arguments. To change inputs,
edit the defaults in the script or call ``surfradpy.scripts.surfrad_mfrsr_raw2netcdf.run``
directly from Python with your own parameters.

produce_aod2netcdf
------------------

Convert John's AOD product to a unified NetCDF output and rsync results if any
files were processed.

Defaults (from ``surfradpy.scripts.aod2netcdf``):

- ``site``: ``all``
- ``path2basefld_in``: ``/nfs/grad/surfrad/aod/``
- ``path2basefld_out``: ``/nfs/grad/surfrad/products_level2/aod_netcdf/v{version}/``
- Processes the last 3000 workplan entries.
- Rsync target: ``/nfs/iftp/aftp/g-rad/surfrad/aod_netcdf/v{version}/``

This entry point does not parse command-line arguments. To change inputs,
edit the defaults in the script or call ``surfradpy.scripts.aod2netcdf.run``
directly from Python.

qcrad2ncei
----------

Process QCRad data to NCEI requirements, including NetCDF conversion, archives
(tar), and manifest generation. This command provides CLI arguments.

Configuration
^^^^^^^^^^^^^

The script reads defaults from ``~/.SurfRadPy/config.ini`` (and creates it if
missing). You can update defaults with ``--update_defaults``.

Key arguments
^^^^^^^^^^^^^

- ``-i / --folder_in``: input QCRad base folder
- ``-o / --folder_out``: output NetCDF base folder
- ``-a / --folder_tar``: output archive base folder
- ``-s / --station``: station abbreviation (or all)
- ``-y / --year``: year to process
- ``-m / --month``: month to process
- ``-w / --overwrite``: overwrite NetCDF outputs
- ``-t / --test``: dry-run
- ``-v / --verbose``: verbose output
- ``--suppress_netcdf`` / ``--suppress_archive`` / ``--suppress_manifest``

Run ``qcrad2ncei --help`` to see the full CLI.

aodinv1
-------

Run AOD inversion (version 1.0) for a limited set of sites with fixed defaults.

Defaults (from ``surfradpy.scripts.aodinversion.produce_aodinversion1_01``):

- ``channels``: ``[415, 500, 670, 870, 1625]``
- ``sites``: ``['tbl']``
- ``start``: ``20160101``
- ``p2fldaod``: ``/nfs/grad/surfrad/products_level2/aod_netcdf/v1.0/``
- ``p2fldout``: ``/nfs/grad/surfrad/products_level2/aodinversion/1.0/``
- Workplan is sorted descending and truncated to 1 entry.

No CLI arguments. Edit the script or call the function directly for other
parameters.

aodinv1cu
---------

Catch-up run of AOD inversion (version 1.0) for multiple sites.

Defaults (from ``surfradpy.scripts.aodinversion.produce_aodinversion1_01_catchup``):

- ``channels``: ``[415, 500, 670, 870, 1625]``
- ``sites``: ``['bon', 'dra', 'fpk', 'gwn', 'sxf', 'psu']``
- ``start``: ``20160101``
- ``p2fldaod``: ``/nfs/grad/surfrad/products_level2/aod_netcdf/v1.0/``
- ``p2fldout``: ``/nfs/grad/surfrad/products_level2/aodinversion/1.0/``

No CLI arguments. Edit the script or call the function directly for other
parameters.

srf_aodrealtime
---------------

Run the realtime AOD pipeline for all sites.

Behavior (from ``surfradpy.scripts.produce_realtime_aod.routine``):

- Uses ``surfradpy.realtime_aod.mfrsr_AOD_lev0()`` and runs ``process_all()``.
- Logs via ``productomator`` automation.

produce_aod_realtime_catchup
----------------------------

Run the realtime AOD pipeline for a subset of sites (catch-up).

Behavior (from ``surfradpy.scripts.produce_realtime_aod.catchup``):

- Uses ``surfradpy.realtime_aod.mfrsr_AOD_lev0(site=['psu'])`` and runs
  ``process_all()``.
- Logs via ``productomator`` automation.

radflux2netcdf
--------------

Convert RadFlux inputs to daily NetCDF outputs.

Defaults (from ``surfradpy.scripts.radflux2netcdf.main``):

- ``path2fld_in``: ``/nfs/iftp/aftp/g-rad/surfrad/RadFlux/``
- ``path2fld_out``: ``/nfs/grad/surfrad/products_level4/radflux/v{version}/``
- ``sites``: ``['tbl', 'dra', 'fpk', 'gwn', 'psu', 'sxf', 'bon']``
- ``start``: ``180 days``
- ``overwrite``: ``False``

No CLI arguments. Edit the script or call the function directly for other
parameters.

radiation2netcdf
----------------

Generate NetCDF files from SURFRAD radiation data.

Defaults (from ``surfradpy.scripts.radiation2netcdf.main``):

- ``p2fld``: ``/nfs/iftp/aftp/data/radiation/surfrad/``
- ``p2fldout``: ``/nfs/grad/surfrad/products_level1/radiation_netcdf/``

No CLI arguments. Edit the script or call the function directly for other
parameters.

srf_langleys
------------

Run the Langley processing steps for realtime AOD data.

Behavior (from ``surfradpy.scripts.createlangleys.main``):

- Calls ``process_langleys`` with ``raise_error=True`` and ``verbose=True``.
- Then runs ``process_langley_concat`` to merge outputs.
