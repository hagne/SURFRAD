Installation
============

Conda (recommended)
-------------------

You will need to install two python package from their respective github repositories:
- atmPy: `https://github.com/htelg/atmPy`
- productomator: `https://github.com/htelg/productomator`

Use your currently active conda environment (create one if you want, but it is
your choice), then install the dependencies into it.

Core dependencies:

.. code-block:: console

   conda install -c conda-forge numpy pandas xarray netcdf4 setproctitle
   pip install .

Use -e keyword for editable/development install.

Alternatively, update an existing environment from the YAML file:

.. code-block:: console

   conda env update --name <env-name> -f environment.yml

Documentation dependencies
--------------------------

If you want to build the docs in the same environment, add the doc packages:

.. code-block:: console

   conda install -c conda-forge sphinx sphinx-rtd-theme

Or update an existing environment from the docs YAML file:

.. code-block:: console

   conda env update --name <env-name> -f environment-docs.yml

Then build the docs from the ``doc/`` directory:

.. code-block:: console

   make html

How to install the processing scripts on a server
--------------------------------------------------
There are many ways to install the processing scripts on a server. Here is one way to do it.
1) Create a designated conda environment on the server for each script that is supposed to be run. This way you create a stable environment that does not break when development of other parts of the code base happens.
2) Install the required packages into the environment as described above.
   a) Make sure that atmPy and productomator are installed with "pip install ."! Do not install in editable/development mode!
3) Install the surfradpy package into the environment with "pip install ." Do not install in editable/development mode!
4) Test that the script works
5) Set up a cron job and run with full path to executable of the conda environment, e.g.:
   
   .. code-block:: console

      /path/to/conda/envs/surfrad_mfrsr_raw2netcdf_env/bin/surfrad_mfrsr_raw2netcdf

