Installation
============

Conda (recommended)
-------------------

Use your currently active conda environment (create one if you want, but it is
your choice), then install the dependencies into it.

Core dependencies:

.. code-block:: console

   conda install -c conda-forge numpy pandas xarray netcdf4 setproctitle
   pip install atmPy productomator
   pip install -e .

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
