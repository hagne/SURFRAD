Configuration
=============

This project uses a small config file for local, machine-specific settings.

Database path
-------------

Set the default SURFRAD database path in a config file at:

``~/.config/surfradpy/config.ini``

Example:

.. code-block:: ini

   [database]
   path = /path/to/surfrad_database.db

You can also override the config file location or the database path via
environment variables:

- ``SURFRAD_CONFIG_PATH``: path to the config file
- ``SURFRAD_DB_PATH``: direct path to the database file
