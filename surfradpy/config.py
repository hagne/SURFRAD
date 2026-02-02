"""
Configuration helpers for SURFRADPY.
"""

import configparser
import os
import pathlib as pl


def get_config_path() -> pl.Path:
    """
    Return the default config file path.
    """
    env_path = os.environ.get("SURFRAD_CONFIG_PATH")
    if env_path:
        return pl.Path(env_path).expanduser()
    return pl.Path.home() / ".config" / "surfradpy" / "config.ini"


def load_config(path: pl.Path | None = None) -> configparser.ConfigParser:
    """
    Load config from disk. Returns an empty ConfigParser if the file is missing.
    """
    cfg = configparser.ConfigParser()
    cfg.read((path or get_config_path()))
    return cfg


def get_db_path() -> pl.Path | None:
    """
    Resolve the SURFRAD database path from environment or config.

    Precedence:
    1) SURFRAD_DB_PATH environment variable
    2) [database] path in config.ini
    """
    env_path = os.environ.get("SURFRAD_DB_PATH")
    if env_path:
        return pl.Path(env_path).expanduser()
    cfg = load_config()
    if cfg.has_option("database", "path"):
        return pl.Path(cfg.get("database", "path")).expanduser()
    return None
