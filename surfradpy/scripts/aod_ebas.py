#!/usr/bin/env python3
"""Command line interface for SURFRAD AOD EBAS generation."""

import argparse
import inspect
import logging

from surfradpy.products import aod_ebas


def run(
    year: int,
    stations: list[str] | None = None,
    output: str | None = None,
    source: str = str(aod_ebas.DEFAULT_SOURCE),
    archive: str = str(aod_ebas.DEFAULT_ARCHIVE),
    debug: bool = False,
) -> None:
    """Generate SURFRAD AOD EBAS NASA Ames files and optionally upload them."""
    if debug:
        logging.basicConfig(level=logging.DEBUG)

    if not aod_ebas.is_available():
        raise RuntimeError(
            "EBAS library not available, install it from "
            "https://git.nilu.no/ebas/ebas-io/-/tree/master"
        )

    aod_ebas.process_aod_ebas(
        year=year,
        stations=stations,
        output=output,
        source=source,
        archive=archive,
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=inspect.getdoc(run) or "SURFRAD AOD EBAS converter and uploader.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--debug", dest="debug", action="store_true", help="enable debug output")
    parser.add_argument("year", help="year to process", type=int)
    parser.add_argument("station", help="station code", nargs="*")
    parser.add_argument("--output", help="set the output directory instead of uploading")
    parser.add_argument(
        "--source",
        help="set the source file directory",
        default=str(aod_ebas.DEFAULT_SOURCE),
    )
    parser.add_argument(
        "--archive",
        help="set the archive directory",
        default=str(aod_ebas.DEFAULT_ARCHIVE),
    )
    args = parser.parse_args(argv)

    try:
        run(
            year=args.year,
            stations=args.station or None,
            output=args.output,
            source=args.source,
            archive=args.archive,
            debug=args.debug,
        )
    except (RuntimeError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
