#!/usr/bin/env python3
"""
Show program info attributes stored in NetCDF files for a test's reference and results.

Usage:
    python3 compare_program_info.py <test_folder>

This mirrors automated_testing/compare_program_info.m: it reads a fixed list of
global attributes from the first NetCDF file found in reference/ and results/ and
prints both values side by side. Missing attributes are allowed and shown as empty
strings to preserve compatibility with older references.
"""

import os
import sys

try:
    import netCDF4 as nc
except ImportError:
    print("ERROR: netCDF4 package is required. Install with: pip install netCDF4")
    sys.exit(2)


INFOS = [
    "git commit hash",
    "PETSc version",
    "NetCDF version",
    "OpenMPI version",
    "compiler version",
    "compiler flags",
]


def find_first_netcdf_file_in_dir(foldername: str) -> str:
    for filename_short in os.listdir(foldername):
        if filename_short.endswith(".nc"):
            return os.path.join(foldername, filename_short)
    raise FileNotFoundError(f"Could not find a NetCDF file in {foldername}")


def read_program_info(foldername: str, fieldnames: list[str]) -> dict[str, str]:
    filename = find_first_netcdf_file_in_dir(foldername)
    program_info: dict[str, str] = {}

    with nc.Dataset(filename, "r") as dataset:
        for attname in fieldnames:
            # Older references might not have all the program info; allow this.
            try:
                attval = dataset.getncattr(attname)
            except AttributeError:
                attval = ""

            fieldname = attname.replace(" ", "_")
            if isinstance(attval, bytes):
                program_info[fieldname] = attval.decode(errors="replace")
            else:
                program_info[fieldname] = str(attval)

    return program_info


def show_both_program_infos(foldername: str, program_info_ref: dict[str, str], program_info_mod: dict[str, str]) -> None:
    print()
    print(f"Program info for {foldername}")
    print()

    for fieldname in program_info_ref:
        field_ref = program_info_ref[fieldname]
        field_mod = program_info_mod[fieldname]

        print(f" {fieldname}")
        print(f"  Reference: {field_ref}")
        print(f"  Results  : {field_mod}")


def compare_program_info(foldername: str) -> None:
    program_info_mod = read_program_info(os.path.join(foldername, "results"), INFOS)
    program_info_ref = read_program_info(os.path.join(foldername, "reference"), INFOS)
    show_both_program_infos(foldername, program_info_ref, program_info_mod)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <test_folder>")
        sys.exit(2)

    try:
        compare_program_info(sys.argv[1])
    except Exception as exc:
        print(f"ERROR: {exc}")
        sys.exit(1)

    sys.exit(0)