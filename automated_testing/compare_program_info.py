#!/usr/bin/env python3
"""
Show program info attributes stored in NetCDF files for a test's reference and results.

Usage:
    python3 compare_program_info.py <test_folder> [reference_subdir]

This reads a fixed list of global attributes from the first NetCDF file found in
reference/ and results/ and prints both values side by side.

Unlike earlier behavior, all listed attributes are required in both files.
If any required attribute is missing, the script raises an error so CI fails.
"""

import os
import sys

try:
    import netCDF4 as nc
except ImportError:
    print("ERROR: netCDF4 package is required. Install with: pip install netCDF4")
    sys.exit(2)


INFOS = [
    "git_commit_hash",
    "PETSc_version",
    "NetCDF_version",
    "OpenMPI_version",
    "compiler_version",
    "compiler_flags",
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
            try:
                attval = dataset.getncattr(attname)
            except AttributeError as exc:
                raise KeyError(f'Missing required global attribute "{attname}" in {filename}') from exc

            if isinstance(attval, bytes):
                program_info[attname] = attval.decode(errors="replace")
            else:
                program_info[attname] = str(attval)

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


def compare_program_info(foldername: str, reference_subdir: str = "reference") -> None:
    program_info_mod = read_program_info(os.path.join(foldername, "results"), INFOS)
    program_info_ref = read_program_info(os.path.join(foldername, reference_subdir), INFOS)
    show_both_program_infos(foldername, program_info_ref, program_info_mod)


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        print(f"Usage: {sys.argv[0]} <test_folder> [reference_subdir]")
        sys.exit(2)

    try:
        reference_subdir = sys.argv[2] if len(sys.argv) == 3 else "reference"
        compare_program_info(sys.argv[1], reference_subdir)
    except Exception as exc:
        print(f"ERROR: {exc}")
        sys.exit(1)

    sys.exit(0)