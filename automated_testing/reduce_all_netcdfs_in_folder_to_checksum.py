#!/usr/bin/env python3
"""
Reduce NetCDF files to checksum-only NetCDF files.

This script is the Python equivalent of:
- reduce_all_netcdfs_in_folder_to_checksum.m
- reduce_netcdf_to_checksum.m

Usage:
    python3 reduce_all_netcdfs_in_folder_to_checksum.py <foldername>

Behavior:
- Finds all *.nc files in <foldername>
- Excludes *_checksum.nc and resource_tracking.nc
- For each file:
  - Creates <name>_checksum.nc containing one 4-value checksum variable per
    original variable (including variables in groups)
  - Deletes the original file
"""

import os
import sys
from typing import List, Tuple

import numpy as np

try:
    import netCDF4 as nc
except ImportError:
    print("ERROR: netCDF4 package is required. Install with: pip install netCDF4")
    sys.exit(2)


CHECKSUM_DIM_NAME = "checksum"
CHECKSUM_DIM_LEN = 4
FILL_VALUE = 9e9


def _collect_variables(group, parent_name: str = "") -> List[Tuple[str, object, str]]:
    """
    Collect all variables in a dataset/group recursively.

    Returns tuples of:
    - original_name: original path used for reading, e.g. "grp/subgrp/var"
    - var_obj: netCDF4.Variable object from source file
    - new_name: flattened checksum variable name, e.g. "grp_subgrp_var"
    """
    found = []

    for var_name, var in group.variables.items():
        if parent_name:
            original_name = f"{parent_name}/{var_name}"
            new_name = f"{parent_name}_{var_name}"
        else:
            original_name = var_name
            new_name = var_name
        found.append((original_name, var, new_name))

    for grp_name, grp in group.groups.items():
        group_path = f"{parent_name}_{grp_name}" if parent_name else grp_name
        found.extend(_collect_variables(grp, group_path))

    return found


def reduce_netcdf_to_checksum(filename: str) -> str:
    """Reduce data in a NetCDF file to checksums. Returns checksum filename."""
    if not os.path.isfile(filename):
        raise FileNotFoundError(f'Could not find file "{filename}"')

    filename_redux = f"{filename[:-3]}_checksum.nc" if filename.endswith(".nc") else f"{filename}_checksum.nc"

    if os.path.isfile(filename_redux):
        os.remove(filename_redux)

    with nc.Dataset(filename, "r") as src, nc.Dataset(filename_redux, "w") as dst:
        # Keep original global dimensions and add checksum dimension (Matlab behavior).
        for dim_name, dim in src.dimensions.items():
            dst.createDimension(dim_name, None if dim.isunlimited() else len(dim))
        dst.createDimension(CHECKSUM_DIM_NAME, CHECKSUM_DIM_LEN)

        variables = _collect_variables(src)

        for original_name, src_var, new_name in variables:
            dst_var = dst.createVariable(
                new_name,
                "f8",
                (CHECKSUM_DIM_NAME,),
                fill_value=FILL_VALUE,
                zlib=False,
                shuffle=False,
                chunksizes=(CHECKSUM_DIM_LEN,),
            )

            # Copy original variable attributes to preserve metadata.
            for attr_name in src_var.ncattrs():
                if attr_name == "_FillValue":
                    continue
                dst_var.setncattr(attr_name, src_var.getncattr(attr_name))

            data = np.ma.filled(src[original_name][:], fill_value=0.0)
            flat = np.asarray(data).reshape(-1)

            has_non_finite = np.any(~np.isfinite(flat))
            if has_non_finite:
                checksum_sum = FILL_VALUE
                checksum_sum_abs = FILL_VALUE
                finite_flat = flat[np.isfinite(flat)]
                if finite_flat.size > 0:
                    checksum_min = np.min(finite_flat)
                    checksum_max = np.max(finite_flat)
                else:
                    checksum_min = np.nan
                    checksum_max = np.nan
            else:
                checksum_sum = np.sum(flat)
                checksum_sum_abs = np.sum(np.abs(flat))
                checksum_min = np.min(flat)
                checksum_max = np.max(flat)

            # Keep min/max behavior, but avoid noisy warnings for non-finite
            # reductions like +inf + -inf.
            with np.errstate(invalid="ignore", over="ignore"):
                checksum_data = np.array(
                    [
                        checksum_sum,
                        checksum_sum_abs,
                        checksum_min,
                        checksum_max,
                    ],
                    dtype=np.float64,
                )

            if has_non_finite:
                checksum = np.ma.array(checksum_data, mask=[True, True, False, False])
                dst_var[:] = checksum
            else:
                dst_var[:] = checksum_data

    return filename_redux


def reduce_all_netcdfs_in_folder_to_checksum(foldername: str) -> None:
    """Reduce all non-checksum NetCDF files in a folder and delete originals."""
    if not os.path.isdir(foldername):
        raise FileNotFoundError(f'Could not find folder "{foldername}"')

    files = sorted(
        name
        for name in os.listdir(foldername)
        if name.endswith(".nc")
        and not name.endswith("_checksum.nc")
        and name.lower() != "resource_tracking.nc"
    )

    for i, name in enumerate(files, start=1):
        print(f"Reducing file {i}/{len(files)}: {name}")

        filename = os.path.join(foldername, name)
        reduce_netcdf_to_checksum(filename)
        os.remove(filename)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <foldername>")
        sys.exit(2)

    try:
        reduce_all_netcdfs_in_folder_to_checksum(sys.argv[1])
    except Exception as exc:
        print(f"ERROR: {exc}")
        sys.exit(1)

    sys.exit(0)
