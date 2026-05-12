#!/usr/bin/env python3
"""
Reduce NetCDF files to checksum-only NetCDF files.

This script is the Python equivalent of:
- reduce_all_netcdfs_in_folder_to_checksum.m
- reduce_netcdf_to_checksum.m

Usage:
        python3 reduce_all_netcdfs_in_folder_to_checksum.py <test_folder>

Behavior:
- Reads all *.nc files from <test_folder>/results
- Excludes *_checksum.nc and resource_tracking.nc
- Creates <test_folder>/results_checksum if needed
- For each file:
    - Creates results_checksum/<name>_checksum.nc containing:
        - one 4-value stats checksum variable per original variable
            [sum_finite, sum_abs_finite, min_finite, max_finite]
        - one 3-value integer count checksum variable per original variable
            [n_NaN, n_Inf, n_mInf]
        Original variables in groups are flattened with underscores in the name.
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
CHECKSUM_COUNT_DIM_NAME = "checksum_count"
CHECKSUM_COUNT_DIM_LEN = 3
FILL_VALUE = 9e9


def _collect_variables(group, parent_name: str = "") -> List[Tuple[str, object, str]]:
    """
    Collect all variables in a dataset/group recursively.

    Returns tuples of:
    - original_name: original path used for reading, e.g. "grp/subgrp/var"
    - var_obj: netCDF4.Variable object from source file
    - new_name: flattened checksum variable base name, e.g. "grp_subgrp_var"
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

    # Intentional recursion: NetCDF groups can be nested arbitrarily deep.
    for grp_name, grp in group.groups.items():
        group_path = f"{parent_name}_{grp_name}" if parent_name else grp_name
        found.extend(_collect_variables(grp, group_path))

    return found


def reduce_netcdf_to_checksum(filename: str, output_dir: str) -> str:
    """Reduce data in a NetCDF file to checksums. Returns checksum filename."""
    if not os.path.isfile(filename):
        raise FileNotFoundError(f'Could not find file "{filename}"')

    basename = os.path.basename(filename)
    if basename.endswith(".nc"):
        basename_redux = f"{basename[:-3]}_checksum.nc"
    else:
        basename_redux = f"{basename}_checksum.nc"
    filename_redux = os.path.join(output_dir, basename_redux)

    if os.path.isfile(filename_redux):
        os.remove(filename_redux)

    with nc.Dataset(filename, "r") as src, nc.Dataset(filename_redux, "w") as dst:
        # Keep original global dimensions and add checksum dimension (Matlab behavior).
        for dim_name, dim in src.dimensions.items():
            dst.createDimension(dim_name, None if dim.isunlimited() else len(dim))
        dst.createDimension(CHECKSUM_DIM_NAME, CHECKSUM_DIM_LEN)
        dst.createDimension(CHECKSUM_COUNT_DIM_NAME, CHECKSUM_COUNT_DIM_LEN)

        for attr_name in src.ncattrs():
            dst.setncattr(attr_name, src.getncattr(attr_name))

        variables = _collect_variables(src)

        for original_name, src_var, new_name in variables:
            is_integer_source = np.issubdtype(src_var.dtype, np.integer)

            stats_dtype = "i8" if is_integer_source else "f8"
            dst_var = dst.createVariable(
                new_name,
                stats_dtype,
                (CHECKSUM_DIM_NAME,),
                fill_value=FILL_VALUE,
                zlib=False,
                shuffle=False,
                chunksizes=(CHECKSUM_DIM_LEN,),
            )

            counts_var = dst.createVariable(
                f"{new_name}_counts",
                "i8",
                (CHECKSUM_COUNT_DIM_NAME,),
                zlib=False,
                shuffle=False,
                chunksizes=(CHECKSUM_COUNT_DIM_LEN,),
            )

            # Copy original variable attributes to preserve metadata.
            for attr_name in src_var.ncattrs():
                if attr_name == "_FillValue":
                    continue
                dst_var.setncattr(attr_name, src_var.getncattr(attr_name))

            data = src[original_name][:]
            flat = np.asarray(data).reshape(-1)

            # Masked entries should not contribute to checksums; map them to NaN.
            if np.ma.isMaskedArray(data):
                mask = np.ma.getmaskarray(data).reshape(-1)
                if np.any(mask):
                    flat = flat.astype(np.float64, copy=False)
                    flat[mask] = np.nan

            finite_flat = flat[np.isfinite(flat)]

            # Existing checksums, all computed over finite values only.
            checksum_sum = np.sum(finite_flat)
            checksum_sum_abs = np.sum(np.abs(finite_flat))
            if finite_flat.size > 0:
                checksum_min = np.min(finite_flat)
                checksum_max = np.max(finite_flat)
            else:
                checksum_min = np.nan
                checksum_max = np.nan

            # Additional checksums that track non-finite values in original data.
            checksum_n_nan = np.sum(np.isnan(flat))
            checksum_n_inf = np.sum(np.isposinf(flat))
            checksum_n_minf = np.sum(np.isneginf(flat))

            # Guard against warning noise if upstream arrays contain non-finite values.
            with np.errstate(invalid="ignore", over="ignore"):
                stats_data = np.array(
                    [
                        checksum_sum,
                        checksum_sum_abs,
                        checksum_min,
                        checksum_max,
                    ],
                    dtype=np.float64,
                )

            counts_data = np.array(
                [checksum_n_nan, checksum_n_inf, checksum_n_minf],
                dtype=np.int64,
            )

            if is_integer_source:
                # Integer source variables get integer-typed sum/sum_abs/min/max.
                stats_out = np.array(np.rint(stats_data), dtype=np.int64)
            else:
                stats_out = stats_data

            dst_var[:] = stats_out
            counts_var[:] = counts_data

    return filename_redux


def reduce_all_netcdfs_in_folder_to_checksum(test_folder: str) -> None:
    """Reduce all non-checksum NetCDF files in test_folder/results to test_folder/results_checksum."""
    if not os.path.isdir(test_folder):
        raise FileNotFoundError(f'Could not find test folder "{test_folder}"')

    results_dir = os.path.join(test_folder, "results")
    if not os.path.isdir(results_dir):
        raise FileNotFoundError(f'Could not find results folder "{results_dir}"')

    results_checksum_dir = os.path.join(test_folder, "results_checksum")
    if os.path.exists(results_checksum_dir):
        raise FileExistsError(f'Folder already exists: "{results_checksum_dir}"')
    os.makedirs(results_checksum_dir)

    files = sorted(
        name
        for name in os.listdir(results_dir)
        if name.endswith(".nc")
        and not name.endswith("_checksum.nc")
        and name.lower() != "resource_tracking.nc"
    )

    for i, name in enumerate(files, start=1):
        print(f"Reducing file {i}/{len(files)}: {name}")

        filename = os.path.join(results_dir, name)
        reduce_netcdf_to_checksum(filename, results_checksum_dir)


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
