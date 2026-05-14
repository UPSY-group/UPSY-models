#!/usr/bin/env python3
"""
Compare all checksum NetCDF files in a test folder's reference/ vs results_checksum/ sub-directories.

Usage:
    python compare_all_netcdfs_in_test_folder.py <foldername>

Exits with code 0 if all files match, 1 otherwise.
"""

import sys
import os
import numpy as np

try:
    import netCDF4 as nc
except ImportError:
    print("ERROR: netCDF4 package is required. Install with: pip install netCDF4")
    sys.exit(2)


TOL = 1e-12


# ---------------------------------------------------------------------------
# compare_netcdf
# ---------------------------------------------------------------------------

def compare_netcdf(filename_ref: str, filename_mod: str) -> bool:
    """Compare two NetCDF files and print any differences. Returns True if identical."""

    basename_ref = os.path.basename(filename_ref)
    basename_mod = os.path.basename(filename_mod)
    if basename_ref.lower() != basename_mod.lower():
        print(f'  Mismatching filenames: reference = "{basename_ref}", model = "{basename_mod}"')
        return False

    with nc.Dataset(filename_ref, "r") as ds_ref, nc.Dataset(filename_mod, "r") as ds_mod:
        info_ok = _compare_group(basename_ref, ds_ref, ds_mod, compare_attrs=False)
        if not info_ok:
            return False
        data_ok = _compare_data_group(filename_ref, filename_mod, basename_ref, ds_ref, ds_mod)

    return data_ok


# ---------------------------------------------------------------------------
# Structure comparison (dimensions, variables, attributes, groups)
# ---------------------------------------------------------------------------

def _compare_group(parent_name: str, grp_ref, grp_mod, compare_attrs: bool = True) -> bool:
    ok = True
    ok = _compare_dimensions(parent_name, grp_ref, grp_mod) and ok
    ok = _compare_variables(parent_name, grp_ref, grp_mod) and ok
    if compare_attrs:
        ok = _compare_attributes(parent_name, grp_ref, grp_mod) and ok
    ok = _compare_groups(parent_name, grp_ref, grp_mod) and ok
    return ok


def _compare_dimensions(parent_name: str, grp_ref, grp_mod) -> bool:
    dims_ref = grp_ref.dimensions
    dims_mod = grp_mod.dimensions

    # Variables expose .dimensions as a tuple of name strings; groups expose a dict
    if isinstance(dims_ref, tuple):
        if dims_ref != dims_mod:
            print(f'  Mismatching dimensions for "{parent_name}": '
                  f'reference = {dims_ref}, model = {dims_mod}')
            return False
        return True

    if len(dims_ref) != len(dims_mod):
        print(f"  Mismatching number of dimensions in {parent_name}: "
              f"reference = {len(dims_ref)}, model = {len(dims_mod)}")
        return False

    ok = True
    for name_ref, dim_ref in dims_ref.items():
        if name_ref not in dims_mod:
            print(f'  Dimension "{parent_name}/{name_ref}" missing from model')
            ok = False
            continue
        dim_mod = dims_mod[name_ref]

        if dim_ref.size != dim_mod.size:
            ok = False
            print(f'  Mismatching dimension "{parent_name}/{name_ref}" length: '
                  f'reference = {dim_ref.size}, model = {dim_mod.size}')

        if dim_ref.isunlimited() != dim_mod.isunlimited():
            ok = False
            print(f'  Mismatching dimension "{parent_name}/{name_ref}" unlimitedness: '
                  f'reference = {dim_ref.isunlimited()}, model = {dim_mod.isunlimited()}')

    return ok


def _compare_attributes(parent_name: str, grp_ref, grp_mod) -> bool:
    atts_ref = {a: grp_ref.getncattr(a) for a in grp_ref.ncattrs()}
    atts_mod = {a: grp_mod.getncattr(a) for a in grp_mod.ncattrs()}

    if len(atts_ref) != len(atts_mod):
        print(f"  Mismatching number of attributes in {parent_name}: "
              f"reference = {len(atts_ref)}, model = {len(atts_mod)}")
        return False

    ok = True
    for name_ref, val_ref in atts_ref.items():
        if name_ref not in atts_mod:
            print(f'  Attribute "{parent_name}/{name_ref}" missing from model')
            ok = False
            continue
        val_mod = atts_mod[name_ref]

        if isinstance(val_ref, str) or isinstance(val_mod, str):
            if str(val_ref).lower() != str(val_mod).lower():
                ok = False
                print(f'  Mismatching attribute "{parent_name}/{name_ref}": '
                      f'reference = "{val_ref}", model = "{val_mod}"')
        else:
            if np.any(np.asarray(val_ref) != np.asarray(val_mod)):
                ok = False
                print(f'  Mismatching attribute "{parent_name}/{name_ref}": '
                      f'reference = {val_ref}, model = {val_mod}')

    return ok


def _compare_variables(parent_name: str, grp_ref, grp_mod) -> bool:
    vars_ref = grp_ref.variables
    vars_mod = grp_mod.variables

    if len(vars_ref) != len(vars_mod):
        print(f"  Mismatching number of variables in {parent_name}: "
              f"reference = {len(vars_ref)}, model = {len(vars_mod)}")
        return False

    ok = True
    for name_ref, var_ref in vars_ref.items():
        if name_ref not in vars_mod:
            print(f'  Variable "{parent_name}/{name_ref}" missing from model')
            ok = False
            continue
        var_mod = vars_mod[name_ref]
        var_path = f"{parent_name}/{name_ref}"

        # Shape / size
        if var_ref.shape != var_mod.shape:
            ok = False
            print(f'  Mismatching variable "{var_path}" size: '
                  f'reference = {list(var_ref.shape)}, model = {list(var_mod.shape)}')

        # Datatype
        if str(var_ref.dtype).lower() != str(var_mod.dtype).lower():
            ok = False
            print(f'  Mismatching variable "{var_path}" datatype: '
                  f'reference = "{var_ref.dtype}", model = "{var_mod.dtype}"')

        # Dimensions
        ok = _compare_dimensions(var_path, var_ref, var_mod) and ok

        # Attributes
        ok = _compare_attributes(var_path, var_ref, var_mod) and ok

        # Chunking
        chunking_ref = var_ref.chunking()
        chunking_mod = var_mod.chunking()
        if chunking_ref != chunking_mod:
            ok = False
            print(f'  Mismatching variable "{var_path}" ChunkSize: '
                  f'reference = {chunking_ref}, model = {chunking_mod}')

        # Fill value
        fill_ref = var_ref._FillValue if hasattr(var_ref, '_FillValue') else None
        fill_mod = var_mod._FillValue if hasattr(var_mod, '_FillValue') else None
        if fill_ref != fill_mod:
            ok = False
            print(f'  Mismatching variable "{var_path}" FillValue: '
                  f'reference = {fill_ref}, model = {fill_mod}')

        # Compression (deflate level, shuffle) - netCDF4-python API
        filters_ref = var_ref.filters() if hasattr(var_ref, 'filters') else {}
        filters_mod = var_mod.filters() if hasattr(var_mod, 'filters') else {}
        if filters_ref.get('complevel') != filters_mod.get('complevel'):
            ok = False
            print(f'  Mismatching variable "{var_path}" DeflateLevel: '
                  f'reference = {filters_ref.get("complevel")}, model = {filters_mod.get("complevel")}')
        if filters_ref.get('shuffle') != filters_mod.get('shuffle'):
            ok = False
            print(f'  Mismatching variable "{var_path}" Shuffle: '
                  f'reference = {filters_ref.get("shuffle")}, model = {filters_mod.get("shuffle")}')

    return ok


def _compare_groups(parent_name: str, grp_ref, grp_mod) -> bool:
    groups_ref = grp_ref.groups
    groups_mod = grp_mod.groups

    if len(groups_ref) == 0 and len(groups_mod) == 0:
        return True

    if len(groups_ref) != len(groups_mod):
        print(f"  Mismatching number of groups in {parent_name}: "
              f"reference = {len(groups_ref)}, model = {len(groups_mod)}")
        return False

    ok = True
    for name_ref, sub_ref in groups_ref.items():
        if name_ref not in groups_mod:
            print(f'  Group "{parent_name}/{name_ref}" missing from model')
            ok = False
            continue
        sub_mod = groups_mod[name_ref]
        group_path = f"{parent_name}/{name_ref}"
        ok = _compare_group(group_path, sub_ref, sub_mod) and ok

    return ok


# ---------------------------------------------------------------------------
# Data comparison
# ---------------------------------------------------------------------------

def _compare_data_group(filename_ref: str, filename_mod: str,
                        parent_name: str, grp_ref, grp_mod) -> bool:
    ok = True

    for var_name in grp_ref.variables:
        var_path = f"{parent_name}/{var_name}"
        d_ref = grp_ref.variables[var_name][:].flatten()
        d_mod = grp_mod.variables[var_name][:].flatten()
        ok = _compare_data_variable(var_path, d_ref, d_mod, filename_ref) and ok

    for grp_name, sub_ref in grp_ref.groups.items():
        sub_mod = grp_mod.groups[grp_name]
        sub_path = f"{parent_name}/{grp_name}"
        ok = _compare_data_group(filename_ref, filename_mod, sub_path, sub_ref, sub_mod) and ok

    return ok


def _compare_data_variable(var_path: str, d_ref: np.ndarray,
                            d_mod: np.ndarray, filename_ref: str) -> bool:
    is_checksum = '_checksum' in os.path.basename(filename_ref)
    is_count_checksum = var_path.endswith('_counts')

    def _both_nan(a: float, b: float) -> bool:
        return np.isnan(a) and np.isnan(b)

    def _scalar_or_nan(x) -> float:
        if np.ma.is_masked(x):
            return np.nan
        return float(x)

    def _scalar_int_or_none(x):
        y = _scalar_or_nan(x)
        if np.isnan(y):
            return None
        return int(round(y))

    if is_checksum and is_count_checksum and len(d_ref) == 3:
        vals_ref = [_scalar_int_or_none(v) for v in d_ref]
        vals_mod = [_scalar_int_or_none(v) for v in d_mod]

        if vals_ref != vals_mod:
            print(f'  Mismatching data in {var_path}')
            print(f'    Reference: {vals_ref}')
            print(f'    Model    : {vals_mod}')
            return False
        return True

    if is_checksum and (not is_count_checksum) and len(d_ref) == 4:
        # Values are: [sum_finite, sum_abs_finite, min_finite, max_finite]
        vals_ref = [_scalar_or_nan(v) for v in d_ref]
        vals_mod = [_scalar_or_nan(v) for v in d_mod]

        ref_sum, ref_sum_abs, ref_min, ref_max = vals_ref
        mod_sum, mod_sum_abs, mod_min, mod_max = vals_mod

        ref_min_max_abs = max(abs(ref_min), abs(ref_max))

        mismatches = []

        if _both_nan(ref_sum, mod_sum):
            pass
        elif np.isnan(ref_sum) or np.isnan(mod_sum):
            mismatches.append(('sum', np.nan))
        elif ref_sum_abs != 0:
            diff_sum = abs((mod_sum - ref_sum) / ref_sum_abs)
            if diff_sum > TOL:
                mismatches.append(('sum', diff_sum))
        elif abs(mod_sum - ref_sum) > TOL:
            mismatches.append(('sum', abs(mod_sum - ref_sum)))

        if _both_nan(ref_sum_abs, mod_sum_abs):
            pass
        elif np.isnan(ref_sum_abs) or np.isnan(mod_sum_abs):
            mismatches.append(('sum_abs', np.nan))
        elif ref_sum_abs != 0:
            diff_sum_abs = abs(1.0 - mod_sum_abs / ref_sum_abs)
            if diff_sum_abs > TOL:
                mismatches.append(('sum_abs', diff_sum_abs))
        elif abs(mod_sum_abs - ref_sum_abs) > TOL:
            mismatches.append(('sum_abs', abs(mod_sum_abs - ref_sum_abs)))

        if _both_nan(ref_min, mod_min):
            pass
        elif np.isnan(ref_min) or np.isnan(mod_min):
            mismatches.append(('min', np.nan))
        elif ref_min_max_abs != 0:
            diff_min = abs((mod_min - ref_min) / ref_min_max_abs)
            if diff_min > TOL:
                mismatches.append(('min', diff_min))
        elif abs(mod_min - ref_min) > TOL:
            mismatches.append(('min', abs(mod_min - ref_min)))

        if _both_nan(ref_max, mod_max):
            pass
        elif np.isnan(ref_max) or np.isnan(mod_max):
            mismatches.append(('max', np.nan))
        elif ref_min_max_abs != 0:
            diff_max = abs((mod_max - ref_max) / ref_min_max_abs)
            if diff_max > TOL:
                mismatches.append(('max', diff_max))
        elif abs(mod_max - ref_max) > TOL:
            mismatches.append(('max', abs(mod_max - ref_max)))

        if mismatches:
            print(f'  Mismatching data in {var_path}')
            print(f'    Reference: {list(d_ref)}')
            print(f'    Model    : {list(d_mod)}')
            return False
        return True

    else:
        # General comparison
        with np.errstate(divide='ignore', invalid='ignore'):
            dd = np.where(d_ref != 0, np.abs(1.0 - d_mod / d_ref), np.abs(d_mod - d_ref))
        if not np.all(dd < TOL):
            print(f'  Mismatching data in {var_path}')
            return False
        return True


# ---------------------------------------------------------------------------
# compare_all_netcdfs_in_test_folder
# ---------------------------------------------------------------------------

def compare_all_netcdfs_in_test_folder(foldername: str) -> bool:
    """Compare all *_checksum.nc files in reference/ vs results_checksum/. Returns True if all match."""

    if not os.path.isdir(foldername):
        print(f'ERROR: Could not find test "{foldername}"')
        return False

    foldername_ref = os.path.join(foldername, 'reference')
    if not os.path.isdir(foldername_ref):
        print(f'ERROR: Could not find reference for test "{foldername}"')
        return False

    foldername_mod = os.path.join(foldername, 'results_checksum')
    if not os.path.isdir(foldername_mod):
        print(f'ERROR: Could not find results_checksum for test "{foldername}"')
        return False

    checksum_files = sorted(
        f for f in os.listdir(foldername_ref)
        if f.endswith('_checksum.nc')
    )

    all_match = True
    n = len(checksum_files)
    for i, fname in enumerate(checksum_files, start=1):
        print(f'Comparing netcdf file {i}/{n}: {fname}')

        filename_ref = os.path.join(foldername_ref, fname)
        filename_mod = os.path.join(foldername_mod, fname)

        if not os.path.isfile(filename_mod):
            print(f'ERROR: File "{fname}" does not exist in results_checksum of test "{foldername}"')
            all_match = False
            continue

        if not compare_netcdf(filename_ref, filename_mod):
            all_match = False

    return all_match


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f'Usage: {sys.argv[0]} <foldername>')
        sys.exit(2)

    success = compare_all_netcdfs_in_test_folder(sys.argv[1])

    if not success:
        print('ERROR: Not all files are identical')
        sys.exit(1)

    print('All files are identical')
    sys.exit(0)
