"""Remap data between UPSY meshes and x/y- or lon/lat-grids.

Remapping weights are calculated by the Fortran program
UPSY_calculate_remapping_weights and stored as sparse (CSR) matrices in
NetCDF files. As with the Atlas inside UPSY, every set of weights is
calculated only once: maps are cached in memory, and weights files are kept
in a directory of pre-calculated weights (default: UPSY_remapping_weights/ in
the repository root, or $UPSY_REMAPPING_WEIGHTS_DIR) so they can be reused
across sessions.

Weights files are identified by a hash of the variables that define the
source and destination grids (not by their filenames), plus the Arakawa
grids and the remapping method.

Usage:

    from upsy.remap import remap
    d_mesh = remap(d_grid, 'grid.nc', 'mesh.nc', 'first_order_conservative')

Data on x/y- and lon/lat-grids is passed with the grid in the last two axes,
ordered (y, x) or (lat, lon) as it is read from a NetCDF file by netCDF4 or
xarray. Data on a mesh is passed with the vertices/triangles/edges in the
last axis. Any leading axes (time, depth, ...) are remapped in one go.

Destination cells for which the weights contain no entries (i.e. for which
UPSY_calculate_remapping_weights found no overlapping source cells) are set
to NaN.
"""

import hashlib
import os
import re
import shutil
import subprocess
import uuid
import warnings
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import netCDF4
import numpy as np
import scipy.sparse as sp

METHODS = ('nearest_neighbour', 'linear', 'first_order_conservative', 'second_order_conservative')
ARAKAWA_GRIDS = ('a', 'b', 'c', 'cx', 'cy')

EXECUTABLE_NAME = 'UPSY_calculate_remapping_weights'

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_WEIGHTS_DIR = REPO_ROOT / 'UPSY_remapping_weights'

# Same name options as in src/UPSY/io/netcdf_basic/netcdf_field_name_options.f90
NAME_OPTIONS_X   = ['x', 'X', 'x1', 'X1', 'nx', 'NX', 'x-coordinate', 'X-coordinate', 'easting', 'Easting']
NAME_OPTIONS_Y   = ['y', 'Y', 'y1', 'Y1', 'ny', 'NY', 'y-coordinate', 'Y-coordinate', 'northing', 'Northing']
NAME_OPTIONS_LON = ['lon', 'Lon', 'long', 'Long', 'longitude', 'Longitude']
NAME_OPTIONS_LAT = ['lat', 'Lat', 'latitude', 'Latitude']

MESH_DIMS = ['vi', 'ti', 'ci', 'two', 'three']
MESH_VARS = ['V', 'nC', 'C', 'niTri', 'iTri', 'VBI', 'Tri', 'Tricc', 'TriC']

# Mesh variables that define the mesh (everything else is derived from these)
MESH_HASH_VARS = ['V', 'Tri', 'xmin', 'xmax', 'ymin', 'ymax', 'lambda_M', 'phi_M', 'beta_stereo']


@dataclass(frozen=True)
class GridInfo:
    """What remap needs to know about a grid or mesh read from a NetCDF file."""
    kind:        str    # 'xy_grid', 'lonlat_grid', or 'mesh'
    hash:        str    # Hash of the variables that define the grid
    shape:       tuple  # Python-order shape of gridded data: (ny,nx) or (nlat,nlon); () for a mesh
    dim_names:   tuple  # Names of those dimensions in the file: ('y','x') or ('lat','lon'); () for a mesh
    mesh_sizes:  tuple  # (nV, nTri, nE) for a mesh; () for a grid

    def n(self, Arakawa_grid):
        """Length of the UPSY data vector on this grid."""
        if self.kind == 'mesh':
            nV, nTri, nE = self.mesh_sizes
            sizes = {'a': nV, 'b': nTri, 'c': nE}
            if Arakawa_grid not in sizes:
                raise ValueError(f"Arakawa grid '{Arakawa_grid}' is not supported on a mesh (only a, b, c)")
            return sizes[Arakawa_grid]
        if Arakawa_grid != 'a':
            raise ValueError(f"Arakawa grid '{Arakawa_grid}' is not supported on an {self.kind} (only a)")
        return int(np.prod(self.shape))


def read_grid_info(filename):
    """Find out what kind of grid/mesh a NetCDF file contains, mirroring the
    inquire_xy_grid / inquire_lonlat_grid / inquire_mesh routines in UPSY."""
    filename = os.path.realpath(filename)
    st = os.stat(filename)
    return _read_grid_info_cached(filename, st.st_mtime_ns, st.st_size)


@lru_cache(maxsize=64)
def _read_grid_info_cached(filename, mtime_ns, size):

    with netCDF4.Dataset(filename) as ds:
        ds.set_auto_mask(False)

        def find(options):
            dim = next((o for o in options if o in ds.dimensions), None)
            var = next((o for o in options if o in ds.variables), None)
            return (dim, var) if dim is not None and var is not None else None

        x, y = find(NAME_OPTIONS_X), find(NAME_OPTIONS_Y)
        lon, lat = find(NAME_OPTIONS_LON), find(NAME_OPTIONS_LAT)
        has_xy_grid     = x is not None and y is not None
        has_lonlat_grid = lon is not None and lat is not None
        has_mesh        = all(d in ds.dimensions for d in MESH_DIMS) and all(v in ds.variables for v in MESH_VARS)

        if has_xy_grid + has_lonlat_grid + has_mesh > 1:
            raise ValueError(f'file {filename} contains more than one of an x/y-grid, a lon/lat-grid, and a mesh')

        if has_xy_grid:
            kind, hash_vars = 'xy_grid', [x[1], y[1]]
            dim_names = (y[0], x[0])
            shape = (len(ds.dimensions[y[0]]), len(ds.dimensions[x[0]]))
            mesh_sizes = ()
        elif has_lonlat_grid:
            kind, hash_vars = 'lonlat_grid', [lon[1], lat[1]]
            dim_names = (lat[0], lon[0])
            shape = (len(ds.dimensions[lat[0]]), len(ds.dimensions[lon[0]]))
            mesh_sizes = ()
        elif has_mesh:
            kind, hash_vars = 'mesh', [v for v in MESH_HASH_VARS if v in ds.variables]
            dim_names, shape = (), ()
            nE = len(ds.dimensions['ei']) if 'ei' in ds.dimensions else None
            mesh_sizes = (len(ds.dimensions['vi']), len(ds.dimensions['ti']), nE)
        else:
            raise ValueError(f'could not find a valid grid or mesh in file {filename}')

        h = hashlib.sha256(kind.encode())
        for name in hash_vars:
            a = np.asarray(ds.variables[name][:])
            a = np.ascontiguousarray(a, dtype=np.int64 if np.issubdtype(a.dtype, np.integer) else np.float64)
            h.update(name.encode())
            h.update(str(a.shape).encode())
            h.update(a.tobytes())

    return GridInfo(kind, h.hexdigest()[:16], shape, dim_names, mesh_sizes)


class Map:
    """A remapping operator between two grids/meshes: d_dst = M @ d_src."""

    def __init__(self, M, grid_src, Arakawa_grid_src, grid_dst, Arakawa_grid_dst, method, attrs=None):
        self.M                = M
        self.grid_src         = grid_src
        self.Arakawa_grid_src = Arakawa_grid_src
        self.grid_dst         = grid_dst
        self.Arakawa_grid_dst = Arakawa_grid_dst
        self.method           = method
        self.attrs            = attrs or {}  # Global attributes of the weights file (program info etc.)
        self.empty_rows       = np.diff(M.indptr) == 0  # Destination cells without any source contributions

    def __repr__(self):
        return (f'Map({self.grid_src.kind} ({self.Arakawa_grid_src}) -> '
                f'{self.grid_dst.kind} ({self.Arakawa_grid_dst}), {self.method}, '
                f'{self.M.shape[0]}x{self.M.shape[1]}, nnz={self.M.nnz})')

    def apply(self, d_src):
        """Apply the map to an array (or xarray.DataArray) of source data.
        Returns a numpy array with the same leading axes as d_src. Destination
        cells without any source contributions are set to NaN."""

        d = _to_numpy_with_grid_dims_last(d_src, self.grid_src)

        if self.grid_src.kind == 'mesh':
            n_src = self.M.shape[1]
            if d.shape[-1] != n_src:
                raise ValueError(f'last axis of source data has length {d.shape[-1]}; expected {n_src}')
            lead = d.shape[:-1]
            d_vec = d.reshape(-1, n_src)
        else:
            if d.shape[-2:] != self.grid_src.shape:
                raise ValueError(f'last two axes of source data have shape {d.shape[-2:]}; '
                                 f'expected {self.grid_src.shape} {self.grid_src.dim_names}')
            lead = d.shape[:-2]
            # UPSY orders grid vectors as n = (i-1)*ny + j (y fastest), i.e. C-order of (x,y)
            d_vec = np.swapaxes(d, -1, -2).reshape(-1, self.M.shape[1])

        d_dst_vec = (self.M @ d_vec.T).T
        if self.empty_rows.any():
            d_dst_vec = d_dst_vec.astype(np.result_type(d_dst_vec, np.float64), copy=False)
            d_dst_vec[:, self.empty_rows] = np.nan

        if self.grid_dst.kind == 'mesh':
            return d_dst_vec.reshape(*lead, self.M.shape[0])
        ny, nx = self.grid_dst.shape
        return np.swapaxes(d_dst_vec.reshape(*lead, nx, ny), -1, -2)


def _to_numpy_with_grid_dims_last(d, grid):
    """For an xarray.DataArray on a grid, make sure the grid dims come last in (y,x) order."""
    if grid.kind != 'mesh' and hasattr(d, 'dims') and hasattr(d, 'transpose'):
        if all(name in d.dims for name in grid.dim_names):
            d = d.transpose(..., *grid.dim_names)
    return np.asarray(d)


class Remapper:
    """Creates, caches, and applies remapping maps.

    weights_dir: directory with pre-calculated weights files
                 (default: $UPSY_REMAPPING_WEIGHTS_DIR, or UPSY_remapping_weights/ in the
                 repository root); created if it does not exist yet
    executable:  path to UPSY_calculate_remapping_weights (default: $UPSY_CALCULATE_REMAPPING_WEIGHTS,
                 then the PATH, then build/src/UPSY/ in this repository)
    n_procs:     number of MPI processes to calculate new weights with
    mpiexec:     MPI launcher command
    timeout:     maximum time (in seconds) for calculating a single set of weights (default: no limit)
    """

    def __init__(self, weights_dir=None, executable=None, n_procs=1, mpiexec='mpiexec', verbose=False,
                 timeout=None):
        self.weights_dir = Path(weights_dir or os.environ.get('UPSY_REMAPPING_WEIGHTS_DIR') or DEFAULT_WEIGHTS_DIR)
        self.weights_dir.mkdir(parents=True, exist_ok=True)
        self.executable  = executable
        self.n_procs     = n_procs
        self.mpiexec     = mpiexec
        self.verbose     = verbose
        self.timeout     = timeout
        self.atlas       = {}

    def remap(self, d_src, filename_grid_src, filename_grid_dst, method,
              Arakawa_grid_src='a', Arakawa_grid_dst='a'):
        """Remap data from the grid/mesh in one NetCDF file to that in another."""
        return self.get_map(filename_grid_src, filename_grid_dst, method,
                            Arakawa_grid_src, Arakawa_grid_dst).apply(d_src)

    def get_map(self, filename_grid_src, filename_grid_dst, method,
                Arakawa_grid_src='a', Arakawa_grid_dst='a', recalculate=False):
        """Get a map from memory, from the weights directory, or by calculating it (in that order)."""

        if method not in METHODS:
            raise ValueError(f"unknown remapping method '{method}'; should be one of {METHODS}")
        for Arakawa_grid in (Arakawa_grid_src, Arakawa_grid_dst):
            if Arakawa_grid not in ARAKAWA_GRIDS:
                raise ValueError(f"unknown Arakawa grid '{Arakawa_grid}'; should be one of {ARAKAWA_GRIDS}")

        grid_src = read_grid_info(filename_grid_src)
        grid_dst = read_grid_info(filename_grid_dst)
        n_src = grid_src.n(Arakawa_grid_src)
        n_dst = grid_dst.n(Arakawa_grid_dst)

        key = (grid_src.kind, grid_src.hash, Arakawa_grid_src, grid_dst.kind, grid_dst.hash, Arakawa_grid_dst, method)

        if not recalculate and key in self.atlas:
            return self.atlas[key]

        filename_weights = self.weights_dir / (
            f'{grid_src.kind}_{grid_src.hash}_{Arakawa_grid_src}__'
            f'{grid_dst.kind}_{grid_dst.hash}_{Arakawa_grid_dst}__{method}.nc')

        M = None
        if not recalculate and filename_weights.exists():
            try:
                M, attrs = read_weights_file(filename_weights, Arakawa_grid_src, Arakawa_grid_dst, method, n_src, n_dst)
                if self.verbose:
                    print(f'Read remapping weights from {filename_weights}')
            except (OSError, KeyError, ValueError) as e:
                warnings.warn(f'invalid remapping weights file {filename_weights} ({e}); recalculating')

        if M is None:
            self._calculate_weights(filename_grid_src, Arakawa_grid_src,
                                    filename_grid_dst, Arakawa_grid_dst, method, filename_weights)
            M, attrs = read_weights_file(filename_weights, Arakawa_grid_src, Arakawa_grid_dst, method, n_src, n_dst)

        self.atlas[key] = Map(M, grid_src, Arakawa_grid_src, grid_dst, Arakawa_grid_dst, method, attrs)
        return self.atlas[key]

    def _calculate_weights(self, filename_grid_src, Arakawa_grid_src,
                           filename_grid_dst, Arakawa_grid_dst, method, filename_weights):
        """Run UPSY_calculate_remapping_weights, writing to a temporary file that is
        then moved into place, so that no incomplete weights files are left behind."""

        filename_tmp = filename_weights.with_name(f'.{filename_weights.stem}.{uuid.uuid4().hex}.tmp.nc')

        cmd = [self.mpiexec, '-n', str(self.n_procs), self._find_executable(),
               os.path.abspath(filename_grid_src), Arakawa_grid_src,
               os.path.abspath(filename_grid_dst), Arakawa_grid_dst,
               method, str(filename_tmp)]
        if self.verbose:
            print('Calculating remapping weights: ' + ' '.join(cmd))

        try:
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=self.timeout)
            except subprocess.TimeoutExpired as e:
                raise RuntimeError(f'{EXECUTABLE_NAME} did not finish within {self.timeout} s: '
                                   + ' '.join(cmd)) from e
            if result.returncode != 0 or not filename_tmp.exists():
                log = _strip_ansi(result.stdout + result.stderr).strip().splitlines()
                errors = [line.strip() for line in log if 'ERROR' in line and 'ERROR STOP' not in line]
                raise RuntimeError(f'{EXECUTABLE_NAME} failed (exit code {result.returncode}): '
                                   + ' '.join(cmd) + '\n' + '\n'.join(errors or log[-20:]))
            os.replace(filename_tmp, filename_weights)
        finally:
            filename_tmp.unlink(missing_ok=True)

    def _find_executable(self):
        candidates = [
            self.executable,
            os.environ.get('UPSY_CALCULATE_REMAPPING_WEIGHTS'),
            shutil.which(EXECUTABLE_NAME),
            REPO_ROOT / 'build' / 'src' / 'UPSY' / EXECUTABLE_NAME]
        for c in candidates:
            if c and os.path.isfile(c) and os.access(c, os.X_OK):
                return str(c)
        raise FileNotFoundError(
            f'could not find {EXECUTABLE_NAME}; pass executable=..., set '
            f'$UPSY_CALCULATE_REMAPPING_WEIGHTS, or add it to the PATH')


def read_weights_file(filename, Arakawa_grid_src, Arakawa_grid_dst, method, n_src, n_dst):
    """Read a weights file written by UPSY_calculate_remapping_weights as a scipy CSR matrix,
    checking that it matches what was asked for. Also returns the file's global attributes."""

    with netCDF4.Dataset(filename) as ds:
        ds.set_auto_mask(False)

        for attr, expected in [('Arakawa_grid_src', Arakawa_grid_src),
                               ('Arakawa_grid_dst', Arakawa_grid_dst),
                               ('method', method)]:
            if ds.getncattr(attr) != expected:
                raise ValueError(f"{attr} = '{ds.getncattr(attr)}', expected '{expected}'")

        m, n = len(ds.dimensions['m']), len(ds.dimensions['n'])
        if (m, n) != (n_dst, n_src):
            raise ValueError(f'matrix is {m}x{n}, expected {n_dst}x{n_src}')

        # Fortran CSR: 1-based row pointers and column indices
        ptr = ds.variables['ptr'][:].astype(np.int64) - 1
        ind = ds.variables['ind'][:].astype(np.int64) - 1
        val = ds.variables['val'][:]
        attrs = {name: ds.getncattr(name) for name in ds.ncattrs()}

    if len(ptr) != m + 1 or ptr[0] != 0 or ptr[-1] != len(ind):
        raise ValueError('inconsistent CSR row pointers')

    return sp.csr_matrix((val, ind, ptr), shape=(m, n)), attrs


def _strip_ansi(s):
    return re.sub(r'\x1b\[[0-9;]*m', '', s)


_default_remapper = None


def get_remapper():
    """The Remapper used by remap() and get_map()."""
    global _default_remapper
    if _default_remapper is None:
        _default_remapper = Remapper()
    return _default_remapper


def remap(d_src, filename_grid_src, filename_grid_dst, method,
          Arakawa_grid_src='a', Arakawa_grid_dst='a'):
    """Remap data from the grid/mesh in one NetCDF file to that in another,
    using (and caching) remapping weights as needed. See the module docstring."""
    return get_remapper().remap(d_src, filename_grid_src, filename_grid_dst, method,
                                Arakawa_grid_src, Arakawa_grid_dst)


def get_map(filename_grid_src, filename_grid_dst, method,
            Arakawa_grid_src='a', Arakawa_grid_dst='a', recalculate=False):
    """Get a Map object (from memory, from the weights directory, or newly calculated)."""
    return get_remapper().get_map(filename_grid_src, filename_grid_dst, method,
                                  Arakawa_grid_src, Arakawa_grid_dst, recalculate)
