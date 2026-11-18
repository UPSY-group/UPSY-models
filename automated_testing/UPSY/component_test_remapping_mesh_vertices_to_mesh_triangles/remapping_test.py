"""Component test for the Python remapping wrapper (tools/python/upsy/remap.py) and
the Fortran program UPSY_calculate_remapping_weights behind it.

Remaps an analytical test field from every mesh (vertices) to every mesh (triangles) in the
test_meshes_and_grids directory, for every supported remapping method.

Usage:
    python3 remapping_test.py [test_meshes_and_grids_dir] [results_dir] [weights_dir] [--figures FIGURES_DIR]

If --figures is given, a figure comparing the source and destination fields is
saved for every gcombination (for visual inspection only; the figures are not checked).
"""

import argparse
import os
import sys
import time
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'tools' / 'python'))
from upsy.remap import Remapper

# Mesh Arakawa grids to test (c-grid is not yet supported)
MESH_ARAKAWA_GRIDS = {'a': 'vertices', 'b': 'triangles'}

# Remapping methods to test (linear and second-order conservative are not yet implemented)
METHODS = ['nearest_neighbour', 'first_order_conservative']

# Number of processes to calculate weights with, and maximum time per set of weights
N_PROCS = 2
TIMEOUT = 600


def test_field(x, y):
    """A smooth field with some structure at different scales (x, y in m)."""
    return (np.cos(2 * np.pi * x / 1.5e6) * np.sin(2 * np.pi * y / 2.0e6)
            + 0.5 * np.exp(-((x - 1.0e6)**2 + (y + 0.5e6)**2) / (4.0e5)**2)
            + x / 3.0e6)


def read_mesh(filename, Arakawa_grid):
    """Coordinates of the mesh vertices (a) or triangle geometric centres (b),
    and the polygons of their Voronoi cells/triangles (for plotting)."""
    with netCDF4.Dataset(filename) as ds:
        if Arakawa_grid == 'a':
            xy    = ds['V'    ][:].data
            Vor   = ds['Vor'  ][:].data
            nVVor = ds['nVVor'][:].data
            VVor  = ds['VVor' ][:].data - 1
            polygons = [Vor[:, VVor[:nVVor[vi], vi]].T for vi in range(xy.shape[1])]
        elif Arakawa_grid == 'b':
            xy  = ds['TriGC'][:].data
            V   = ds['V'    ][:].data
            Tri = ds['Tri'  ][:].data - 1
            polygons = list(np.transpose(V[:, Tri], (2, 1, 0)))
        else:
            raise ValueError(f"unsupported Arakawa grid '{Arakawa_grid}'")
    return xy[0], xy[1], polygons


def write_result(filename, d, dims, attrs):
    """Write a remapped field to a NetCDF file. The global attributes of the
    weights file (program info) are copied along, so the test suite can
    check which build of UPSY_calculate_remapping_weights produced them."""
    with netCDF4.Dataset(filename, 'w') as ds:
        ds.setncatts(attrs)
        for name, size in zip(dims, d.shape):
            ds.createDimension(name, size)
        ds.createVariable('d', 'f8', dims)[:] = d


def plot_comparison(filename, title, xm_src, ym_src, polygons_src, xm_dst, ym_dst, polygons_dst, d_mesh_src, d_mesh_dst):
    """1x2 figure: source -> destination."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    cmap = plt.get_cmap('viridis')
    # nanmin/nanmax: destination cells without any overlapping source cells are NaN
    norm = plt.Normalize(min(np.nanmin(d_mesh_dst), np.nanmin(d_mesh_src)), max(np.nanmax(d_mesh_dst), np.nanmax(d_mesh_src)))

    def plot_mesh(ax, polygons, d):
        return ax.add_collection(PolyCollection(polygons, array=d, cmap=cmap, norm=norm,
                                                edgecolors='face', linewidths=0.2))

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
    plot_mesh(axs[0], polygons_src, d_mesh_src)
    pc = plot_mesh(axs[1], polygons_dst, d_mesh_dst)

    axs[0].set_title('source')
    axs[1].set_title('destination')

    xlim = (min(xm_src.min(), xm_dst.min()), max(xm_src.max(), xm_dst.max()))
    ylim = (min(ym_src.min(), ym_dst.min()), max(ym_src.max(), ym_dst.max()))
    for ax in axs.flat:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(pc, ax=axs, shrink=0.6, label='test field')
    fig.suptitle(title)
    fig.savefig(filename, dpi=100)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('test_meshes_and_grids_dir')
    parser.add_argument('results_dir')
    parser.add_argument('weights_dir')
    parser.add_argument('--figures', metavar='FIGURES_DIR')
    args = parser.parse_args()

    meshes = sorted(Path(args.test_meshes_and_grids_dir).glob('mesh_*.nc'))
    if not meshes:
        raise FileNotFoundError(f'no meshes found in {args.test_meshes_and_grids_dir}')

    os.makedirs(args.results_dir, exist_ok=True)
    if args.figures:
        os.makedirs(args.figures, exist_ok=True)

    remapper = Remapper(weights_dir=args.weights_dir, n_procs=N_PROCS, timeout=TIMEOUT)

    n_total = len(meshes) * len(meshes) * len(METHODS)
    i = 0

    Arakawa_grid_src = 'a'
    mesh_src_location = 'vertices'

    Arakawa_grid_dst = 'b'
    mesh_dst_location = 'triangles'

    for mesh_src in meshes:

        xm_src, ym_src, polygons_src = read_mesh(mesh_src, Arakawa_grid_src)
        d_mesh_src = test_field(xm_src, ym_src)
        mesh_src_dim = 'vi' if Arakawa_grid_src == 'a' else 'ti'

        for mesh_dst in meshes:

            xm_dst, ym_dst, polygons_dst = read_mesh(mesh_dst, Arakawa_grid_dst)
            d_mesh_dst = test_field(xm_dst, ym_dst)
            mesh_dst_dim = 'vi' if Arakawa_grid_dst == 'a' else 'ti'

            for method in METHODS:
                i += 1
                print(f'[{i:3d}/{n_total}] {mesh_src.stem} ({mesh_src_location}) -> {mesh_dst.stem} ({mesh_dst_location}), {method}',
                      end='', flush=True)
                t0 = time.time()

                map_m2m = remapper.get_map(mesh_src, mesh_dst, method, Arakawa_grid_src, Arakawa_grid_dst)
                d_mesh_dst = map_m2m.apply(d_mesh_src)
                write_result(
                    Path(args.results_dir) / f'res_{mesh_src.stem}_{mesh_src_location}_TO_{mesh_dst.stem}_{mesh_dst_location}_{method}.nc',
                    d_mesh_dst, (mesh_dst_dim,), map_m2m.attrs)

                if args.figures:
                    plot_comparison(
                        Path(args.figures) / f'{mesh_src.stem}_{mesh_src_location}_TO_{mesh_dst.stem}_{mesh_dst_location}_{method}.png',
                        f'{mesh_src.name} ({mesh_src_location})\n{mesh_dst.name} ({mesh_dst_location})\n{method}',
                        xm_src, ym_src, polygons_src, xm_dst, ym_dst, polygons_dst, d_mesh_src, d_mesh_dst)

                print(f': {time.time() - t0:.1f} s', flush=True)


if __name__ == '__main__':
    main()
