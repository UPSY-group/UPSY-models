"""Component test for the Python remapping wrapper (tools/python/upsy/remap.py) and
the Fortran program UPSY_calculate_remapping_weights behind it.

Remaps an analytical test field from every x/y-grid to every x/y-grid in the
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


def read_grid(filename):
    with netCDF4.Dataset(filename) as ds:
        return ds['x'][:].data, ds['y'][:].data


def write_result(filename, d, dims, attrs):
    """Write a remapped field to a NetCDF file. The global attributes of the
    weights file (program info) are copied along, so the test suite can
    check which build of UPSY_calculate_remapping_weights produced them."""
    with netCDF4.Dataset(filename, 'w') as ds:
        ds.setncatts(attrs)
        for name, size in zip(dims, d.shape):
            ds.createDimension(name, size)
        ds.createVariable('d', 'f8', dims)[:] = d


def plot_comparison(filename, title, x_src, y_src, x_dst, y_dst, d_grid_src, d_grid_dst):
    """1x2 figure: source -> destination."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    cmap = plt.get_cmap('viridis')
    # nanmin/nanmax: destination cells without any overlapping source cells are NaN
    norm = plt.Normalize(min(np.nanmin(d_grid_src), np.nanmin(d_grid_dst)), max(np.nanmax(d_grid_src), np.nanmax(d_grid_dst)))

    def plot_grid(ax, x, y, d):
        return ax.pcolormesh(x, y, d, shading='nearest', cmap=cmap, norm=norm)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
    plot_grid(axs[0], x_src, y_src, d_grid_src)
    pc = plot_grid(axs[1], x_dst, y_dst, d_grid_dst)

    axs[0].set_title('source')
    axs[1].set_title('destination')

    dx_src, dy_src = x_src[1] - x_src[0], y_src[1] - y_src[0]
    dx_dst, dy_dst = x_dst[1] - x_dst[0], y_dst[1] - y_dst[0]
    xlim = (min(x_src[0] - dx_src / 2, x_dst[0] - dx_dst / 2), max(x_src[-1] + dx_src / 2, x_dst[-1] + dx_dst / 2))
    ylim = (min(y_src[0] - dy_src / 2, y_dst[0] - dy_dst / 2), max(y_src[-1] + dy_src / 2, y_dst[-1] + dy_dst / 2))
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

    grids  = sorted(Path(args.test_meshes_and_grids_dir).glob('grid_*.nc'))
    if not grids:
        raise FileNotFoundError(f'no grids found in {args.test_meshes_and_grids_dir}')

    os.makedirs(args.results_dir, exist_ok=True)
    if args.figures:
        os.makedirs(args.figures, exist_ok=True)

    remapper = Remapper(weights_dir=args.weights_dir, n_procs=N_PROCS, timeout=TIMEOUT)

    n_total = len(grids) * len(grids) * len(METHODS)
    i = 0

    for grid_src in grids:

        x_src, y_src = read_grid(grid_src)
        d_grid_src = test_field(x_src[None, :], y_src[:, None])

        for grid_dst in grids:

            x_dst, y_dst = read_grid(grid_dst)

            for method in METHODS:
                i += 1
                print(f'[{i:3d}/{n_total}] {grid_src.stem} -> {grid_dst.stem}, {method}',
                      end='', flush=True)
                t0 = time.time()

                map_g2g = remapper.get_map(grid_src, grid_dst, method, 'a', 'a')
                d_grid_dst = map_g2g.apply(d_grid_src)
                write_result(
                    Path(args.results_dir) / f'res_{grid_src.stem}_TO_{grid_dst.stem}_{method}.nc',
                    d_grid_dst, ('y', 'x'), map_g2g.attrs)

                if args.figures:
                    plot_comparison(
                        Path(args.figures) / f'{grid_src.stem}_TO_{grid_dst.stem}_{method}.png',
                        f'{grid_src.name}\n{grid_dst.name}\n{method}',
                        x_src, y_src, x_dst, y_dst, d_grid_src, d_grid_dst)

                print(f': {time.time() - t0:.1f} s', flush=True)


if __name__ == '__main__':
    main()
