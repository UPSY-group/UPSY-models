#!/usr/bin/env python3
"""Create figure(s) for the Ant_init_full_01 UFEMISM integrated test."""

import argparse
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


def resolve_results_dir(path):
    """Accept either a results folder or a parent folder containing results/."""

    candidate_file = os.path.join(path, 'main_output_ANT_00001.nc')
    if os.path.isfile(candidate_file):
        return path

    nested_results_dir = os.path.join(path, 'results')
    nested_candidate_file = os.path.join(nested_results_dir, 'main_output_ANT_00001.nc')
    if os.path.isfile(nested_candidate_file):
        return nested_results_dir

    return path


def build_voronoi_polygons_km(ds):
    """Construct Voronoi-cell polygons in km from UFEMISM mesh connectivity."""

    n_vvor = np.asarray(ds['nVVor'][:], dtype=int)
    vvor = np.asarray(ds['VVor'][:], dtype=int)
    vor = np.asarray(ds['Vor'][:], dtype=float)

    polygons = []
    for vi in range(len(n_vvor)):
        n_vertices = n_vvor[vi]
        # Connectivity is 1-based in UFEMISM output.
        vor_indices = vvor[:n_vertices, vi] - 1
        x_km = vor[0, vor_indices] / 1e3
        y_km = vor[1, vor_indices] / 1e3
        polygons.append(np.column_stack((x_km, y_km)))

    return polygons


def plot_hs_panel(ax, dataset):
    """Plot final-step Hs on the unstructured Voronoi mesh."""

    if 'Hs' not in dataset.variables:
        raise KeyError('Hs not found in NetCDF file')

    hs = np.asarray(dataset['Hs'][:], dtype=float)
    if hs.ndim != 2:
        raise ValueError(f'Unexpected Hs shape: {hs.shape}')

    hs_final = hs[-1, :]
    polygons = build_voronoi_polygons_km(dataset)

    finite_values = hs_final[np.isfinite(hs_final)]
    if finite_values.size == 0:
        raise ValueError('Hs contains no finite values at the final time step')

    poly = PolyCollection(
        polygons,
        array=hs_final,
        cmap='viridis',
        edgecolors='none',
        linewidths=0.0,
    )
    poly.set_clim(vmin=np.nanmin(finite_values), vmax=np.nanmax(finite_values))
    ax.add_collection(poly)

    xmin = float(np.asarray(dataset['xmin'][:])) / 1e3
    xmax = float(np.asarray(dataset['xmax'][:])) / 1e3
    ymin = float(np.asarray(dataset['ymin'][:])) / 1e3
    ymax = float(np.asarray(dataset['ymax'][:])) / 1e3

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Final ice surface elevation (Hs)', fontsize=11)
    ax.grid(False)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('left', size='4.5%', pad=0.08)
    cbar = plt.colorbar(poly, cax=cax)
    cbar.set_label('Hs (m)', fontsize=10)
    cax.yaxis.set_ticks_position('left')
    cax.yaxis.set_label_position('left')


def plot_dhi_panel(ax, dataset):
    """Plot final-step dHi on the unstructured Voronoi mesh."""

    if 'dHi' not in dataset.variables:
        raise KeyError('dHi not found in NetCDF file')

    dhi = np.asarray(dataset['dHi'][:], dtype=float)
    if dhi.ndim != 2:
        raise ValueError(f'Unexpected dHi shape: {dhi.shape}')

    dhi_final = dhi[-1, :]
    polygons = build_voronoi_polygons_km(dataset)

    poly = PolyCollection(
        polygons,
        array=dhi_final,
        cmap='RdBu_r',
        edgecolors='none',
        linewidths=0.0,
    )
    poly.set_clim(vmin=-300.0, vmax=300.0)
    ax.add_collection(poly)

    xmin = float(np.asarray(dataset['xmin'][:])) / 1e3
    xmax = float(np.asarray(dataset['xmax'][:])) / 1e3
    ymin = float(np.asarray(dataset['ymin'][:])) / 1e3
    ymax = float(np.asarray(dataset['ymax'][:])) / 1e3

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Final ice thickness difference (dHi)', fontsize=11)
    ax.grid(False)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='4.5%', pad=0.08)
    cbar = plt.colorbar(poly, cax=cax)
    cbar.set_label('dHi (m)', fontsize=10)


def plot_uabs_surf_panel(ax, dataset):
    """Plot final-step uabs_surf on the unstructured triangle mesh."""

    if 'uabs_surf' not in dataset.variables:
        raise KeyError('uabs_surf not found in NetCDF file')

    uabs_surf = np.asarray(dataset['uabs_surf'][:], dtype=float)
    if uabs_surf.ndim != 2:
        raise ValueError(f'Unexpected uabs_surf shape: {uabs_surf.shape}')

    uabs_final = uabs_surf[-1, :]

    if 'Hi' not in dataset.variables:
        raise KeyError('Hi not found in NetCDF file')
    hi = np.asarray(dataset['Hi'][:], dtype=float)
    if hi.ndim != 2:
        raise ValueError(f'Unexpected Hi shape: {hi.shape}')
    hi_final = hi[-1, :]

    tri = np.asarray(dataset['Tri'][:], dtype=int) - 1
    vertices = np.asarray(dataset['V'][:], dtype=float) / 1e3
    triangles = [
        np.column_stack((vertices[0, tri[:, ti]], vertices[1, tri[:, ti]]))
        for ti in range(tri.shape[1])
    ]

    # Keep triangles where at least one corner vertex has substantial ice.
    hi_on_tri_vertices = hi_final[tri]
    is_ice_covered_triangle = np.any(hi_on_tri_vertices > 0.1, axis=0)

    # Log normalization requires strictly positive values.
    uabs_plot = np.where(
        is_ice_covered_triangle & np.isfinite(uabs_final) & (uabs_final > 0.0),
        uabs_final,
        np.nan,
    )
    finite_values = uabs_plot[np.isfinite(uabs_plot)]
    if finite_values.size == 0:
        raise ValueError('uabs_surf contains no positive finite values at the final time step')

    poly = PolyCollection(
        triangles,
        array=uabs_plot,
        cmap='plasma',
        norm=LogNorm(vmin=np.nanmin(finite_values), vmax=np.nanmax(finite_values)),
        edgecolors='none',
        linewidths=0.0,
    )
    ax.add_collection(poly)

    xmin = float(np.asarray(dataset['xmin'][:])) / 1e3
    xmax = float(np.asarray(dataset['xmax'][:])) / 1e3
    ymin = float(np.asarray(dataset['ymin'][:])) / 1e3
    ymax = float(np.asarray(dataset['ymax'][:])) / 1e3

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Final surface velocity (uabs_surf)', fontsize=11)
    ax.grid(False)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('left', size='4.5%', pad=0.08)
    cbar = plt.colorbar(poly, cax=cax)
    cbar.set_label('uabs_surf (m yr$^{-1}$)', fontsize=10)
    cax.yaxis.set_ticks_position('left')
    cax.yaxis.set_label_position('left')


def plot_ice_volume_af_panel(ax, dataset):
    """Plot the time series of volume above floatation (ice_volume_af)."""

    if 'time' not in dataset.variables:
        raise KeyError('time not found in scalar NetCDF file')
    if 'ice_volume_af' not in dataset.variables:
        raise KeyError('ice_volume_af not found in scalar NetCDF file')

    time = np.asarray(dataset['time'][:], dtype=float)
    ice_volume_af = np.asarray(dataset['ice_volume_af'][:], dtype=float)

    if time.ndim != 1:
        raise ValueError(f'Unexpected time shape: {time.shape}')
    if ice_volume_af.ndim != 1:
        raise ValueError(f'Unexpected ice_volume_af shape: {ice_volume_af.shape}')
    if time.shape[0] != ice_volume_af.shape[0]:
        raise ValueError('time and ice_volume_af must have the same length')

    valid = np.isfinite(time) & np.isfinite(ice_volume_af) & (np.abs(ice_volume_af) < 1e30)
    if not np.any(valid):
        raise ValueError('ice_volume_af contains no valid finite values')

    time_valid = time[valid]
    vaf_valid = ice_volume_af[valid]

    ax.plot(time_valid, vaf_valid, color='k', linewidth=1.5)
    ax.set_title('Volume above floatation (ice_volume_af)', fontsize=11)
    ax.set_xlabel('Time (yr)', fontsize=10)
    ax.set_ylabel('VAF (m s.l.e.)', fontsize=10)
    ax.grid(True, alpha=0.3)


def main():
    """Create a 2x2 diagnostic figure for integrated_test_Ant_init_full_01."""

    parser = argparse.ArgumentParser(
        description='Create Ant_init_full_01 figures from UFEMISM results.'
    )
    parser.add_argument(
        '--results-root',
        default='.',
        help='Path to the folder containing UFEMISM output or its results/ subfolder.',
    )
    parser.add_argument(
        '--output-dir',
        default=None,
        help='Directory where output figure PNG files are written.',
    )
    args = parser.parse_args()

    results_root = resolve_results_dir(os.path.abspath(args.results_root))

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, '../../..'))

    if args.output_dir is None:
        output_dir = os.path.join(repo_root, 'automated_testing', 'figures')
    else:
        output_dir = os.path.abspath(args.output_dir)

    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, 'Fig_integrated_test_Ant_init_full_01.png')
    nc_file = os.path.join(results_root, 'main_output_ANT_00001.nc')
    scalar_nc_file = os.path.join(results_root, 'scalar_output_ANT_00001.nc')

    if not os.path.isfile(nc_file):
        raise FileNotFoundError(f'Required UFEMISM output not found: {nc_file}')
    if not os.path.isfile(scalar_nc_file):
        raise FileNotFoundError(f'Required UFEMISM scalar output not found: {scalar_nc_file}')

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 9.0), dpi=100)

    with netCDF4.Dataset(nc_file, 'r') as ds:
        plot_hs_panel(axes[0, 0], ds)
        plot_dhi_panel(axes[0, 1], ds)
        plot_uabs_surf_panel(axes[1, 0], ds)

    with netCDF4.Dataset(scalar_nc_file, 'r') as ds_scalar:
        plot_ice_volume_af_panel(axes[1, 1], ds_scalar)

    fig.tight_layout()
    fig.savefig(output_file, dpi=100)
    print(f'Figure saved to {output_file}')
    plt.close(fig)


if __name__ == '__main__':
    main()
