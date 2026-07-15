#!/usr/bin/env python3
"""Create a velocity comparison plot for the SSA icestream full integrated test."""

import argparse
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4
import numpy as np


SIMULATIONS = [
    ('32 km', 'results_32km', 'tab:blue'),
    ('16 km', 'results_16km', 'tab:orange'),
    ('8 km', 'results_8km', 'tab:green'),
    ('4 km', 'results_4km', 'tab:red'),
]


def Schoof2006_icestream(A, H, tantheta, L, m, y):
    """Reproduce the Schoof (2006) ice-stream analytical solution."""

    ice_density = 910.0
    grav = 9.81

    y = np.asarray(y, dtype=float)

    f = -ice_density * grav * H * tantheta
    B = A ** (-1.0 / 3.0)
    W = L * (m + 1.0) ** (1.0 / m)

    ua = -2.0 * f**3 * L**4 / (B**3 * H**3)
    ub = 0.25 * ((y / L) ** 4 - (m + 1.0) ** (4.0 / m))
    uc = (-3.0 / ((m + 1.0) * (m + 4.0))) * (
        np.abs(y / L) ** (m + 4.0) - (m + 1.0) ** (1.0 + (4.0 / m))
    )
    ud = (3.0 / ((m + 1.0) ** 2 * (2.0 * m + 4.0))) * (
        np.abs(y / L) ** (2.0 * m + 4.0) - (m + 1.0) ** (2.0 + (4.0 / m))
    )
    ue = (-1.0 / ((m + 1.0) ** 3 * (3.0 * m + 4.0))) * (
        np.abs(y / L) ** (3.0 * m + 4.0) - (m + 1.0) ** (3.0 + (4.0 / m))
    )

    u = ua * (ub + uc + ud + ue)
    u = np.where(np.abs(y) > W, 0.0, u)

    return u


def read_surface_velocity(filename):
    """Read the transect coordinate and surface velocity from one NetCDF file."""

    if not os.path.exists(filename):
        raise FileNotFoundError(f'Required UFEMISM output not found: {filename}')

    with netCDF4.Dataset(filename, 'r') as ds:
        v = ds['V'][:]
        if v.ndim != 2 or v.shape[1] < 2:
            raise ValueError(f'Unexpected V shape in {filename}: {v.shape}')

        distance_km = np.asarray(v[1, :], dtype=float) / 1e3

        if 'u_3D' not in ds.variables:
            raise KeyError(f'u_3D not found in {filename}')

        u_3d = ds['u_3D'][:]

        if u_3d.ndim != 3:
            raise ValueError(f'Unexpected velocity shape in {filename}: {u_3d.shape}')

        surface_velocity = np.asarray(u_3d[-1, 0, :], dtype=float)

    return distance_km, surface_velocity


def main():
    parser = argparse.ArgumentParser(
        description='Create the SSA icestream full velocity comparison figure.'
    )
    parser.add_argument(
        '--results-root',
        default='.',
        help='Path to the folder containing UFEMISM results_* directories.',
    )
    parser.add_argument(
        '--output-dir',
        default=None,
        help='Directory where the output figure is written.',
    )
    args = parser.parse_args()

    results_root = os.path.abspath(args.results_root)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, '../../..'))

    if args.output_dir is None:
        output_dir = os.path.join(repo_root, 'automated_testing', 'figures')
    else:
        output_dir = os.path.abspath(args.output_dir)

    os.makedirs(output_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.0, 4.2), dpi=100)
    ax.grid(True, alpha=0.3)

    x_limits = []
    y_limits = []
    for label, folder, color in SIMULATIONS:
        filename = os.path.join(results_root, folder, 'transect_southnorth.nc')
        distance_km, surface_velocity = read_surface_velocity(filename)

        order = np.argsort(distance_km)
        distance_km = distance_km[order]
        surface_velocity = surface_velocity[order]

        ax.plot(distance_km, surface_velocity, color=color, linewidth=2.0, label=label)

        x_limits.extend([np.nanmin(distance_km), np.nanmax(distance_km)])
        y_limits.extend([np.nanmin(surface_velocity), np.nanmax(surface_velocity)])

    analytic_x_m = np.linspace(min(x_limits) * 1e3, max(x_limits) * 1e3, 1000)
    analytic_velocity = Schoof2006_icestream(1e-18, 2000.0, -3e-4, 150e3, 1.0, analytic_x_m)
    ax.plot(
        analytic_x_m / 1e3,
        analytic_velocity,
        color='k',
        linestyle='--',
        linewidth=2.0,
        label='Analytical solution',
    )

    y_limits.extend([np.nanmin(analytic_velocity), np.nanmax(analytic_velocity)])

    ax.set_xlabel('Distance x (km)', fontsize=11)
    ax.set_ylabel('Surface velocity u (m/yr)', fontsize=11)
    ax.set_xlim(min(x_limits), max(x_limits))
    ax.set_ylim(min(y_limits), max(y_limits))
    ax.legend(loc='best', fontsize=10, frameon=False)

    output_file = os.path.join(output_dir, 'Fig_integrated_test_SSA_icestream_full_velocity.png')
    fig.savefig(output_file, dpi=100, bbox_inches='tight')
    print(f'Figure saved to {output_file}')
    plt.close(fig)


if __name__ == '__main__':
    main()