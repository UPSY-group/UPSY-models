#!/usr/bin/env python3
"""Create Halfar dome benchmark plots for the UFEMISM weekly benchmarks."""

import argparse
import glob
import os
import re

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.lines import Line2D


HALFAR_SIMULATIONS = [
    ('40 km', 'results_Halfar_40km', 'tab:blue'),
    ('20 km', 'results_Halfar_20km', 'tab:orange'),
    ('10 km', 'results_Halfar_10km', 'tab:green'),
    ('5 km', 'results_Halfar_5km', 'tab:red'),
]

HALFAR_STATIC_SIMULATIONS = [
    ('40 km', 'results_Halfar_static_40km', 'tab:blue'),
    ('20 km', 'results_Halfar_static_20km', 'tab:orange'),
    ('10 km', 'results_Halfar_static_10km', 'tab:green'),
    ('5 km', 'results_Halfar_static_5km', 'tab:red'),
]

SECONDS_PER_YEAR = 31556943.36


def parse_config_number(config_file, key):
    """Read one numeric value from a UFEMISM config file."""
    key_pattern = re.compile(r'^\s*' + re.escape(key) + r'\s*=\s*(.*?)\s*(?:!.*)?$')

    with open(config_file, 'r', encoding='utf-8') as f:
        for raw_line in f:
            line = raw_line.rstrip('\n')
            match = key_pattern.match(line)
            if not match:
                continue

            value_str = match.group(1).strip().strip("'").strip('"')
            value_str = value_str.replace('D', 'E').replace('d', 'e')
            return float(value_str)

    raise KeyError(f'Key {key} not found in {config_file}')


def read_halfar_parameters(config_file):
    """Read Halfar analytical-solution parameters from one config file."""
    return {
        'A_flow': parse_config_number(config_file, 'uniform_Glens_flow_factor_config'),
        'n_flow': parse_config_number(config_file, 'Glens_flow_law_exponent_config'),
        'H0': parse_config_number(config_file, 'refgeo_idealised_Halfar_H0_config'),
        'R0': parse_config_number(config_file, 'refgeo_idealised_Halfar_R0_config'),
        'time': parse_config_number(config_file, 'end_time_of_run_config'),
    }


def halfar_solution(A_flow, n_flow, H0, R0, x, y, time):
    """Compute the Halfar analytical ice-thickness solution."""
    ice_density = 917.0
    grav = 9.81

    gamma = (2.0 / 5.0) * (A_flow / SECONDS_PER_YEAR) * (ice_density * grav) ** n_flow
    t0 = (
        1.0
        / ((5.0 * n_flow + 3.0) * gamma)
        * ((2.0 * n_flow + 1.0) / (n_flow + 1.0)) ** n_flow
        * (R0 ** (n_flow + 1.0))
        / (H0 ** (2.0 * n_flow + 1.0))
    )

    tp = (time * SECONDS_PER_YEAR) + t0
    r = np.sqrt(np.asarray(x, dtype=float) ** 2 + np.asarray(y, dtype=float) ** 2)

    f1 = (t0 / tp) ** (2.0 / (5.0 * n_flow + 3.0))
    f2 = (t0 / tp) ** (1.0 / (5.0 * n_flow + 3.0))
    f3 = r / R0

    inside = 1.0 - (f2 * f3) ** ((n_flow + 1.0) / n_flow)
    H = H0 * f1 * np.maximum(0.0, inside) ** (n_flow / (2.0 * n_flow + 1.0))
    return H


def list_transect_files(results_folder):
    """Return all transect files in one UFEMISM results folder."""
    pattern = os.path.join(results_folder, 'transect_*.nc')
    transect_files = sorted(glob.glob(pattern))
    if not transect_files:
        raise FileNotFoundError(f'No transect files found in {results_folder}')
    return transect_files


def read_final_thickness(transect_file):
    """Read final-step ice thickness and radial distance from one transect file."""
    with netCDF4.Dataset(transect_file, 'r') as ds:
        if 'Hi' not in ds.variables:
            raise KeyError(f'Hi not found in {transect_file}')
        if 'V' not in ds.variables:
            raise KeyError(f'V not found in {transect_file}')

        time = ds['time'][:]
        if len(time) == 0:
            raise ValueError(f'Empty time axis in {transect_file}')

        ti_end = len(time) - 1
        thickness = np.asarray(ds['Hi'][ti_end, :], dtype=float)

        vertices = np.asarray(ds['V'][:], dtype=float)
        if vertices.shape[0] != 2:
            raise ValueError(f'Unexpected V shape in {transect_file}: {vertices.shape}')

        x = vertices[0, :]
        y = vertices[1, :]
        radius_km = np.sqrt(x**2 + y**2) / 1e3

    order = np.argsort(radius_km)
    return radius_km[order], thickness[order]


def plot_group(ax, results_root, simulations, panel_title, analytical_time_mode='end'):
    """Plot all final-thickness transects for one simulation group on one axis."""
    x_analytic_m = np.linspace(0.0, 600e3, 1201)
    analytical_plotted = False

    for label, folder, color in simulations:
        folder_path = os.path.join(results_root, folder)
        transect_files = list_transect_files(folder_path)

        for transect_file in transect_files:
            radius_km, thickness = read_final_thickness(transect_file)
            ax.plot(radius_km, thickness, color=color, linewidth=1.1, alpha=0.65)

        config_name = folder.replace('results_', 'config_') + '.cfg'
        config_file = os.path.join(results_root, config_name)
        if not os.path.exists(config_file):
            raise FileNotFoundError(f'Config file not found for {folder}: {config_file}')

        if not analytical_plotted:
            params = read_halfar_parameters(config_file)
            if analytical_time_mode == 'zero':
                analytical_time = 0.0
            else:
                analytical_time = params['time']

            H_analytic = halfar_solution(
                A_flow=params['A_flow'],
                n_flow=params['n_flow'],
                H0=params['H0'],
                R0=params['R0'],
                x=x_analytic_m,
                y=0.0,
                time=analytical_time,
            )
            ax.plot(x_analytic_m / 1e3, H_analytic, color='k', linestyle='--', linewidth=2.0)
            analytical_plotted = True

    ax.set_title(panel_title, fontsize=12)
    ax.set_xlabel('Distance from dome centre (km)', fontsize=11)
    ax.grid(True, alpha=0.3)
    legend_handles = [
        Line2D([], [], color='k', linestyle='-', linewidth=2.0, label='UFEMISM transects'),
        Line2D([], [], color='k', linestyle='--', linewidth=2.0, label='Analytical'),
    ]
    for label, _, color in simulations:
        legend_handles.append(Line2D([], [], color=color, linestyle='-', linewidth=2.0, label=label))

    ax.legend(handles=legend_handles, loc='best', fontsize=9, frameon=False)


def main():
    parser = argparse.ArgumentParser(
        description='Create Halfar dome full benchmark figure from UFEMISM results.'
    )
    parser.add_argument(
        '--results-root',
        default='.',
        help='Path to the folder containing UFEMISM results_* directories.',
    )
    parser.add_argument(
        '--output-dir',
        default=None,
        help='Directory where output figure PNG files are written.',
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

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), dpi=100, sharex=True, sharey=True)

    plot_group(
        ax=axes[0],
        results_root=results_root,
        simulations=HALFAR_SIMULATIONS,
        panel_title='Halfar',
        analytical_time_mode='end',
    )
    axes[0].set_ylabel('Final ice thickness (m)', fontsize=11)

    plot_group(
        ax=axes[1],
        results_root=results_root,
        simulations=HALFAR_STATIC_SIMULATIONS,
        panel_title='Halfar_static',
        analytical_time_mode='zero',
    )

    axes[0].set_xlim(0.0, 600.0)

    output_file = os.path.join(output_dir, 'Fig_integrated_test_Halfar_dome_full.png')
    fig.tight_layout()
    fig.savefig(output_file, dpi=100)
    print(f'Figure saved to {output_file}')
    plt.close(fig)


if __name__ == '__main__':
    main()
