#!/usr/bin/env python3
"""
Plot UFEMISM model output for the MISMIP_mod experiment.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4
import os
import glob
import argparse


def list_transect_files(foldername):
    """Return all transect NetCDF files in a results folder."""
    pattern = os.path.join(foldername, 'transect_*.nc')
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No transect files found in: {foldername}")
    return files


def read_transect_results(filename):
    """Read one UFEMISM transect file and return end-state geometry and xGL(t)."""
    with netCDF4.Dataset(filename, 'r') as ds:
        time = ds['time'][:]
        xgl = ds['grounding_line_distance_from_start'][:]

        V = ds['V'][:]
        x = V[0, :]
        y = V[1, :]
        d = np.sqrt(x**2 + y**2)

        ti_end = len(time) - 1
        hs_end = ds['Hs'][ti_end, :]
        hib_end = ds['Hib'][ti_end, :]
        hb_end = ds['Hb'][ti_end, :]

    return {
        'transect_name': os.path.splitext(os.path.basename(filename))[0],
        'time': time,
        'xGL': xgl,
        'd': d,
        'Hs_end': hs_end,
        'Hib_end': hib_end,
        'Hb_end': hb_end,
    }


def read_simulation_results(foldername):
    """Read all transect files for one simulation folder."""
    transect_files = list_transect_files(foldername)
    return [read_transect_results(fn) for fn in transect_files]


def setup_multipanel_figure(wa, ha, margins_hor, margins_ver):
    """
    Set up a multi-panel figure with specified layout.

    Parameters:
    -----------
    wa : int or list
        Width(s) of axes in pixels
    ha : int or list
        Height(s) of axes in pixels
    margins_hor : list
        Horizontal margins [left, ...inner..., right] in pixels
    margins_ver : list
        Vertical margins [top, ...inner..., bottom] in pixels

    Returns:
    --------
    fig : matplotlib.figure.Figure
        The figure object
    axes : dict
        Dictionary of axes indexed as (i, j)
    """

    # Handle scalar inputs for wa and ha
    if isinstance(wa, int):
        wa = [wa] * (len(margins_hor) - 1)
    if isinstance(ha, int):
        ha = [ha] * (len(margins_ver) - 1)

    nax = len(margins_hor) - 1
    nay = len(margins_ver) - 1

    # Figure size
    wf = sum(margins_hor) + sum(wa)
    hf = sum(margins_ver) + sum(ha)

    # Create figure (convert pixels to inches, assuming 100 dpi)
    dpi = 100
    fig = plt.figure(figsize=(wf / dpi, hf / dpi), dpi=dpi)
    fig.patch.set_facecolor('white')

    # Create axes
    axes = {}
    for i in range(nay):
        for j in range(nax):
            x = (sum(margins_hor[:j+1]) + sum(wa[:j])) / wf
            y_top = (sum(margins_ver[:i+1]) + sum(ha[:i])) / hf
            w = wa[j] / wf
            h = ha[i] / hf
            y = 1.0 - y_top - h  # Flip y-coordinate for matplotlib

            ax = fig.add_axes([x, y, w, h])
            ax.tick_params(labelsize=10)
            ax.grid(True, alpha=0.3)
            axes[(i, j)] = ax

    return fig, axes


def main():
    """Create MISMIP_mod geometry and grounding-line plots."""

    parser = argparse.ArgumentParser(
        description='Create MISMIP_mod benchmark figures from UFEMISM results.'
    )
    parser.add_argument(
        '--results-root',
        default='.',
        help='Path to the folder containing UFEMISM results_* directories.'
    )
    parser.add_argument(
        '--output-dir',
        default=None,
        help='Directory where output figure PNG files are written.'
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

    simulations = [
        ('spinup_10km', 'results_spinup_10km', 'r'),
        ('advance_10km', 'results_advance_10km', 'g'),
        ('retreat_10km', 'results_retreat_10km', 'b'),
    ]

    all_results = {}
    for sim_name, sim_folder, _ in simulations:
        all_results[sim_name] = read_simulation_results(os.path.join(results_root, sim_folder))

    # Set up figure
    wa = [500, 500]  # Width of two panels in pixels
    ha = 300  # Height in pixels
    margins_hor = [100, 35, 100]  # Left, middle, right margins
    margins_ver = [25, 80]  # Top, bottom margins

    fig, axes = setup_multipanel_figure(wa, ha, margins_hor, margins_ver)

    # ========== Left panel: Ice geometry ==========
    ax = axes[(0, 0)]
    ax.set_xlim([0, 1000])
    ax.set_ylim([-500, 4000])
    ax.set_xlabel('Distance from ice divide (km)', fontsize=11)
    ax.set_ylabel('z (m)', fontsize=11)

    # Left-panel legend handles
    for sim_name, _, color in simulations:
        ax.plot([], [], color=color, linestyle='-', linewidth=2, label=sim_name)
    ax.plot([], [], color='k', linestyle='-', linewidth=2)
    ax.plot([], [], color='k', linestyle='-', linewidth=2)
    ax.plot([], [], color='k', linestyle='-', linewidth=2)

    # Plot end-state geometry for all transects in each simulation
    for sim_name, _, color in simulations:
        for tr in all_results[sim_name]:
            d_km = tr['d'] / 1e3
            ax.plot(d_km, tr['Hs_end' ], color=color, linestyle='-', linewidth=1.5)
            ax.plot(d_km, tr['Hib_end'], color=color, linestyle='-', linewidth=1.5)
            ax.plot(d_km, tr['Hb_end' ], color='k'  , linestyle='-', linewidth=1.5)

    ax.legend(loc='upper right', fontsize=9, ncol=2)

    # ========== Right panel: Grounding-line position over time ==========
    ax = axes[(0, 1)]
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.set_xlim([20, 55])
    ax.set_ylim([550, 850])
    ax.set_xlabel('Time (kyr)', fontsize=11)
    ax.set_ylabel('Grounding-line position (km)', fontsize=11)

    # Grounding-line position over time for all transects and simulations
    for sim_name, _, _ in simulations:
        for tr in all_results[sim_name]:
            ax.plot(tr['time'] / 1e3, tr['xGL'] / 1e3, color='k', linewidth=1.25, alpha=0.75)

    # Save figure
    output_file = os.path.join(output_dir, 'Fig_integrated_test_MISMIP_mod_full.png')
    fig.savefig(output_file, dpi=100, bbox_inches='tight')
    print(f"Figure saved to {output_file}")
    plt.close(fig)


if __name__ == '__main__':
    main()
