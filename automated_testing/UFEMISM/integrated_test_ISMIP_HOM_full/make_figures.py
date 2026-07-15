#!/usr/bin/env python3
"""Create ISMIP-HOM comparison figures for UFEMISM integrated test results."""

import argparse
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import netCDF4
import numpy as np


EXPERIMENTS = ['A', 'B', 'C', 'D']
LENGTH_SCALES = ['160', '80', '40', '20', '10', '5']
STOKES_APPROXS = ['SIASSA', 'DIVA', 'BPA']

LENGTH_SCALES_PLOT = [
    ['005', '010', '020'],
    ['040', '080', '160'],
]

COLOUR_BY_STOKES = {
    'SIASSA': (0.8, 0.1, 0.5),
    'DIVA': (0.8, 0.8, 0.1),
    'BPA': (0.2, 0.9, 0.9),
}

# Copied from external/data/model_ensembles/ISMIP-HOM/list_ISMIP_HOM_u_limits.m
U_LIMITS = {
    'A': {
        '160': (0.0, 120.0),
        '080': (0.0, 100.0),
        '040': (0.0, 70.0),
        '020': (0.0, 50.0),
        '010': (10.0, 35.0),
        '005': (11.0, 18.0),
    },
    'B': {
        '160': (0.0, 120.0),
        '080': (0.0, 120.0),
        '040': (0.0, 80.0),
        '020': (0.0, 60.0),
        '010': (5.0, 30.0),
        '005': (4.0, 14.0),
    },
    'C': {
        '160': (0.0, 200.0),
        '080': (0.0, 70.0),
        '040': (10.0, 35.0),
        '020': (13.0, 20.0),
        '010': (13.5, 17.0),
        '005': (6.0, 18.0),
    },
    'D': {
        '160': (0.0, 300.0),
        '080': (0.0, 140.0),
        '040': (10.0, 50.0),
        '020': (14.0, 24.0),
        '010': (15.0, 18.0),
        '005': (6.0, 18.0),
    },
}


def trim_leading_zeros(length_scale):
    """Trim leading zeros while preserving a valid numeric string."""

    trimmed = length_scale.lstrip('0')
    return trimmed if trimmed else '0'


def list_ismip_hom_u_limits(experiment, length_scale):
    """Return plotting limits matching the Pattyn et al. (2008) screenshots."""

    try:
        return U_LIMITS[experiment][length_scale]
    except KeyError as exc:
        raise ValueError(
            f'Invalid experiment/length-scale combination: {experiment}, {length_scale}'
        ) from exc


def read_ufe_results(results_dir, experiments, length_scales, stokes_approxs):
    """Read UFEMISM ISMIP-HOM transect data for all experiment combinations."""

    results = []

    for experiment in experiments:
        for length_scale in length_scales:
            for stokes_approx in stokes_approxs:
                filename = os.path.join(
                    results_dir,
                    f'transect_{experiment}_{length_scale}_{stokes_approx}.nc',
                )

                if not os.path.exists(filename):
                    raise FileNotFoundError(f'Required UFEMISM output not found: {filename}')

                with netCDF4.Dataset(filename, 'r') as ds:
                    x = np.asarray(ds['V'][0, :], dtype=float)
                    x = (x - np.min(x)) / (np.max(x) - np.min(x))
                    u_surf = np.asarray(ds['u_3D'][0, 0, :], dtype=float)

                results.append(
                    {
                        'experiment': experiment,
                        'length_scale': length_scale,
                        'Stokes_approx': stokes_approx,
                        'x': x,
                        'u_surf': u_surf,
                    }
                )

    return results


def resolve_results_dir(path):
    """Accept either a results folder or the parent test folder containing results/."""

    if os.path.isfile(os.path.join(path, 'transect_A_160_SIASSA.nc')):
        return path

    nested_results_dir = os.path.join(path, 'results')
    if os.path.isfile(os.path.join(nested_results_dir, 'transect_A_160_SIASSA.nc')):
        return nested_results_dir

    return path


def setup_multipanel_figure(wa, ha, margins_hor, margins_ver):
    """Set up figure/axes using the same pixel geometry as the MATLAB script."""

    nax = len(margins_hor) - 1
    nay = len(margins_ver) - 1

    if len(wa) != nax:
        raise ValueError('len(wa) must equal len(margins_hor) - 1')
    if len(ha) != nay:
        raise ValueError('len(ha) must equal len(margins_ver) - 1')

    wf = sum(margins_hor) + sum(wa)
    hf = sum(margins_ver) + sum(ha)

    dpi = 100
    fig = plt.figure(figsize=(wf / dpi, hf / dpi), dpi=dpi)
    fig.patch.set_facecolor('white')

    axes = {}
    for i in range(nay):
        for j in range(nax):
            x_px = sum(margins_hor[: j + 1]) + sum(wa[:j])
            y_top_px = sum(margins_ver[: i + 1]) + sum(ha[:i])

            x = x_px / wf
            w = wa[j] / wf
            h = ha[i] / hf
            y = 1.0 - (y_top_px / hf) - h

            ax = fig.add_axes([x, y, w, h])
            ax.tick_params(labelsize=14)
            ax.grid(True, alpha=0.35)
            axes[(i, j)] = ax

    return fig, axes


def main():
    """Create ISMIP-HOM benchmark figures from UFEMISM results."""

    parser = argparse.ArgumentParser(
        description='Create ISMIP-HOM benchmark figures from UFEMISM results.'
    )
    parser.add_argument(
        '--results-root',
        default='.',
        help='Path to the folder containing UFEMISM transect_*.nc output files.',
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

    ismip_hom_dir = os.path.join(repo_root, 'external', 'data', 'model_ensembles', 'ISMIP-HOM')
    screenshot_dir = os.path.join(ismip_hom_dir, 'screenshots')

    ufe_results = read_ufe_results(results_root, EXPERIMENTS, LENGTH_SCALES, STOKES_APPROXS)

    for experiment in EXPERIMENTS:
        figure_filename = os.path.join(output_dir, f'Fig_integrated_test_ISMIP_HOM_full_{experiment}.png')

        wa = [350, 350, 350]
        ha = [250, 250]
        margins_hor = [80, 60, 60, 25]
        margins_ver = [50, 40, 80]

        fig, axes = setup_multipanel_figure(wa, ha, margins_hor, margins_ver)

        # Empty lines to build legend labels exactly once.
        ax_legend = axes[(0, 0)]
        ax_legend.plot([], [], color=COLOUR_BY_STOKES['SIASSA'], linewidth=3, label='SIA/SSA')
        ax_legend.plot([], [], color=COLOUR_BY_STOKES['DIVA'], linewidth=3, label='DIVA')
        ax_legend.plot([], [], color=COLOUR_BY_STOKES['BPA'], linewidth=3, label='BPA')

        for i in range(len(LENGTH_SCALES_PLOT)):
            for j in range(len(LENGTH_SCALES_PLOT[0])):
                ax = axes[(i, j)]

                length_scale = LENGTH_SCALES_PLOT[i][j]
                length_scale_trim = trim_leading_zeros(length_scale)
                u_lim = list_ismip_hom_u_limits(experiment, length_scale)

                ax.set_xlim(0.0, 1.0)
                ax.set_xticks(np.arange(0.0, 1.01, 0.2))
                ax.set_ylim(u_lim[0], u_lim[1])
                ax.set_title(f'L = {length_scale_trim} km', fontsize=14)

                if i == 0:
                    ax.set_xticklabels([])
                if j == 0:
                    ax.set_ylabel('u (m yr$^{-1}$)', fontsize=14)
                if i == 1:
                    ax.set_xlabel('x / L', fontsize=14)

                screenshot_file = os.path.join(
                    screenshot_dir,
                    f'im_ISMIP_HOM_{experiment}_{length_scale}.png',
                )
                if not os.path.exists(screenshot_file):
                    raise FileNotFoundError(
                        f'Required ISMIP-HOM screenshot not found: {screenshot_file}'
                    )

                image = mpimg.imread(screenshot_file)
                ax.imshow(
                    image,
                    extent=(0.0, 1.0, u_lim[0], u_lim[1]),
                    aspect='auto',
                    interpolation='nearest',
                    zorder=0,
                )

                for run in ufe_results:
                    if (
                        run['experiment'] == experiment
                        and run['length_scale'] == length_scale_trim
                    ):
                        stokes = run['Stokes_approx']
                        if stokes not in COLOUR_BY_STOKES:
                            raise ValueError(f'Invalid Stokes approximation: {stokes}')

                        ax.plot(
                            run['x'],
                            run['u_surf'],
                            color=COLOUR_BY_STOKES[stokes],
                            linewidth=3,
                            zorder=2,
                        )

        ax_legend.legend(loc='best', fontsize=12)

        fig.savefig(figure_filename, dpi=100, bbox_inches='tight')
        print(f'Figure saved to {figure_filename}')
        plt.close(fig)


if __name__ == '__main__':
    main()
