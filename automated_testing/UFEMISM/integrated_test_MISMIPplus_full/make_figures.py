#!/usr/bin/env python3
"""
Plot UFEMISM model output compared to MISMIP+ ensemble results.

This script reads model output from NetCDF files and creates a two-panel figure:
- Left panel: Ice geometry (surface and base) at initial and final states
- Right panel: Grounding-line position over time
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import os
import argparse


def read_mismip_ensemble(foldername):
    """
    Read and process the MISMIP+ model ensemble from Cornford et al. (2020)

    Cornford, S. L., Seroussi, H., Asay-Davis, X. S., Gudmundsson, G. H.,
    Arthern, R., Borstad, C., Christmann, J., Dias dos Santos, T.,
    Feldmann, J., Goldberg, D., Hoffman, M. J., Humbert, A., Kleiner, T.,
    Leguy, G., Lipscomb, W. H., Merino, N., Durand, G., Morlighem, M.,
    Pollard, D., Rückamp, M., Williams, C. R., and Yu, H.:
    Results of the third Marine Ice Sheet Model Intercomparison Project
    (MISMIP+), The Cryosphere 14, 2283--2301, 2020.
    """

    c2020 = []
    foldername = os.path.join(foldername, 'submission_data')
    experiment = 'Ice1r'

    list_of_models = list_models_for_experiment(foldername, experiment)

    for model_name in list_of_models:
        time, xgl_mid = read_xgl_mid_experiment_model(foldername, experiment, model_name)
        model = {
            'name': model_name,
            'time': time,
            'xGL': xgl_mid
        }
        c2020.append(model)

    return c2020


def list_models_for_experiment(foldername, experiment):
    """List all models for a given experiment in the folder."""

    list_of_models = []

    if not os.path.exists(foldername):
        return list_of_models

    for filename in os.listdir(foldername):
        if filename.startswith(experiment + '_') and filename.endswith('.nc'):
            model_name = filename[len(experiment) + 1:-3]
            list_of_models.append(model_name)

    return sorted(list_of_models)


def read_xgl_mid_experiment_model(foldername, experiment, model_name):
    """
    Read grounding-line position at model mid-stream for a specific model.
    """

    filename = os.path.join(foldername, f'{experiment}_{model_name}.nc')

    if not os.path.exists(filename):
        return np.array([]), np.array([])

    with netCDF4.Dataset(filename, 'r') as ds:
        # Read time and process to get relevant timeframes
        time = ds['time'][:]

        # Find end of relevant time frames (some models have lots of empty frames after)
        ti_end = 0
        for i in range(len(time)):
            if time[i] < 1e5:
                ti_end = i + 1
            else:
                break

        time = time[:ti_end]

        # Read grounding line coordinates
        # xGL and yGL have dimensions (nPointGL, nTime)
        xgl = ds['xGL'][:]
        ygl = ds['yGL'][:]

        # Extract midstream grounding line position
        y_mid = 4e4
        xgl_mid = np.zeros(ti_end)

        for ti in range(ti_end):
            # Get all points on grounding line at this time step
            xgl_ti = xgl[:, ti]
            ygl_ti = ygl[:, ti]

            # Remove NaNs and unrealistic values (large fill values)
            mask = ~np.isnan(xgl_ti) & ~np.isnan(ygl_ti) & (xgl_ti < 1e30)
            xgl_ti = xgl_ti[mask]
            ygl_ti = ygl_ti[mask]

            if len(xgl_ti) > 0:
                # Find point closest to y_mid
                j_mid = np.argmin(np.abs(ygl_ti - y_mid))
                xgl_mid[ti] = xgl_ti[j_mid]

        # Some models submitted data in km instead of m
        if np.max(xgl_mid) < 1e3:
            xgl_mid = xgl_mid * 1e3

    return time, xgl_mid


def read_ufemism_results(foldername):
    """Read UFEMISM model output from transect NetCDF file."""

    filename = os.path.join(foldername, 'transect_westeast.nc')

    ufe = {}

    if not os.path.exists(filename):
        raise FileNotFoundError(f"Required UFEMISM output not found: {filename}")

    with netCDF4.Dataset(filename, 'r') as ds:
        # Grounding-line position over time
        ufe['time'] = ds['time'][:]
        ufe['xGL'] = ds['grounding_line_distance_from_start'][:]

        # Coordinate and bedrock
        V = ds['V'][:]  # V has shape (2, n_spatial) where row 0 is x, row 1 is y
        ufe['x'] = V[0, :]  # x-coordinates (first row)
        ufe['Hb'] = ds['Hb'][0, :]  # Bedrock at first time step (or any time, should be same)

        # Initial geometry (first time step, all spatial points)
        ufe['Hs_init'] = ds['Hs'][0, :]  # Surface elevation at t=0
        ufe['Hib_init'] = ds['Hib'][0, :]  # Ice base at t=0

        # Final geometry (last time step, all spatial points)
        ti_end = len(ufe['time']) - 1
        ufe['Hs_end'] = ds['Hs'][ti_end, :]  # Surface elevation at final time
        ufe['Hib_end'] = ds['Hib'][ti_end, :]  # Ice base at final time

    return ufe


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
    """Main function to create the comparison figure."""

    parser = argparse.ArgumentParser(
        description='Create MISMIPplus benchmark figures from UFEMISM results.'
    )
    parser.add_argument(
        '--results-root',
        default='.',
        help='Path to the folder containing UFEMISM results_* directories.'
    )
    args = parser.parse_args()

    results_root = os.path.abspath(args.results_root)

    # Build paths from the script location to avoid CWD-dependent behavior.
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Read ensemble results
    foldername = os.path.abspath(
        os.path.join(script_dir, '../../../external/data/model_ensembles/MISMIPplus/')
    )
    c2020 = read_mismip_ensemble(foldername)

    # Read UFEMISM results
    ufe_5km = read_ufemism_results(os.path.join(results_root, 'results_5km_ice1r'))
    ufe_4km = read_ufemism_results(os.path.join(results_root, 'results_4km_ice1r'))

    # Set up figure
    wa = [500, 500]  # Width of two panels in pixels
    ha = 300  # Height in pixels
    margins_hor = [100, 25, 100]  # Left, middle, right margins
    margins_ver = [25, 80]  # Top, bottom margins

    fig, axes = setup_multipanel_figure(wa, ha, margins_hor, margins_ver)

    # ========== Left panel: Ice geometry ==========
    ax = axes[(0, 0)]
    ax.set_xlim([0, 800])
    ax.set_ylim([-1000, 1800])
    ax.set_xlabel('x (km)', fontsize=11)
    ax.set_ylabel('z (m)', fontsize=11)

    # Plot lines for legend
    ax.plot([], [], 'r-', linewidth=2, label='5 km')
    ax.plot([], [], 'b-', linewidth=2, label='4 km')
    ax.plot([], [], 'k-', linewidth=2, label='t = 0')
    ax.plot([], [], 'k--', linewidth=2, label='t = 100 yr')

    # 5 km results - initial geometry
    ax.plot(ufe_5km['x'] / 1e3, ufe_5km['Hs_init'], 'r-', linewidth=2)
    ax.plot(ufe_5km['x'] / 1e3, ufe_5km['Hib_init'], 'r-', linewidth=2)

    # 5 km results - final geometry
    ax.plot(ufe_5km['x'] / 1e3, ufe_5km['Hs_end'], 'r--', linewidth=2)
    ax.plot(ufe_5km['x'] / 1e3, ufe_5km['Hib_end'], 'r--', linewidth=2)

    # 4 km results - initial geometry
    ax.plot(ufe_4km['x'] / 1e3, ufe_4km['Hs_init'], 'b-', linewidth=2)
    ax.plot(ufe_4km['x'] / 1e3, ufe_4km['Hib_init'], 'b-', linewidth=2)

    # 4 km results - final geometry
    ax.plot(ufe_4km['x'] / 1e3, ufe_4km['Hs_end'], 'b--', linewidth=2)
    ax.plot(ufe_4km['x'] / 1e3, ufe_4km['Hib_end'], 'b--', linewidth=2)

    # Bedrock (plotted last so it's visible)
    ax.plot(ufe_4km['x'] / 1e3, ufe_4km['Hb'], 'k-', linewidth=3)

    ax.legend(loc='upper right', fontsize=10)

    # ========== Right panel: Grounding-line position over time ==========
    ax = axes[(0, 1)]
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.set_xlim([0, 100])
    ax.set_ylim([360, 470])
    ax.set_xlabel('Time (yr)', fontsize=11)
    ax.set_ylabel('Grounding-line position (km)', fontsize=11)

    # MISMIP+ ensemble
    for model in c2020:
        ax.plot(model['time'], model['xGL'] / 1e3, 'k-', linewidth=1, alpha=0.7)

    # UFEMISM results
    ax.plot(ufe_5km['time'], ufe_5km['xGL'] / 1e3, 'r-', linewidth=3, label='5 km')
    ax.plot(ufe_4km['time'], ufe_4km['xGL'] / 1e3, 'b-', linewidth=3, label='4 km')

    # Save figure
    output_file = 'Fig_integrated_test_MISMIPplus_full.png'
    fig.savefig(output_file, dpi=100, bbox_inches='tight')
    print(f"Figure saved to {output_file}")

    plt.show()


if __name__ == '__main__':
    main()
