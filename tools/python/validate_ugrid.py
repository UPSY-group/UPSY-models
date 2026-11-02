#!/usr/bin/env python
"""Validate and plot 2D UGRID NetCDF files (e.g. those written by UPSY's save_ugrid_as_netcdf).

Files are opened with xugrid's standard reader, so a file that does not follow the
UGRID conventions fails to load; fix the writer in that case, not this script.
After loading, a number of planar checks that xugrid itself does not do are run:

  - every face has >= 3 distinct nodes
  - no duplicate (coinciding) nodes, no unused nodes
  - all faces are oriented anticlockwise (positive signed area)
  - all faces are convex (required by the convex polygon clipper used for remapping)
  - every edge is shared by at most two faces (manifold), and the boundary edges
    form closed loops (watertight)
  - if the file has a face_face_connectivity variable, it matches the
    connectivity xugrid derives from the face-node connectivity

Usage:
  python validate_ugrid.py file1.nc [file2.nc ...] [--plot] [--save-dir DIR]

Exits with status 1 if any file fails to load or fails a check.
"""

import argparse
import os
import sys

import warnings

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import xugrid as xu
from xugrid.constants import FILL_VALUE


def signed_face_areas(grid):
    """Signed area of every face (shoelace formula); positive = anticlockwise."""
    xy = grid.face_node_coordinates  # (n_face, n_max_node, 2), NaN-padded
    x, y = xy[..., 0], xy[..., 1]
    n = grid.n_node_per_face
    idx = np.arange(xy.shape[1])
    # Index of the next node of each face, wrapping around at the face's own node count
    nxt = np.where(idx[None, :] + 1 < n[:, None], idx[None, :] + 1, 0)
    x2 = np.take_along_axis(x, nxt, axis=1)
    y2 = np.take_along_axis(y, nxt, axis=1)
    valid = idx[None, :] < n[:, None]
    return 0.5 * np.sum(np.where(valid, x * y2 - x2 * y, 0.0), axis=1)


def nonconvex_faces(grid, tol):
    """Indices of faces with a reflex (clockwise) corner, assuming anticlockwise orientation."""
    xy = grid.face_node_coordinates
    n = grid.n_node_per_face
    idx = np.arange(xy.shape[1])
    nxt = np.where(idx[None, :] + 1 < n[:, None], idx[None, :] + 1, 0)
    prv = np.where(idx[None, :] > 0, idx[None, :] - 1, n[:, None] - 1)
    p = xy
    q = np.take_along_axis(xy, nxt[..., None], axis=1)
    r = np.take_along_axis(xy, prv[..., None], axis=1)
    e1 = p - r
    e2 = q - p
    cross = e1[..., 0] * e2[..., 1] - e1[..., 1] * e2[..., 0]
    # Normalise to a distance so the tolerance has length units
    cross = cross / np.maximum(np.hypot(e1[..., 0], e1[..., 1]), np.finfo(float).tiny)
    valid = idx[None, :] < n[:, None]
    return np.flatnonzero(np.any(valid & (cross < -tol), axis=1))


def check_grid(ds, grid, tol_rel):
    """Run the planar checks on one Ugrid2d. Returns (errors, info, bad_faces)."""
    errors = []
    info = []
    bad_faces = set()

    # xugrid converts all connectivity to 0-based with its own fill value,
    # whatever start_index and _FillValue the file uses
    fnc = grid.face_node_connectivity
    fill = FILL_VALUE
    n_per_face = grid.n_node_per_face

    xmin, ymin, xmax, ymax = grid.bounds
    tol = tol_rel * max(xmax - xmin, ymax - ymin)

    info.append(f'{grid.n_node} nodes, {grid.n_face} faces, {grid.n_edge} edges')
    info.append(f'nodes per face: {n_per_face.min()} - {n_per_face.max()}')

    # Degenerate faces
    few = np.flatnonzero(n_per_face < 3)
    if few.size:
        errors.append(f'{few.size} faces have fewer than 3 nodes, e.g. {few[:5].tolist()}')
        bad_faces.update(few.tolist())
    rep = [i for i in range(grid.n_face)
           if np.unique(fnc[i, :n_per_face[i]]).size < n_per_face[i]]
    if rep:
        errors.append(f'{len(rep)} faces contain the same node more than once, e.g. {rep[:5]}')
        bad_faces.update(rep)

    # Unused and duplicate nodes
    used = np.zeros(grid.n_node, dtype=bool)
    used[fnc[fnc != fill]] = True
    if not used.all():
        unused = np.flatnonzero(~used)
        errors.append(f'{unused.size} nodes are not used by any face, e.g. {unused[:5].tolist()}')
    dist, _ = grid.node_kdtree.query(grid.node_coordinates, k=2)
    dup = np.flatnonzero(dist[:, 1] < tol)
    if dup.size:
        errors.append(f'{dup.size} nodes coincide with another node (within {tol:.3g}), '
                      f'e.g. {dup[:5].tolist()}')

    # Orientation and convexity
    area = signed_face_areas(grid)
    cw = np.flatnonzero(area <= 0)
    if cw.size:
        errors.append(f'{cw.size} faces are not anticlockwise (signed area <= 0), '
                      f'e.g. {cw[:5].tolist()}')
        bad_faces.update(cw.tolist())
    ncv = nonconvex_faces(grid, tol)
    if ncv.size:
        errors.append(f'{ncv.size} faces are not convex, e.g. {ncv[:5].tolist()}')
        bad_faces.update(ncv.tolist())
    info.append(f'total face area: {area.sum():.6g} '
                f'(bounding box area: {(xmax - xmin) * (ymax - ymin):.6g})')

    # Manifold and watertight
    efc = grid.edge_face_connectivity  # (n_edge, 2), second column = fill for boundary edges
    n_edge_faces = np.sum(efc != fill, axis=1)
    if np.any(n_edge_faces > 2):
        errors.append(f'{np.sum(n_edge_faces > 2)} edges are shared by more than two faces')
    boundary_edges = grid.edge_node_connectivity[n_edge_faces == 1]
    counts = np.bincount(boundary_edges.ravel(), minlength=grid.n_node)
    bad_nodes = np.flatnonzero((counts != 0) & (counts != 2))
    if bad_nodes.size:
        errors.append(f'boundary is not a set of closed loops: {bad_nodes.size} boundary nodes '
                      f'have {set(counts[bad_nodes].tolist())} boundary edges instead of 2, '
                      f'e.g. nodes {bad_nodes[:5].tolist()}')
    info.append(f'{len(boundary_edges)} boundary edges')

    # File's own face_face_connectivity (if any) versus the one xugrid derives
    # (xugrid does not keep the file's own connectivity variables, so read it from the
    # raw dataset; xarray's standard decoding turns _FillValue entries into NaN)
    ffc_name = ds[grid.name].attrs.get('face_face_connectivity')
    if ffc_name is not None and ffc_name in ds.variables:
        ffc_var = ds[ffc_name]
        start = int(ffc_var.attrs.get('start_index', 0))
        ffc_file = ffc_var.transpose(grid.face_dimension, ...).values
        derived = grid.face_face_connectivity.tocsr()
        mismatch = []
        for i in range(grid.n_face):
            row = ffc_file[i]
            row = row[~np.isnan(row)] if row.dtype.kind == 'f' else row
            mine = set((row.astype(int) - start).tolist())
            theirs = set(derived.indices[derived.indptr[i]:derived.indptr[i + 1]].tolist())
            if mine != theirs:
                mismatch.append(i)
        if mismatch:
            errors.append(f'{len(mismatch)} faces have a face_face_connectivity that differs from '
                          f'the one derived from face_node_connectivity, e.g. {mismatch[:5]}')
            bad_faces.update(mismatch)
    else:
        info.append('no face_face_connectivity in file; skipped that check')

    return errors, info, sorted(bad_faces)


def plot_grid(grid, bad_faces, title, filename=None):
    fig, ax = plt.subplots(figsize=(8, 8))
    grid.plot(ax=ax, color='k', linewidth=0.3)
    ax.plot(grid.face_x, grid.face_y, '.', color='tab:blue', markersize=1.5)
    if bad_faces:
        for i in bad_faces:
            n = grid.n_node_per_face[i]
            xy = grid.face_node_coordinates[i, :n]
            ax.fill(xy[:, 0], xy[:, 1], color='tab:red', alpha=0.6)
        ax.set_title(f'{title}\n{len(bad_faces)} faces failed checks (red)')
    else:
        ax.set_title(title)
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=200)
        plt.close(fig)


def validate_file(filename, args):
    print(f'\n{filename}')
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            uds = xu.open_dataset(filename)
    except Exception as e:
        print(f'  FAIL: xugrid could not open the file as UGRID: {type(e).__name__}: {e}')
        return False
    # Warnings from xugrid's reader point at missing or non-standard metadata
    messages = [' '.join(str(w.message).split()) for w in caught
                if 'xugrid' in w.filename]
    for message in dict.fromkeys(messages):
        print(f'  WARNING (xugrid): {message}')
    ds = xr.open_dataset(filename)

    grids = [g for g in uds.ugrid.grids if isinstance(g, xu.Ugrid2d)]
    if not grids:
        print('  FAIL: file contains no 2D UGRID mesh topology')
        return False

    ok = True
    for grid in grids:
        errors, info, bad_faces = check_grid(ds, grid, args.tol)
        print(f'  mesh "{grid.name}":')
        for line in info:
            print(f'    {line}')
        for line in errors:
            print(f'    FAIL: {line}')
        if not errors:
            print('    OK')
        ok = ok and not errors

        if args.plot or args.save_dir:
            title = f'{os.path.basename(filename)} : {grid.name}'
            out = None
            if args.save_dir:
                os.makedirs(args.save_dir, exist_ok=True)
                base = os.path.splitext(os.path.basename(filename))[0]
                out = os.path.join(args.save_dir, f'{base}_{grid.name}.png')
            plot_grid(grid, bad_faces, title, out)

    return ok


def main():
    parser = argparse.ArgumentParser(description='Validate and plot 2D UGRID NetCDF files')
    parser.add_argument('files', nargs='+', help='UGRID NetCDF file(s)')
    parser.add_argument('--plot', action='store_true', help='show a plot of each mesh')
    parser.add_argument('--save-dir', help='save a plot of each mesh as PNG in this folder')
    parser.add_argument('--tol', type=float, default=1e-9,
                        help='geometric tolerance, relative to the domain size (default: 1e-9)')
    args = parser.parse_args()

    results = [validate_file(f, args) for f in args.files]

    n_fail = results.count(False)
    print(f'\n{len(results) - n_fail}/{len(results)} files passed')

    if args.plot:
        plt.show()

    sys.exit(1 if n_fail else 0)


if __name__ == '__main__':
    main()
