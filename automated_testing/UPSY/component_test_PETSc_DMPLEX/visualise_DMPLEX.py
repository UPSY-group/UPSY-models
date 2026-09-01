#!/usr/bin/env python3
"""Visualise a two-dimensional PETSc DMPlex saved in an HDF5 file."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np


DEFAULT_FILENAME = Path(__file__).parent / "results" / "PETSc_DMPLEX_output.h5"


class DMPlexFormatError(ValueError):
    """Raised when a file does not use the expected PETSc DMPlex layout."""


def calculate_point_depths(point_cones: dict[int, list[int]]) -> dict[int, int]:
    """Return the DAG depth of each DMPlex point."""

    point_depths: dict[int, int] = {}

    def point_depth(point: int, ancestors: set[int]) -> int:
        if point in point_depths:
            return point_depths[point]
        if point in ancestors:
            raise DMPlexFormatError("the DMPlex topology contains a cycle")
        if point not in point_cones:
            raise DMPlexFormatError(f"point {point} is missing from /topology/order")

        cone = point_cones[point]
        if not cone:
            point_depths[point] = 0
        else:
            point_depths[point] = 1 + max(point_depth(child, ancestors | {point}) for child in cone)
        return point_depths[point]

    for point in point_cones:
        point_depth(point, set())

    return point_depths


def order_boundary_edges(edges: list[list[int]]) -> list[int]:
    """Order a closed sequence of two-vertex edges into a cell boundary."""

    if not edges or any(len(edge) != 2 for edge in edges):
        raise DMPlexFormatError("a two-dimensional cell must be bounded by two-vertex edges")

    polygon = [edges[0][0], edges[0][1]]
    unused_edges = edges[1:]

    while unused_edges:
        current_vertex = polygon[-1]
        matching_edges = [edge for edge in unused_edges if current_vertex in edge]
        if len(matching_edges) != 1:
            raise DMPlexFormatError("could not order the boundary edges of a cell")

        edge = matching_edges[0]
        unused_edges.remove(edge)
        next_vertex = edge[1] if edge[0] == current_vertex else edge[0]

        if not unused_edges:
            if next_vertex != polygon[0]:
                raise DMPlexFormatError("the boundary edges of a cell do not form a closed polygon")
        elif next_vertex == polygon[0]:
            raise DMPlexFormatError("the boundary edges close before all cell edges are used")
        else:
            polygon.append(next_vertex)

    return polygon


def extract_cell_vertices(
    cell: int,
    point_cones: dict[int, list[int]],
    point_depths: dict[int, int],
) -> list[int]:
    """Return the ordered vertex point IDs on a two-dimensional cell boundary."""

    cone = point_cones[cell]
    cone_depths = [point_depths[point] for point in cone]

    if all(depth == 0 for depth in cone_depths):
        return cone
    if all(depth == 1 for depth in cone_depths):
        return order_boundary_edges([point_cones[edge] for edge in cone])

    raise DMPlexFormatError(f"cell point {cell} has an unsupported mixed-depth cone")


def load_dmplex_mesh(filename: Path) -> tuple[list[np.ndarray], dict[int, np.ndarray]]:
    """Read vertices and two-dimensional cell boundaries from a PETSc DMPlex file."""

    try:
        with h5py.File(filename, "r") as h5_file:
            vertices = np.asarray(h5_file["/geometry/vertices"], dtype=float)
            cells_dataset = h5_file["/topology/cells"]
            cell_values = np.asarray(cells_dataset, dtype=int).ravel()
            cone_sizes = np.asarray(h5_file["/topology/cones"], dtype=int).ravel()
            point_ids = np.asarray(h5_file["/topology/order"], dtype=int).ravel()
            cell_dimension = int(cells_dataset.attrs["cell_dim"])
    except KeyError as exception:
        raise DMPlexFormatError(f"missing PETSc DMPlex dataset: {exception}") from exception

    if cell_dimension != 2:
        raise DMPlexFormatError(f"only two-dimensional DMPlex files are supported, got dimension {cell_dimension}")
    if vertices.ndim != 2 or vertices.shape[1] < 2:
        raise DMPlexFormatError("/geometry/vertices must have at least two coordinate columns")
    if len(point_ids) != len(cone_sizes):
        raise DMPlexFormatError("/topology/order and /topology/cones have incompatible lengths")
    if np.any(cone_sizes < 0) or int(cone_sizes.sum()) != len(cell_values):
        raise DMPlexFormatError("/topology/cones does not describe /topology/cells")
    if len(np.unique(point_ids)) != len(point_ids):
        raise DMPlexFormatError("/topology/order contains duplicate point IDs")

    point_cones: dict[int, list[int]] = {}
    offset = 0
    for point, cone_size in zip(point_ids, cone_sizes, strict=True):
        point_cones[int(point)] = [int(child) for child in cell_values[offset : offset + cone_size]]
        offset += cone_size

    point_depths = calculate_point_depths(point_cones)
    vertex_ids = [point for point in point_ids if point_depths[int(point)] == 0]
    if len(vertex_ids) != len(vertices):
        raise DMPlexFormatError(
            "/geometry/vertices does not contain one coordinate row for each zero-depth DMPlex point"
        )

    vertex_coordinates = {
        int(point): coordinates[:2] for point, coordinates in zip(vertex_ids, vertices, strict=True)
    }
    cell_ids = [point for point in point_ids if point_depths[int(point)] == cell_dimension]
    polygons = [
        np.asarray([vertex_coordinates[point] for point in extract_cell_vertices(int(cell), point_cones, point_depths)])
        for cell in cell_ids
    ]

    return polygons, vertex_coordinates


def deduplicate_polygons(polygons: list[np.ndarray]) -> list[np.ndarray]:
    """Remove geometrically identical cells written by replicated MPI ranks."""

    unique_polygons: list[np.ndarray] = []
    polygon_keys: set[tuple[tuple[float, float], ...]] = set()
    for polygon in polygons:
        key = tuple(sorted(tuple(float(value) for value in vertex) for vertex in polygon))
        if key not in polygon_keys:
            polygon_keys.add(key)
            unique_polygons.append(polygon)

    return unique_polygons


def plot_dmplex_mesh(
    polygons: list[np.ndarray], vertex_coordinates: dict[int, np.ndarray], point_labels: bool
) -> plt.Figure:
    """Create a matplotlib figure for a DMPlex mesh."""

    figure, axis = plt.subplots(figsize=(7, 7), constrained_layout=True)
    unique_polygons = deduplicate_polygons(polygons)
    axis.add_collection(
        PolyCollection(
            unique_polygons,
            facecolors="#d8ecef",
            edgecolors="#145a6b",
            linewidths=1.5,
        )
    )

    coordinates = np.asarray(list(vertex_coordinates.values()))
    unique_coordinates = np.unique(coordinates, axis=0)
    axis.scatter(
        unique_coordinates[:, 0],
        unique_coordinates[:, 1],
        color="#d1495b",
        edgecolor="#6d1f2a",
        linewidth=0.5,
        s=35,
        zorder=2,
    )

    if point_labels:
        for point, coordinates in vertex_coordinates.items():
            axis.annotate(str(point), coordinates, xytext=(5, 5), textcoords="offset points")

    axis.autoscale_view()
    axis.margins(0.1)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_title(f"PETSc DMPlex: {len(unique_polygons)} cells, {len(unique_coordinates)} vertices")

    return figure


def main() -> None:
    """Parse command-line arguments and display or save the mesh figure."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "filename",
        nargs="?",
        type=Path,
        default=DEFAULT_FILENAME,
        help=f"PETSc DMPlex HDF5 file (default: {DEFAULT_FILENAME})",
    )
    parser.add_argument("-o", "--output", type=Path, help="write the figure to this file instead of displaying it")
    parser.add_argument("--point-labels", action="store_true", help="label vertices with their DMPlex point IDs")
    arguments = parser.parse_args()

    try:
        polygons, vertex_coordinates = load_dmplex_mesh(arguments.filename)
    except (DMPlexFormatError, OSError) as exception:
        parser.error(str(exception))

    figure = plot_dmplex_mesh(polygons, vertex_coordinates, arguments.point_labels)
    if arguments.output is None:
        plt.show()
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(arguments.output, dpi=180)
        print(f"Wrote {arguments.output}")


if __name__ == "__main__":
    main()