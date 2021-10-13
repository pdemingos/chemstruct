# Copyright 2019 Pedro G. Demingos

"""General tools."""

import numpy as np
from math import floor

from constants import axis_to_dim
from files.csv import Csv


def break_regions(points, cell, gap=6.0, struct="bonds"):
    """
    Takes a list of points in 3D space and divides them in sub-regions.
    System must be periodic to work properly. Cell can be triclinic.

    Parameters
    ----------
    points : list of lists or numpy.arrays
        Points in 3D space.
    cell : list of lists of floats
        Cell (i.e. lattice, 3x3 matrix) for division.
    gap : list of floats, optional
        Minimum linear size of sub-regions, in angstroms. Standard is 6.0.
    struct : str, optional
        Topological structures being computed. Only for printing purposes.
        Expected: "bonds", "angles" or "diehedrals". Standard is "bonds".

    Returns
    -------
    regions : list (4D matrix)
        Regions[x][y][z] is a list with the indices of all points in
        the (x, y, z) region.
    localizers : list (2D matrix)
        Localizers[n] is a list [x, y, z] with the indices of the
        nth point's region.
    moved_points : list of numpy.arrays
        New positions. They're changed if any point if outside the box,
        which doesn't change anything for a periodic system.

    """

    # gets all positions to be positive
    original_points = np.array(points)
    translation = -1 * np.array([min(original_points[:, 0]),
                                 min(original_points[:, 1]),
                                 min(original_points[:, 2])])
    moved_points = []
    for point in points:
        moved_points.append(point + translation)
    moved_points = np.array(moved_points)
    # moved_points = original_points
    cell = np.array(cell)

    # computes the size of sub-regions
    nx = int(cell[0][0] // gap)
    ny = int(cell[1][1] // gap)
    nz = int(cell[2][2] // gap)
    print("Regions for computing {}: ({}, {}, {})".format(
        struct, nx, ny, nz))
    da = cell[0] / nx  # vector
    db = cell[1] / ny  # vector
    dc = cell[2] / nz  # vector
    inv_matrix = np.linalg.inv(np.transpose(np.array([da, db, dc])))

    # instantiates region and localizers
    regions = []
    for x in range(nx):
        regions.append([])
        for y in range(ny):
            regions[x].append([])
            for z in range(nz):
                regions[x][y].append([])
    localizers = np.empty([len(moved_points), 3])

    # puts every point in a region
    for (index, point) in enumerate(moved_points):
        # va = proj(point, da)
        # vb = proj(point, db)
        # vc = proj(point, dc)
        # na = int(va[0] // da[0])
        # nb = int(vb[1] // db[1])
        # nc = int(vc[2] // dc[2])

        # changes the basis
        na, nb, nc = tuple(np.matmul(inv_matrix, point))
        na, nb, nc = floor(na), floor(nb), floor(nc)

        # moves points/atoms into the cell
        while na >= nx:
            moved_points[index] -= cell[0]
            na -= nx
        while na < 0:
            moved_points[index] += cell[0]
            na += nx
        while nb >= ny:
            moved_points[index] -= cell[1]
            nb -= ny
        while nb < 0:
            moved_points[index] += cell[1]
            nb += ny
        while nc >= nz:
            moved_points[index] -= cell[2]
            nc -= nz
        while nc < 0:
            moved_points[index] += cell[2]
            nc += nz
        # something might be needed for angles etc,
        # due to numerical precision of //
        # but so far we're fine

        regions[na][nb][nc].append(index)
        localizers[index][:] = [na, nb, nc]

    regions = np.array(regions)
    return regions, localizers, moved_points


def linear_counter(points: list, axis: str, bins: int, output_path=None):
    """
    Takes a list of points and returns histogram-like x-y information.

    Parameters
    ----------
    points : list
        Points in 3D space.
    axis : str
        Axis to count on, should be 'x', 'y' or 'z'.
    bins : int
        Number of divisions wanted for the axis.
    output_path : str, optional
        Path wanted for output csv file with the histogram-like data.

    Returns
    -------
    grid : list
        Center of each bin (i.e. x information in a histogram).
    stacks : list
        Counting for each bin (i.e. y information in a histogram).

    """
    dim = axis_to_dim(axis)
    if not isinstance(bins, int):
        raise TypeError("bins must be int, got {}".format(type(bins)))
    try:
        values = [p[dim] for p in points]
    except IndexError:
        raise IndexError("at least one point in points has dimension != 3")
    stacks, _ = np.histogram(values, bins)
    stacks = list(stacks)
    max_value, min_value = max(values), min(values)
    v_range = max_value - min_value
    bin_width = v_range / bins
    grid = [min_value + n * bin_width + bin_width / 2 for n in range(bins)]
    if isinstance(output_path, str):
        csv = Csv()
        csv.add_array(grid)
        csv.add_array(stacks)
        csv.write_csv(output_path)
    return grid, stacks

