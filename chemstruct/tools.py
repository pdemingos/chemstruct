"""General tools."""


import numpy as np
from math import floor
from glob import glob

from constants import axis_to_dim
from files.csv import Csv
from files.xyz import Xyz


def proj(b, a):
    # a and b must be numpy.arrays
    # not being used
    """Returns projection of vector b on vector a."""
    return np.dot(a, b) / np.linalg.norm(a)**2 * a


def break_regions(points, cell, gap=2.0, struct="bonds"):
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


def divide_points(points: list, box: list, region_size: list):
    # don't use this: use the break_regions function instead
    """
    Takes a list of points in 3D space and divides them in sub-regions.

    Parameters
    ----------
    points : list of lists or numpy.arrays
        Points in 3D space.
    box : list of floats
        Box sizes i.e. [lx, ly, lz].
    region_size : list of floats
        Minimum sizes of sub-regions.

    Returns
    -------
    regions : list (4D matrix)
        regions[x][y][z] is a list with the indices of all points in
        the (x, y, z) region.
    localizers : list (2D matrix)
        localizers[n] is a list [x, y, z] with the indices of the
        nth point's region.

    """
    original_points = np.array(points)
    # below we get all positions to be positive
    translation = -1 * np.array([min(original_points[:, 0]),
                                 min(original_points[:, 1]),
                                 min(original_points[:, 2])])
    moved_points = []
    for point in points:
        moved_points.append(point + translation)
    moved_points = np.array(moved_points)
    box = np.array(box)
    region_size = np.array(region_size)
    n_regions = box // region_size
    n_regions = n_regions.astype(int)
    region_size = box / n_regions
    regions = []  # I didn't want to do this with numpy.arrays
    for x in range(n_regions[0]):
        regions.append([])
        for y in range(n_regions[1]):
            regions[x].append([])
            for z in range(n_regions[2]):
                regions[x][y].append([])
    localizers = np.empty([len(points), 3])
    for (index, point) in enumerate(moved_points):
        x, y, z = tuple(point // region_size)
        x, y, z = int(x), int(y), int(z)
        try:
            regions[x][y][z].append(index)
            localizers[index][:] = [x, y, z]
        except IndexError:
            if x > n_regions[0] - 1:
                x -= 1
            if y > n_regions[1] - 1:
                y -= 1
            if z > n_regions[2] - 1:
                z -= 1
            # print("x = {}".format(x))
            # print("y = {}".format(y))
            # print("z = {}".format(z))
            # print("len_x = {}".format(len(regions)))
            # print("len_y = {}".format(len(regions[x])))
            # print("len_z = {}".format(len(regions[x][z])))
            regions[x][y][z].append(index)
            localizers[index][:] = [x, y, z]
    regions = np.array(regions)
    return regions, localizers


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


def post_filtration(path_to_cfg_dir: str, membrane_type: str,
                    types_filtrated):
    """
    Reads every LAMMPS CFG file in the directory, counting how many
    atoms of each type has already passed through the membrane;
    also gets the position of every atom inside the membrane.
    Writes a CSV file called "passed.csv" with the filtration countings,
    as well as an additional CSV file for each atom type inside
    the membrane (considering the whole simulation at once).

    Parameters
    ----------
    path_to_cfg_dir : str
        Path to the directory containing the CFG files.
    membrane_type : str
        Atom type (e.g. "C") of the membrane.
    types_filtrated : tuple
        Atom types to be counted to the right of the membrane, as well
        as inside it.

    Notes
    -----
    The filtration is expected to be in the x direction.
    The filtrated atoms are expected to the right of the membrane.
    The membrane is expected not to move during the simulation.

    Example:
    --------
    post_filtration("/mydir/cfg", "C", ("O", "Na", "Cl"))

    """

    passed = dict()
    inside_y = dict()
    inside_z = dict()
    for atom_type in types_filtrated:
        print(atom_type)
        passed[atom_type] = []
        inside_y[atom_type] = []
        inside_z[atom_type] = []

    cfgs = glob(path_to_cfg_dir + "/*.cfg")

    from files.lmp import Cfg
    cfg0 = Cfg(cfgs[0])  # membrane cannot move!
    membrane_xs = []  # filtration along x axis!
    for atom in cfg0.atoms.atoms:
        if atom.type == membrane_type:
            membrane_xs.append(atom.position[0])
    membrane_min = min(membrane_xs)
    membrane_max = max(membrane_xs)

    print(membrane_min, membrane_max)

    steps = []
    length = len(cfgs)
    for (i, cfg_path) in enumerate(cfgs):

        for atom_type in passed.keys():
            passed[atom_type].append(0)

        file_name = cfg_path.split("/")[-1]
        charmm, step, cfg = tuple(file_name.split("."))
        assert charmm == "charmm"
        assert cfg == "cfg"
        step = int(step)
        steps.append(step)

        try:
            cfg = Cfg(cfg_path)
        except IndexError or TypeError or NameError:
            print("WARNING: bad cfg file found")
            continue

        for atom in cfg.atoms.atoms:
            if atom.position[0] > membrane_max:  # to the right!
                try:
                    passed[str(atom.type)][-1] += 1
                except KeyError:
                    continue
            elif atom.position[0] > membrane_min:  # inside the membrane
                try:
                    inside_y[str(atom.type)].append(atom.position[1])
                    inside_z[str(atom.type)].append(atom.position[2])
                except KeyError:
                    continue
        print("Done: {}/{}".format(i, length))

    csv = Csv()
    csv.add_header("step")
    csv.add_array(steps)
    for (atom_type, count) in passed.items():
        csv.add_header(atom_type)
        csv.add_array(count)

    print(type(passed["O"]))
    print(passed["O"][0])

    csv.sort_by(0)  # i.e. step
    csv.write_csv(path_to_cfg_dir.replace("/cfg", "/passed.csv"))

    for atom_type in inside_y.keys():
        csv = Csv()
        csv.add_header(atom_type + "-y")
        csv.add_header(atom_type + "-z")
        csv.add_array(inside_y[atom_type])
        csv.add_array(inside_z[atom_type])
        csv.write_csv(path_to_cfg_dir.replace("/cfg",
                                              "/{}-in.csv".format(atom_type)))


def sort_graphene(xyz_path: str):
    xyz = Xyz(xyz_path)
    pass  # TODO
