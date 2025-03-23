import numpy as np
from numpy.typing import NDArray

from .utilities import is2d

###########################
# Position initialisation #
###########################


def uniform_random(
    amount_of_particles: int, box_dim: NDArray[np.float64], seed: int | None = None
):
    """
    Initialise `amount_of_particles` positions which fit in the specified box

    Parameters
    ----------
    amount_of_particles : int
        The amount of particles to generate
    box_dim : np.ndarray(float)
        The dimensions of the simulation box (x,y,z)
    seed : optional(int)
        The random seed to use, if `None` a seed is generated

    Returns
    -------
    np.ndarray of size (amount_of_particles,3) with positions inside the box
    """
    np.random.seed(seed)
    x_coordinates = np.random.uniform(0, box_dim[0], (amount_of_particles, 1))
    y_coordinates = np.random.uniform(0, box_dim[1], (amount_of_particles, 1))
    z_coordinates = np.random.uniform(0, box_dim[2], (amount_of_particles, 1))

    return np.hstack((x_coordinates, y_coordinates, z_coordinates))


def static(amount_of_particles: int, box_dim: NDArray[np.float64]):
    """
    Take the first `amount_of_particles` from a predefined array of locations
    rescaled to the size of the box.

    Parameters
    ----------
    amount_of_particles : int
        The amount of particles to return
    box_dim : np.ndarray(float)
        The dimensions of the simulation box (x,y,z)

    Returns
    -------
    np.ndarray of size (amount_of_particles,3) with positions inside the box
    """
    particles = np.array(
        [
            [0.1, 0.1, 0.1],
            [0.1, 0.1, 0.5],
            [0.1, 0.1, 0.9],
            [0.1, 0.5, 0.1],
            [0.1, 0.5, 0.5],
            [0.1, 0.5, 0.9],
            [0.1, 0.9, 0.1],
            [0.1, 0.9, 0.5],
            [0.1, 0.9, 0.9],
            [0.5, 0.1, 0.1],
            [0.5, 0.1, 0.5],
            [0.5, 0.1, 0.9],
            [0.5, 0.5, 0.1],
            [0.5, 0.5, 0.5],
            [0.5, 0.5, 0.9],
            [0.5, 0.9, 0.1],
            [0.5, 0.9, 0.5],
            [0.5, 0.9, 0.9],
            [0.9, 0.1, 0.1],
            [0.9, 0.1, 0.5],
            [0.9, 0.1, 0.9],
            [0.9, 0.5, 0.1],
            [0.9, 0.5, 0.5],
            [0.9, 0.5, 0.9],
            [0.9, 0.9, 0.1],
            [0.9, 0.9, 0.5],
            [0.9, 0.9, 0.9],
        ]
    )

    return particles[:amount_of_particles, :] * box_dim


def fcc_lattice(
    num_atoms: int,
    lat_const: float,
    corner_offset: np.ndarray[tuple[int], np.dtype[np.float64]],
) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
    """
    Initializes a system of atoms on an fcc lattice.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system
    lattice_const : float
        The lattice constant for an fcc lattice

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle coordinates
    """
    # we can initialise only bottom corners of the fcc cubes, this means there are 4 particles per each of these positions,
    # and the positions are spread around in 3d space.
    layers = int(np.ceil(np.pow(num_atoms / 4, 1 / 3)))

    # NOTE: I don't know if we want to do this or just start filling the grid and stop when we have enough
    if layers**3 * 4 != num_atoms:
        print(
            f"Illegal amount of atoms {num_atoms}, "
            + f"a full lattice can be constructed with {layers**3 * 4} or {(layers - 1) ** 3 * 4} instead"
        )
        exit(4)

    positions: list[list[float]] = []

    # I don't like these nested loops, but I don't think they matter that much since it'll still be O(n) while the sim is O(n^2)
    for x in range(layers * 2):
        for y in range(layers * 2):
            for z in range(layers * 2):
                if (x + y + z) % 2 == 0:
                    positions.append([x, y, z])

    nppos = np.array(positions, dtype=np.float64) * lat_const / 2 + corner_offset
    # check if the positions are a 2d array (for typing)
    if is2d(nppos):
        return nppos
    else:
        raise TypeError("initial positions array not 2D")


###########################
# Velocity initialisation #
###########################


def zero_speed(amount_of_particles):
    """
    Generate `amount_of_particles` zero-velocities

    Parameters
    ----------
    amount_of_particles: int
        How many particles need to be generated

    Returns
    -------
    np.ndarray
        Array of particle velocities
    """

    return np.zeros((amount_of_particles, 3))


def init_velocity(num_atoms, temp, seed=None):
    """
    Initializes the system with Gaussian distributed velocities.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system.
    temp : float
        The (unitless) temperature of the system.
    seed : int
        The random seed to use for numpy

    Returns
    -------
    vel_vec : np.ndarray
        Array of particle velocities
    """
    np.random.seed(seed)

    # Using Maxwell-Boltzmann distribution, scale = sqrt(T)
    scale = np.sqrt(temp)
    velocities = np.random.normal(loc=0.0, scale=scale, size=(num_atoms, 3))
    # Subtracting mean from the velocities
    velocities -= np.mean(velocities, axis=0)
    return velocities
