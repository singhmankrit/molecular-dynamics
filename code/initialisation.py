import numpy as np

###########################
# Position initialisation #
###########################


def uniform_random(amount_of_particles, box_dim, seed=None):
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


def static(amount_of_particles, box_dim):
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


def fcc_lattice(num_atoms, lat_const):
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

    return


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

    # TODO: I don't know what the scale should be. (for week 4)
    scale = np.sqrt(temp)
    np.random.seed(seed)
    velocities = np.random.normal(loc=0.0, scale=scale, size=(num_atoms, 3))
    return velocities
