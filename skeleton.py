#!/usr/bin/env python
"""
This is a suggestion for structuring your simulation code properly.
However, it is not set in stone. You may modify it if you feel like
you have a good reason to do so.
"""

import numpy as np

init_pos = np.array(
    [
        [0.2, 0.5, 0.0],
        [0.5, 0.8, 0.0],
    ]
)


def simulate(init_pos, init_vel, num_tsteps, timestep, box_dim):
    """
    Molecular dynamics simulation using the Euler or Verlet's algorithms
    to integrate the equations of motion. Calculates energies and other
    observables at each timestep.

    Parameters
    ----------
    init_pos : np.ndarray
        The initial positions of the atoms in Cartesian space
    init_vel : np.ndarray
        The initial velocities of the atoms in Cartesian space
    num_tsteps : int
        The total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : np.ndarray(float)
        Dimensions of the simulation box

    Returns
    -------
    Any quantities or observables that you wish to study.
    """
    positions = [init_pos]
    velocities = [init_vel]

    # atomic weight modified to be in kg
    mass = 39.792 * 1.660539066e-27
    current_positions = init_pos
    current_velocities = init_vel

    for step in np.arange(num_tsteps):
        # create the n-by-n matrix of all the distances and the n-by-n-by-3 matrix of the relative positions
        relative_positions, distances = atomic_distances(current_positions, box_dim)
        # get the n-by-3 matrix of all the total forces on the particles
        forces = lj_force(relative_positions, distances)
        # Euler integration step, we'll have to rewrite this to improve energy conservation
        current_positions = (
            current_positions + current_velocities * timestep
        ) % box_dim
        current_velocities += forces * timestep / mass

        # append the new positions and velocities to the arrays after PR 1 is merged
        # positions.append(current_positions)
        # velocities.append(current_velocities)

    return


def atomic_distances(
    pos: np.typing.NDArray[np.float64], box_dim: np.typing.NDArray[np.float64]
) -> np.typing.NDArray[np.float64]:
    """
    Calculates relative positions and distances between particles.

    parameters
    ----------
    pos : np.ndarray
        The positions of the particles in cartesian space
    box_dim : float
        The dimension of the simulation box

    returns
    -------
    rel_pos : np.ndarray
        Relative positions of particles
    rel_dist : np.ndarray
        The distance between particles
    """
    # create meshgrids to make an n-by-n matrix of distances
    central_x, other_x = np.meshgrid(pos[:, 0], pos[:, 0])
    central_y, other_y = np.meshgrid(pos[:, 1], pos[:, 1])
    central_z, other_z = np.meshgrid(pos[:, 2], pos[:, 2])
    # moving to the coordinate frame of the central particle
    # to find the closest position of those around
    x_dist = (central_x - other_x + box_dim[0] / 2) % box_dim[0] - box_dim[0] / 2
    y_dist = (central_y - other_y + box_dim[1] / 2) % box_dim[1] - box_dim[1] / 2
    z_dist = (central_z - other_z + box_dim[2] / 2) % box_dim[2] - box_dim[2] / 2

    # return the full distance matrix of shape n-by-n
    return (
        np.stack([x_dist, y_dist, z_dist]),
        np.ma.masked_values(
            np.sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist), 0.0
        ),
    )


def lj_force(rel_pos, rel_dist):
    """
    Calculates the net forces on each atom from the matrices containing the positions and distances.

    Parameters
    ----------
    rel_pos : np.ndarray
        Relative particle positions as obtained from atomic_distances
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    np.ndarray
        nx3 array having the net vector force acting on particle i due to all other particles
    """
    epsilon = 119.8 * 1.380649e-23
    sigma = 3.405e-10

    np.fill_diagonal(rel_dist, 0.0)
    
    # Compute force magnitude using the Lennard-Jones force formula
    force_magnitude = np.where(rel_dist!=0,(24 * epsilon / rel_dist**7) * (sigma**6) * (1 - 2 * (sigma**6 / rel_dist**6)),0)

    # Compute force matrix
    force_direction = np.where(rel_dist!=0,rel_pos / rel_dist,0)

    force_matrix = force_magnitude * force_direction
    
    # Sum forces acting on each particle
    net_force = np.sum(force_matrix, axis=1)  
    
    return net_force.T 
    


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


def kinetic_energy(vel):
    """
    Computes the kinetic energy of an atomic system.

    Parameters
    ----------
    vel: np.ndarray
        Velocity of particle

    Returns
    -------
    float
        The total kinetic energy of the system.
    """

    return


def lj_potential(distance):
    epsilon = 119.8 * 1.380649e-23
    sigma = 3.405e-10
    return 4 * epsilon * ((sigma / distance) ** 12 - (sigma / distance) ** 6)


def potential_energy(rel_dist):
    """
    Computes the potential energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    float
        The total potential energy of the system.
    """

    return 1 / 2 * np.sum(lj_potential(rel_dist))


def init_velocity(num_atoms, temp):
    """
    Initializes the system with Gaussian distributed velocities.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system.
    temp : float
        The (unitless) temperature of the system.

    Returns
    -------
    vel_vec : np.ndarray
        Array of particle velocities
    """

    # NOTE: I don't know what the scale should be, it should also be a maxwell-boltzmann distribution to be
    # fully correct. Though that'd require scipy or someone finding how
    scale = temp
    velocities = np.random.normal(loc=0.0, scale=scale, size=(num_atoms, 3))
    return velocities


if __name__ == "main":
    timesteps = 1000
    step_size = 0.01
    temp = 113.7
    simulate(
        init_pos,
        init_velocity(len(init_pos), temp),
        timesteps,
        step_size,
        np.array([1.0, 1.0, 1.0]),
    )
