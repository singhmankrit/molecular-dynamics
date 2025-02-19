#!/usr/bin/env python
"""
This is a suggestion for structuring your simulation code properly.
However, it is not set in stone. You may modify it if you feel like
you have a good reason to do so.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D


debug = True if os.environ.get("DEBUG") is not None else False


def dprint(str):
    """
    Prints the passed string only if the debug environment variable is set.

    Parameters
    ----------
    str : string
        The string to print
    """
    if debug:
        print(str)


amount_of_particles = 3

init_pos = np.random.uniform(0.0, 1.0, (amount_of_particles, 3))
dprint(f"created {amount_of_particles} particles")
# atomic weight modified to be in kg
mass = 39.792 * 1.660539066e-27

epsilon = 119.8 * 1.380649e-23
sigma = 3.405e-10


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
    energies = []

    current_positions = init_pos
    current_velocities = init_vel

    for step in np.arange(num_tsteps):
        dprint(
            f"""
            at step {step} the particles are at {current_positions}
            the particles have velocities {current_velocities}
            """
        )

        # create the n-by-n matrix of all the distances and the n-by-n-by-3 matrix of the relative positions
        relative_positions, distances = atomic_distances(current_positions, box_dim)
        # get the n-by-3 matrix of all the total forces on the particles
        forces = lj_force(relative_positions, distances)
        # get current total energy and append
        energies.append(total_energy(distances, current_velocities))
        # Euler integration step, we'll have to rewrite this to improve energy conservation
        current_positions = (
            current_positions + current_velocities * timestep
        ) % box_dim
        current_velocities += forces * timestep / mass

        # append the new positions and velocities to the arrays
        positions.append(current_positions)
        velocities.append(current_velocities)

    return positions, velocities, energies


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

    relative_positions = np.stack([x_dist, y_dist, z_dist])
    distances = np.ma.masked_values(
        np.sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist),
        0.0,
        rtol=1e-60,
        atol=1e-60,
    )
    dprint(f"there are {np.ma.count_masked(distances)} masked distance values")

    # return the full distance matrix of shape n-by-n
    return (relative_positions, distances)


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
    # Compute force magnitude using the Lennard-Jones force formula
    force_magnitude = (
        24 * epsilon * (sigma / rel_dist) ** 7 - 48 * epsilon * (sigma / rel_dist) ** 13
    )
    dprint(
        f"the minimal force magnitude is {np.min(force_magnitude)} and the maximum force is {np.max(force_magnitude)}"
    )

    # Compute force matrix
    force_direction = rel_pos / rel_dist

    force_matrix = force_magnitude * force_direction

    # Sum forces acting on each particle
    net_force = -np.sum(force_matrix, axis=1)

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
    # Kinetic energy for each particle
    ke_individual = 1 / 2 * mass * np.sum(vel**2, axis=1)
    # Total Kinetic Energy
    ke = np.sum(ke_individual)
    return ke


def lj_potential(distance):
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


def total_energy(rel_dist, vel):
    """
    Computes the total energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances
    vel : np.ndarray
        Velocity of particle

    Returns
    -------
    float
        The total energy of the system.
    """
    return kinetic_energy(vel) + potential_energy(rel_dist)


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


def plot_energy(energies, file_name="energies.png"):
    """
    Plots the energy vs timesteps.

    Parameters
    ----------
    energies : list
        List of energies
    """
    plt.title("Energy vs timesteps")
    plt.xlabel("timesteps")
    plt.ylabel("Energy")
    plt.plot(energies)
    dprint(f"saving the energies plot to {file_name}")
    plt.savefig(file_name)


def create_animation(positions, timesteps, box_size, name="particles.mp4"):
    """
    Creates an animation of the system.

    Parameters
    ----------
    positions : list
        List of positions
    timesteps : int
        Number of timesteps
    box_size : np.ndarray
        Size of the simulation box
    name : str
        Name of the animation file
    """
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlim(0, box_size[0])
    ax.set_ylim(0, box_size[1])
    ax.set_zlim(0, box_size[2])
    positions = np.array(positions)  # Convert to NumPy array if it's a list
    # Scatter plot for particles
    (particles,) = ax.plot([], [], [], "bo", markersize=5)

    # Update function for animation
    def update(frame):
        particles.set_data(positions[frame, :, 0], positions[frame, :, 1])  # X, Y
        particles.set_3d_properties(positions[frame, :, 2])  # Z
        return (particles,)

    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=timesteps, blit=True)

    dprint(f"saving the animation to {name}")
    ani.save(name, writer="ffmpeg", fps=30)


if __name__ == "__main__":
    timesteps = 1000
    step_size = 0.01
    temp = 113.7
    box_size = np.array([1.0, 1.0, 1.0]) * 1e-6
    print(
        f"simulating the particles for {timesteps} timesteps, with a time step size of {step_size}"
    )
    pos, vel, energies = simulate(
        init_pos,
        np.zeros((amount_of_particles, 3)),
        timesteps,
        step_size,
        box_size,
    )
    print("finished simulating, plotting the energies over time")
    plot_energy(energies)
    print("plotted the energies, now creating the animation")
    create_animation(pos, timesteps, box_size)
