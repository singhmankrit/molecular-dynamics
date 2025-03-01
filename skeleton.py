#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import initialisation
import 
from utilities import parse_config, dprint



def plot_energy(kinetic, potential, file_name="energies.png"):
    """
    Plots the energy vs timesteps.

    Parameters
    ----------
    kinetic : list
        List of kinetic energies
    potential : list
        List of potential energies
    file_name : string
        Name of file to save to
    """
    plt.title("Energy vs timesteps")
    plt.xlabel("timesteps")
    plt.ylabel("Energy")
    plt.plot(kinetic, label="kinetic", color="orange")
    plt.plot(potential, label="potential", color="purple")
    plt.plot(np.array(kinetic) + np.array(potential), label="total", color="black")
    plt.legend()
    plt.tight_layout()
    dprint(f"saving the energies plot to {file_name}")
    plt.savefig(file_name)
    plt.close()


def plot_distances(distance_list, particle=0, file_name="distances.png"):
    """
    Plots the distances between particles.

    Parameters
    ----------
    distance_list : list
        List of distances
    particle : int
        Particle to plot distances from
    file_name : string
        Name of file to save to
    """
    plt.title(f"Distances between particle {particle} and other particles")
    plt.xlabel("timesteps")
    plt.ylabel("Distance")
    for i in range(0, len(distance_list[0])):
        if i == particle:
            continue
        plt.plot(np.array(distance_list)[:, particle, i], label=f"Particle {i}")
    plt.tight_layout()
    plt.legend()
    dprint(f"saving the distances plot to {file_name}")
    plt.savefig(file_name)
    plt.close()


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
    (
        amount_of_particles,
        step_size,
        timesteps,
        temperature,
        box_size,
        random_seed,
        position_init_method,
        velocity_init_method,
        simulator_type,
    ) = parse_config("config.json")
    print(
        f"simulating {amount_of_particles} particles for {timesteps} timesteps, with a time step size of {step_size}"
    )
    init_pos = None
    if position_init_method == "uniform":
        init_pos = initialisation.uniform_random(
            amount_of_particles, box_size, random_seed
        )
    elif position_init_method == "static":
        init_pos = initialisation.static(amount_of_particles, box_size)
    else:
        print(
            f"Please select a valid position init method ('uniform', 'static'), currently: {position_init_method}"
        )
        exit(2)
    init_vel = None
    if velocity_init_method == "zero":
        init_vel = initialisation.zero_speed(amount_of_particles)
    elif velocity_init_method == "mbdist":
        init_vel = initialisation.init_velocity(
            amount_of_particles, temperature, random_seed
        )
    else:
        print(
            f"Please select a valid velocity init method ('zero', 'mbdist'), currently: {velocity_init_method}"
        )
        exit(3)
        
    simulator = None
    if simulator_type == "verlet":
        simulator = simulators.verlet
    elif simulator_type == "euler":
        simulator = simulators.euler
    else:
        print(
            f"Please select a valid simulator ('verlet', 'euler'), currently: {simulator_type}"
        )
        exit(4)
    

    pos, vel, kinetic, potential, distance_list = simulator(
        init_pos,
        init_vel,
        timesteps,
        step_size,
        box_size,
    )
    print("finished simulating, plotting the energies over time")
    plot_energy(kinetic, potential)
    print("plotted the energies, now plotting the distances")
    plot_distances(distance_list)
    print("plotted the distances, now creating the animation")
    create_animation(pos, timesteps, box_size)
