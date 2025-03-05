import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from utilities import dprint


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
    print("Now plotting the energies")
    _ = plt.figure()
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
    print("Now plotting the distances")
    _ = plt.figure()
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


def create_animation(positions, timesteps, box_size, file_name="particles.mp4"):
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

    dprint(f"saving the animation to {file_name}")
    ani.save(file_name, writer="ffmpeg", fps=30)

