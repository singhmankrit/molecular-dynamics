import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from utilities import dprint


def plot_energy(kinetic, potential, step_size, timesteps, eq_timestep, file_name="energies.png"):
    """
    Plots the energy vs timesteps.

    Parameters
    ----------
    kinetic : list
        List of kinetic energies
    potential : list
        List of potential energies
    step_size : float
        Step size of simulation
    timesteps : int
        Number of timesteps
    eq_timestep : int
        Equilibrium timestep
    file_name : string
        Name of file to save to
    """
    print("Now plotting the energies")
    time = np.arange(timesteps+1) * step_size
    eq_time = eq_timestep * step_size
    _ = plt.figure(figsize=(10, 7))
    plt.title("Energy vs Time")
    plt.xlabel(r"Time  $\left ( \sqrt{\frac{m\sigma^2}{\epsilon}}\right )$ ")
    plt.ylabel(r"Energy ($\epsilon$)")
    plt.plot(time, kinetic, label="kinetic", color="orange")
    plt.plot(time, potential, label="potential", color="purple")
    plt.plot(time, np.array(kinetic) + np.array(potential),
             label="total", color="black")
    plt.axvline(eq_time, color="green", linestyle="--", label="equilibrium")
    plt.legend()
    plt.tight_layout()
    dprint(f"saving the energies plot to {file_name}")
    plt.savefig(file_name)
    plt.close()


def plot_distances(distance_list, step_size, timesteps, eq_timestep, particle=0, file_name="distances.png"):
    """
    Plots the distances between particles.

    Parameters
    ----------
    distance_list : list
        List of distances
    step_size : float
        Step size of simulation
    timesteps : int
        Number of timesteps
    eq_timestep : int
        Equilibrium timestep
    particle : int
        Particle to plot distances from
    file_name : string
        Name of file to save to
    """
    print("Now plotting the distances")
    time = np.arange(timesteps+1) * step_size
    eq_time = eq_timestep * step_size
    _ = plt.figure(figsize=(10, 7))
    plt.title(f"Distances between particle {particle} and other particles")
    plt.xlabel(r"Time  $\left ( \sqrt{\frac{m\sigma^2}{\epsilon}}\right )$ ")
    plt.ylabel(r"Distance ($\sigma$)")
    for i in range(0, len(distance_list[0])):
        if i == particle:
            continue
        plt.plot(time, np.array(distance_list)[
                 :, particle, i], label=f"Particle {i}")
    plt.axvline(eq_time, color="green", linestyle="--", label="equilibrium")
    plt.tight_layout()
    plt.legend()
    dprint(f"saving the distances plot to {file_name}")
    plt.savefig(file_name)
    plt.close()


def create_animation(
    positions, timesteps, step_size, eq_timestep, box_size, selected=0, file_name="particles.mp4"
):
    """
    Creates an animation of the system.

    Parameters
    ----------
    positions : list
        List of positions
    timesteps : int
        Number of timesteps
    step_size : float
        Step size of simulation
    eq_timestep : int
        Equilibrium timestep
    box_size : np.ndarray
        Size of the simulation box
    name : str
        Name of the animation file
    """
    print("Now creating the animation")
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlim(0, box_size[0])
    ax.set_ylim(0, box_size[1])
    ax.set_zlim(0, box_size[2])
    positions = np.array(positions)  # Convert to NumPy array if it's a list
    n_particles = positions.shape[1]
    # Scatter plot for particles
    particles = ax.scatter(
        positions[0, :, 0],
        positions[0, :, 1],
        positions[0, :, 2],
        c=["r" if i == selected else "b" for i in range(n_particles)],
    )
    title = ax.set_title("Time = 0.00")
    # Update function for animation

    def update(frame):
        time = frame * step_size
        particles._offsets3d = (
            positions[frame, :, 1],
            positions[frame, :, 0],
            positions[frame, :, 2],
        )
        if frame < eq_timestep:
            title.set_text(f"Time = {time:.2f}")
        else:
            title.set_text(
                f"Time = {time:.2f} \n Reached Equilibrium at Time = {eq_timestep * step_size:.2f}")
        return particles, title

    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=timesteps, blit=True)

    dprint(f"saving the animation to {file_name}")
    ani.save(file_name, writer="ffmpeg", fps=30)


def plot_pair_correlation(r_values, g_r, file_name="pair_correlation.png"):
    """
    Plots the pair correlation function.

    Parameters
    ----------
    r_values : list
        List of center of bins
    g_r : list
        List of pair correlation function values
    file_name : string
        Name of file to save to
    """
    print("Now plotting the pair correlation function")

    # Plot the results
    _ = plt.figure(figsize=(10, 7))
    plt.title(r"Pair Correlation Function ($g(r)$) vs Distance")
    plt.xlabel(r"Distance ($\sigma$)")
    plt.ylabel(r"Pair Correlation Function ($g(r)$)")
    plt.plot(r_values, g_r, label="Pair Correlation Function")
    plt.legend()
    plt.tight_layout()
    dprint(f"saving the pair correlation plot to {file_name}")
    plt.savefig(file_name)
    plt.close()
