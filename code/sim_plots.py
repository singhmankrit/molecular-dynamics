import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from alive_progress import alive_bar
from numpy.typing import NDArray
from scipy import optimize as opt
import csv
from tabulate import tabulate

from .utilities import dprint


def plot_energy(
    kinetic: NDArray[np.float64],
    potential: NDArray[np.float64],
    step_size: float,
    timesteps: int,
    eq_timestep: int,
    file_name: str = "energies.png",
):
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
    time = np.arange(timesteps + 1) * step_size
    eq_time = eq_timestep * step_size
    _ = plt.figure(figsize=(10, 7))
    plt.title("Energy vs Time")
    plt.xlabel(r"Time  $\left ( \sqrt{\frac{m\sigma^2}{\epsilon}}\right )$ ")
    plt.ylabel(r"Energy ($\epsilon$)")
    plt.plot(time, kinetic, label="kinetic", color="orange")
    plt.plot(time, potential, label="potential", color="purple")
    plt.plot(
        time, np.array(kinetic) + np.array(potential), label="total", color="black"
    )
    plt.axvline(eq_time, color="green", linestyle="--", label="equilibrium")
    plt.legend()
    plt.tight_layout()
    dprint(f"saving the energies plot to {file_name}")
    plt.savefig(file_name)
    plt.close()


def create_animation(
    positions: NDArray[np.float64],
    timesteps: int,
    step_size: float,
    eq_timestep: int,
    box_size: NDArray[np.float64],
    selected: int = 0,
    file_name: str = "particles.mp4",
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
    with alive_bar(timesteps + 2) as bar:

        def update(frame):
            time = frame * step_size
            particles._offsets3d = (
                positions[frame, :, 0] % box_size[0],
                positions[frame, :, 1] % box_size[1],
                positions[frame, :, 2] % box_size[2],
            )
            if frame < eq_timestep:
                title.set_text(f"Time = {time:.2f}")
            else:
                title.set_text(
                    f"Time = {time:.2f} \n Reached Equilibrium at Time = {eq_timestep * step_size:.2f}"
                )
            bar()
            return particles, title

        # Create animation
        ani = animation.FuncAnimation(fig, update, frames=timesteps, blit=True)

        dprint(f"saving the animation to {file_name}")
        ani.save(file_name, writer="ffmpeg", fps=30)


def plot_pair_correlation(
    r_values: NDArray[np.float64],
    g_r: NDArray[np.float64],
    file_name: str = "pair_correlation.png",
):
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


def plot_MSD(
    msd: NDArray[np.float64],
    time: NDArray[np.float64],
    eq_time: float,
    amp: float,
    exponent: float,
    file_name: str = "MSD.png",
):
    """
    Plots the mean square displacement.

    Parameters
    ----------
    msd : list
        List of mean square displacements
    time : list
        List of time passed
    eq_time: float
        The equilibrium time
    amp: float
        Amplitude of the best fit
    exponent: float
        Exponent of the best fit
    file_name : string
        Name of file to save to
    """
    print("Now plotting the Mean Square Displacement (MSD)")

    _ = plt.figure(figsize=(10, 7))
    plt.title("MSD vs Time")
    plt.xlabel(r"Time  $\left ( \sqrt{\frac{m\sigma^2}{\epsilon}}\right )$ ")
    plt.ylabel(r"MSD ($\sigma^2$)")
    plt.plot(time, msd, label="MSD")
    plt.plot(
        time,
        power_model(time - eq_time, amp, exponent),
        label=f"Fit (exp: {exponent:.2f})",
        linestyle="--",
    )
    plt.legend()
    plt.tight_layout()
    dprint(f"saving the MSD plot to {file_name}")
    plt.savefig(file_name)
    plt.close()


def print_outputs(variables: dict[str, str], simulator_type: str, export_csv: bool):
    """
    Prints the simulation results and optionally exports them to a CSV file.

    Parameters
    ----------
    variables : dict
        A dictionary containing the simulation variables and their values.
    simulator_type : str
        The type of simulator used, utilized in the CSV filename if exporting.
    export_csv : bool
        If True, exports the variables to a CSV file named `results_<simulator_type>.csv`.

    Outputs
    -------
    Prints a formatted table of the variables and their values to the console.
    """
    if export_csv:
        with open(f"results_{simulator_type}.csv", mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Variable", "Value"])
            for key, value in variables.items():
                writer.writerow([key, value])
    print_table: list[tuple[str, str]] = []
    for key, value in variables.items():
        if isinstance(value, (list, np.ndarray)):
            value_str = str(value)
            print_table.append((key, value_str))
        else:
            print_table.append((key, value))

    print(tabulate(print_table, headers=["Variable", "Value"], tablefmt="grid"))


def power_model(t, A, B):
    """
    Simple power model for the fit function to map to

    Parameters
    ----------
    t: float
        Time to fill into the model
    A: float
        Constant scaling factor in the model
    B: float
        Exponent factor in the model
    """
    return A * t**B


def r_squared(y, y_fit):
    """
    Calculates the r^2 metric for the fit y_fit relative to the original data y.

    Parameters
    ----------
    y: list[float]
        A list of datapoints the fit should be compared with
    y_fit: list[float]
        A list of estimates retrieved from the fit

    Returns
    -------
    How good the fit is in relation to the simple mean
    """
    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1 - (ss_res / ss_tot)


def best_fit(
    msd: NDArray[np.float64], t: NDArray[np.float64]
) -> tuple[str, float, float, float]:
    """
    Finds the phase by fitting a power model and then comparing the result
    to values from literature

    Parameters
    ----------
    msd: list[float]
        The mean squared displacement at time t
    t: list[float]
        The times t the msd's are at

    Returns
    -------
    phase: str
        The phase closest to the fit
    exponent: float
        The exponent that appears in the fit
    amplitude: float
        The constant amplitude factor of the fit
    r2_pow: float
        The r2 value of the fit
    """
    try:
        popt_pow, cov = opt.curve_fit(
            power_model, t, msd, bounds=([-0.5, 0.0], [np.inf, 3.0])
        )
    except:
        return "NoConverge", np.nan, np.nan, 0
    msd_fit_pow = power_model(t, *popt_pow)
    r2_pow = r_squared(msd, msd_fit_pow)

    phase = "Unknown"
    if popt_pow[1] < np.sqrt(2) / 2:
        phase = "Solid"
    elif popt_pow[1] < np.sqrt(2):
        phase = "Liquid"
    elif popt_pow[1] < 3.0:
        phase = "Gas"
    return phase, float(popt_pow[1]), float(popt_pow[0]), float(r2_pow)
