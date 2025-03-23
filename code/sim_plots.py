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
    msd: NDArray[np.float64], time: NDArray[np.float64], file_name: str = "MSD.png"
):
    """
    Plots the mean square displacement.

    Parameters
    ----------
    msd : list
        List of mean square displacements
    time : list
        List of time passed
    file_name : string
        Name of file to save to
    """
    print("Now plotting the Mean Square Displacement (MSD)")

    _ = plt.figure(figsize=(10, 7))
    plt.title("MSD vs Time")
    plt.xlabel(r"Time  $\left ( \sqrt{\frac{m\sigma^2}{\epsilon}}\right )$ ")
    plt.ylabel(r"MSD ($\sigma^2$)")
    plt.plot(time, msd, label="MSD")
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


def linear_model(t: float, D: float) -> float:
    """
    A linear fit model for fitting MSD, fits best when the MSD is Diffusive (liquid)

    Parameters
    ----------
    t: float
        The time to evaluate at
    D: float
        The slope of the linear model

    Returns
    -------
    The model output at the time t
    """
    return D * t


def quadratic_model(t: float, A: float) -> float:  # Ballistic (gas)
    """
    A quadratic fit model for fitting MSD, fits best when the MSD is Ballistic (gas)

    Parameters
    ----------
    t: float
        The time to evaluate at
    A: float
        The scaling factor for the quadratic model

    Returns
    -------
    The model output at the time t
    """
    return A * t**2


def constant_model(t: float, C: float, D: float) -> float:  # Localized (solid)
    """
    An exponential decay to constant fit model for fitting MSD, fits best when the MSD is Localized (solid)

    Parameters
    ----------
    t: float
        The time to evaluate at
    A: float
        The scaling factor for the quadratic model

    Returns
    -------
    The model output at the time t
    """
    return C * (1 - np.exp(-D * t))


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


def best_fit(msd, t):
    popt_lin, _ = opt.curve_fit(linear_model, t, msd)
    popt_quad, _ = opt.curve_fit(quadratic_model, t, msd)
    msd_fit_lin = linear_model(t, *popt_lin)
    msd_fit_quad = quadratic_model(t, *popt_quad)

    r2_lin = r_squared(msd, msd_fit_lin)
    r2_quad = r_squared(msd, msd_fit_quad)
    try:
        popt_const, _ = opt.curve_fit(constant_model, t, msd, p0=[2,0.1], bounds=((0,0),(np.inf,1)))
        msd_fit_const = constant_model(t, *popt_const)
        r2_const = r_squared(msd, msd_fit_const)
    except:
        r2_const = -np.inf
    # # Print results
    # print(f"Linear Fit (Liquid): D = {popt_lin[0]:.5f}, R² = {r2_lin:.5f}")
    # print(f"Quadratic Fit (Gas): A = {popt_quad[0]:.5f}, R² = {r2_quad:.5f}")
    # print(
    #     f"Constant Fit (Solid): C = {popt_const[0]:.5f}, R² = {r2_const:.5f}")

    # Determine best fit
    best_fit = max((r2_lin, "Liquid"), (r2_quad, "Gas"),
                   (r2_const, "Solid"))[1]
    # print(f"Best fit suggests the system behaves as a {best_fit}")

    # # Plot results
    # plt.scatter(t, msd, label="MSD Data", color="black", s=10)
    # plt.plot(t, msd_fit_lin,
    #          label=f"Linear Fit (R²={r2_lin:.2f})", linestyle="--")
    # plt.plot(t, msd_fit_quad,
    #          label=f"Quadratic Fit (R²={r2_quad:.2f})", linestyle=":")
    # plt.plot(t, msd_fit_const,
    #          label=f"Constant Fit (R²={r2_const:.2f})", linestyle="-.")
    # plt.xlabel("Time")
    # plt.ylabel("MSD")
    # plt.legend()
    # plt.savefig("best_fit.png")
    # plt.close()

    return best_fit, r2_lin, r2_quad, r2_const
