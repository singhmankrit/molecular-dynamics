#!/usr/bin/env python
import multiprocessing
import pickle
import numpy as np
import matplotlib.pyplot as plt
from code import initialisation, sim_plots, simulators, observables

num_particles = 500
seed = 21921
timesteps = 3000
step_size = 0.005
equilibrium_steps = 25
temperature_tolerance = 0.1
equilibrium_stable_check = 3
simulator_type = "verlet"
bin_size = 0.1


def do_simulation(inputs):
    temp, vol_per_particle = inputs
    box_side_length = np.power(num_particles * vol_per_particle, 1 / 3)
    box_size = np.array([box_side_length, box_side_length, box_side_length])
    lattice_const = box_side_length / np.power(num_particles / 4, 1 / 3)
    print(
        f"am in simulation with temp {temp}, box size {box_size} and lattice constant of {lattice_const}"
    )
    corner_offset = np.array([lattice_const, lattice_const, lattice_const]) / 2

    init_pos = initialisation.fcc_lattice(
        num_particles, lattice_const, corner_offset=corner_offset
    )
    init_vel = initialisation.init_velocity(num_particles, temp, seed)
    pos, vel, kinetic, potential, virials, histograms, eq_timestep, avg_temp = (
        simulators.simulate(
            init_pos.copy(),
            init_vel.copy(),
            timesteps,
            step_size,
            box_size,
            equilibrium_steps,
            temp,
            temperature_tolerance,
            equilibrium_stable_check,
            simulator_type,
            bin_size,
            alive_params={"disable": True},
        )
    )
    print(
        f"the simulation with temp {temp}, box size {box_size} and lattice constant of {lattice_const} finished"
    )
    if eq_timestep > 0:
        msd = observables.compute_msd(pos, eq_timestep)
        time = np.arange(eq_timestep, timesteps + 1, 1) * step_size
        fit = sim_plots.best_fit(msd, time)
        print(f"the fit for {vol_per_particle}, {temp} is {fit}")
        return fit
    else:
        return ("no equilibrium", 0, 0, 0)


if __name__ == "__main__":
    temperatures = np.linspace(0.1, 20, 31)
    rhos = np.logspace(-2, 1, 31)
    vol_per_particles = 1 / rhos

    tgrid, vgrid = np.meshgrid(temperatures, vol_per_particles)
    tgrid, rhogrid = np.meshgrid(temperatures, rhos)

    with multiprocessing.Pool(20) as pool:
        fits = pool.map(do_simulation, zip(tgrid.ravel(), vgrid.ravel()))
    print(np.array(fits)[:, 0].reshape((len(temperatures), len(vol_per_particles))))
    print(np.array(fits)[:, 1].reshape((len(temperatures), len(vol_per_particles))))
    print(np.array(fits)[:, 2].reshape((len(temperatures), len(vol_per_particles))))
    print(np.array(fits)[:, 3].reshape((len(temperatures), len(vol_per_particles))))

    fits2d = np.array(fits)[:, 0].reshape((len(temperatures), len(vol_per_particles)))
    with open("phases.pkl", "wb") as file:
        pickle.dump((tgrid, rhogrid, vgrid, fits), file)

    fig, ax = plt.subplots()

    solid_idxs = np.where(fits2d == "Solid")
    ax.scatter(tgrid[solid_idxs], rhogrid[solid_idxs], c="r", label="Solid")
    liquid_idxs = np.where(fits2d == "Liquid")
    ax.scatter(tgrid[liquid_idxs], rhogrid[liquid_idxs], c="b", label="Liquid")
    gas_idxs = np.where(fits2d == "Gas")
    ax.scatter(tgrid[gas_idxs], rhogrid[gas_idxs], c="g", label="Gas")

    ax.set_xlabel("temperature")
    ax.set_ylabel("$rho$")
    ax.set_title("phase diagram of argon")
    ax.legend()

    fig.tight_layout()
    fig.savefig("phase_diagram.png")
    plt.show()
