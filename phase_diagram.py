#!/usr/bin/env python
import multiprocessing
from os.path import isfile
import pickle
import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt
from code import initialisation, sim_plots, simulators, observables
import tqdm.contrib.concurrent as conc

num_particles = 500
seed = 21921
np.random.seed(seed)
timesteps = 4000
step_size = 0.005
equilibrium_steps = 25
temperature_tolerance = 0.1
equilibrium_stable_check = 3
simulator_type = "verlet"
bin_size = 0.1


def save_to_cache(temp, lattice_const, to_save):
    with open(f"cache/phase_{temp}_{lattice_const}.pkl", "wb") as file:
        pickle.dump(to_save, file)


def load_from_cache(temp, lattice_const):
    path = f"cache/phase_{temp}_{lattice_const}.pkl"
    if not isfile(path):
        return None
    with open(path, "rb") as file:
        try:
            return pickle.load(file)
        except:
            print(f"could not load {path}, due to {repr(sys.exception())}")
            return None


def do_simulation(temp, vol_per_particle):
    box_side_length = np.power(num_particles * vol_per_particle, 1 / 3)
    box_size = np.array([box_side_length, box_side_length, box_side_length])
    lattice_const = box_side_length / np.power(num_particles / 4, 1 / 3)
    corner_offset = np.array([lattice_const, lattice_const, lattice_const]) / 2
    cached = load_from_cache(temp, lattice_const)
    pos, eq_timestep, kinetic, potential = None, None, None, None
    if cached is not None:
        pos, eq_timestep, kinetic, potential = cached
    else:
        init_pos = initialisation.fcc_lattice(
            num_particles, lattice_const, corner_offset=corner_offset
        )
        init_vel = initialisation.init_velocity(
            num_particles, temp, np.random.randint(0, 100000)
        )
        pos, _, kinetic, potential, _, _, eq_timestep, _ = simulators.simulate(
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
        try:
            save_to_cache(temp, lattice_const, (pos, eq_timestep, kinetic, potential))
        except:
            print(repr(sys.exception()))
    if eq_timestep < 0:
        return ("NoEquilibrium", np.nan, np.nan, 0)
    msd = observables.compute_msd(pos, eq_timestep)
    if (
        np.abs(
            (kinetic[eq_timestep] + potential[eq_timestep])
            / (kinetic[-1] + potential[-1])
            - 1
        )
        > 0.1
    ):
        return ("Explosion", np.nan, np.nan, 0)
    time = np.arange(eq_timestep, timesteps + 1, 1) * step_size
    fit = sim_plots.best_fit(msd, time)
    return fit


if __name__ == "__main__":
    MIN_TEMP_EXP = -2.5
    MAX_TEMP_EXP = 0.5
    TEMP_POINTS_CALC = 41
    TEMP_POINTS_DISPLAY = 400
    MIN_PRESSURE_EXP = -2
    MAX_PRESSURE_EXP = 1
    PRESSURE_POINTS_CALC = 41
    PRESSURE_POINTS_DISPLAY = 400

    temperatures = np.logspace(MIN_TEMP_EXP, MAX_TEMP_EXP, TEMP_POINTS_CALC)
    pressures = np.logspace(MIN_PRESSURE_EXP, MAX_PRESSURE_EXP, PRESSURE_POINTS_CALC)

    tgrid, pressure_grid = np.meshgrid(temperatures, pressures)
    rhogrid = pressure_grid / tgrid
    vgrid = 1 / rhogrid

    fits = conc.process_map(
        do_simulation, tgrid.ravel(), vgrid.ravel(), max_workers=20, chunksize=1
    )

    with open("phases.pkl", "wb") as file:
        pickle.dump((tgrid, rhogrid, vgrid, fits), file)

    pgrid = rhogrid * tgrid
    print(fits)
    fits2d = np.array(fits)[:, 0].reshape((len(temperatures), len(pressures)))
    pows2d = (
        np.array(fits)[:, 1]
        .reshape((len(temperatures), len(pressures)))
        .astype(np.float64)
    )
    mask = np.isnan(pows2d)
    tgrid_flat = tgrid[~mask]
    pgrid_flat = pgrid[~mask]
    pows_flat = pows2d[~mask]

    disp_temps = np.logspace(MIN_TEMP_EXP, MAX_TEMP_EXP, TEMP_POINTS_DISPLAY)
    disp_pres = np.logspace(MIN_PRESSURE_EXP, MAX_PRESSURE_EXP, PRESSURE_POINTS_DISPLAY)
    dtgrid, dpgrid = np.meshgrid(disp_temps, disp_pres)
    pows_interpolated = scipy.interpolate.griddata(
        (tgrid_flat, pgrid_flat), pows_flat, (dtgrid, dpgrid), method="linear"
    )

    fig, ax = plt.subplots()
    contour = ax.contourf(
        dtgrid,
        dpgrid,
        pows_interpolated,
        cmap="managua",
        levels=[-0.5, np.sqrt(2) / 2, np.sqrt(2), 2.5],
    )

    ax.set_xlabel("temperature")
    ax.set_ylabel(r"$\rho \dot t$")
    ax.set_yscale("log")
    ax.set_title("phase diagram of argon")

    fig.colorbar(contour)
    fig.tight_layout()
    fig.savefig("phase_diagram_contour.png")

    fig, ax = plt.subplots()
    solid_idxs = np.where(fits2d == "Solid")
    ax.scatter(
        tgrid[solid_idxs],
        pgrid[solid_idxs],
        c="r",
        label="Solid",
        marker="s",
    )
    liquid_idxs = np.where(fits2d == "Liquid")
    ax.scatter(
        tgrid[liquid_idxs], pgrid[liquid_idxs], c="b", label="Liquid", marker="o"
    )
    gas_idxs = np.where(fits2d == "Gas")
    ax.scatter(tgrid[gas_idxs], pgrid[gas_idxs], c="g", label="Gas", marker="x")

    ax.set_xlabel("temperature")
    ax.set_ylabel(r"$\rho \dot t$")
    ax.set_yscale("log")
    ax.set_title("phase diagram of argon")

    ax.legend()
    fig.tight_layout()
    fig.savefig("phase_diagram_points.png")
