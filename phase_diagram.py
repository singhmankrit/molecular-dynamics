#!/usr/bin/env python
import multiprocessing
from time import sleep
import numpy as np
from code import initialisation, sim_plots, simulators, utilities, observables


num_particles = 500
seed = 21921
timesteps = 2000
step_size = 0.005
equilibrium_steps = 25
temperature_tolerance = 0.1
equilibrium_stable_check = 3
simulator_type = "verlet"
bin_size = 0.1


def do_simulation(inputs):
    print("am in simulation")
    vol_per_particle, temp = inputs
    box_side_length = np.power(num_particles * vol_per_particle, 1 / 3)
    box_size = np.array([box_side_length, box_side_length, box_side_length])
    lattice_const = box_side_length / np.power(num_particles / 4, 1 / 3)
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
            alive_params={"force_tty": False},
        )
    )
    print(f"The simulation for {vol_per_particle}, {temp} finished")


if __name__ == "__main__":
    temperatures = np.linspace(0.1, 30, 40)
    vol_per_particles = np.logspace(0, 8, 20)

    tgrid, vgrid = np.meshgrid(temperatures, vol_per_particles)

    with multiprocessing.Pool(4) as pool:
        print("am in pool")
        pool.map(do_simulation, zip(tgrid.ravel(), vgrid.ravel()))

    print("dispatched all, waiting")
    sleep(1000)
