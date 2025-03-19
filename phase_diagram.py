#!/usr/bin/env python
import numpy as np
from code import initialisation, sim_plots, simulators, utilities, observables


temperatures = np.linspace(0.1, 30, 40)
vol_per_particles = np.logspace(0, 8, 20)
num_particles = 500
seed = 21921
timesteps = 2000
step_size = 0.005
equilibrium_steps = 25
temperature_tolerance = 0.1
equilibrium_stable_check = 3
simulator_type = "verlet"
bin_size = 0.1

for vol_per_particle in vol_per_particles:
    for temp in temperatures:
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
            )
        )
