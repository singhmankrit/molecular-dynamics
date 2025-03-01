#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import initialisation
import simulators
import sim_plots
from utilities import parse_config, dprint


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
        plots,
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
    
    print("Starting simulation")
    pos, vel, kinetic, potential, distance_list = simulator(
        init_pos,
        init_vel,
        timesteps,
        step_size,
        box_size,
    )
    print("Finished simulation")
    
    if "energies" in plots:
        sim_plots.plot_energy(kinetic, potential)
    if "distances" in plots:
        sim_plots.plot_distances(distance_list)
    if "animation" in plots:
        sim_plots.create_animation(pos, timesteps, box_size)
