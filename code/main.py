#!/usr/bin/env python


import initialisation
import simulators
import sim_plots
import utilities


# Read the configuration file
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
) = utilities.parse_config("config.json")
print(
    f"simulating {amount_of_particles} particles for {timesteps} timesteps, with a time step size of {step_size}"
)

# generate the initial positions of the particles using the specified method
init_pos = None
if position_init_method == "uniform":
    init_pos = initialisation.uniform_random(amount_of_particles, box_size, random_seed)
elif position_init_method == "static":
    init_pos = initialisation.static(amount_of_particles, box_size)
else:
    print(
        f"Please select a valid position init method ('uniform', 'static'), currently: {position_init_method}"
    )
    exit(2)

# generate the initial velocities of the particles using the specified method
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

# select the simulator to use
simulator = None
if simulator_type == "verlet":
    simulator = simulators.verlet
elif simulator_type == "euler":
    simulator = simulators.euler
elif simulator_type == "leapfrog":
    simulator = simulators.leapfrog
else:
    print(
        f"Please select a valid simulator ('leapfrog', 'verlet', 'euler'), currently: {simulator_type}"
    )
    exit(4)

# run the simulation and put the results into variables
print("Starting simulation")
pos, vel, kinetic, potential, distance_list = simulator(
    init_pos,
    init_vel,
    timesteps,
    step_size,
    box_size,
)
print("Finished simulation")

# generate plots from the results
if "energies" in plots:
    sim_plots.plot_energy(kinetic, potential)
if "distances" in plots:
    sim_plots.plot_distances(distance_list)
if "animation" in plots:
    sim_plots.create_animation(pos, timesteps, box_size)
