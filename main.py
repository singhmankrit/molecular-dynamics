#!/usr/bin/env python


import hashlib
import pickle
from os.path import isfile

import numpy as np
from code import initialisation, sim_plots, simulators, utilities, observables

# Read the configuration file
(
    amount_of_particles,
    step_size,
    timesteps,
    equilibrium_steps,
    temperature,
    temperature_tolerance,
    equilibrium_stable_check,
    box_size,
    random_seed,
    position_init_method,
    velocity_init_method,
    simulator_types,
    outputs,
    enable_cache,
    lat_const,
    corner_offset,
    bin_size,
    export_csv,
) = utilities.parse_config("config.json")
utilities.dprint(
    f"simulating {amount_of_particles} particles for {timesteps} timesteps, with a time step size of {step_size}"
)

if not utilities.is1d(box_size):
    raise TypeError("box_size is not 1D")
# generate the initial positions of the particles using the specified method
init_pos = None
if position_init_method == "uniform":
    init_pos = initialisation.uniform_random(amount_of_particles, box_size, random_seed)
elif position_init_method == "static":
    init_pos = initialisation.static(amount_of_particles, box_size)
elif position_init_method == "fcc":
    cornernp = np.array(corner_offset)
    if not utilities.is1d(cornernp):
        raise TypeError("corner_offet is not 1D")
    init_pos = initialisation.fcc_lattice(
        amount_of_particles, lat_const, cornernp
    )  # TODO: lattice constant (array??)
else:
    print(
        f"Please select a valid position init method ('uniform', 'static', 'fcc'), currently: {position_init_method}"
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

for simulator_type in simulator_types:
    pos, vel, kinetic, potential, distance_list = None, None, None, None, None
    # check if we have a cached run already
    hash = hashlib.sha1(
        "{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}".format(
            amount_of_particles,
            step_size,
            timesteps,
            equilibrium_steps,
            temperature,
            temperature_tolerance,
            equilibrium_stable_check,
            list(box_size),
            random_seed,
            position_init_method,
            velocity_init_method,
            simulator_type,
            lat_const,
            corner_offset,
            bin_size,
        ).encode()
    )
    path = f"cache/{hash.hexdigest()}.pkl"
    if isfile(path) and enable_cache:
        print(f"found cached results at {path}")
        with open(path, "rb") as f:
            pos, vel, kinetic, potential, virials, histograms, eq_timestep, avg_temp = (
                pickle.load(f)
            )
    else:
        # run the simulation and put the results into variables
        print(f"Starting {simulator_type} simulation")
        pos, vel, kinetic, potential, virials, histograms, eq_timestep, avg_temp = (
            simulators.simulate(
                init_pos.copy(),
                init_vel.copy(),
                timesteps,
                step_size,
                box_size,
                equilibrium_steps,
                temperature,
                temperature_tolerance,
                equilibrium_stable_check,
                simulator_type,
                bin_size,
            )
        )
        if eq_timestep == -1:
            print(f"Equilibrium not reached in simulation")
        print(f"Finished {simulator_type} simulation")

        if enable_cache:
            with open(path, "wb") as f:
                pickle.dump(
                    (
                        pos,
                        vel,
                        kinetic,
                        potential,
                        virials,
                        histograms,
                        eq_timestep,
                        avg_temp,
                    ),
                    f,
                )
    temp_error = np.abs(avg_temp / temperature - 1) * 100
    variables = {
        "Number of Particles": amount_of_particles,
        "Step Size": step_size,
        "Number of Timesteps": timesteps,
        "Box Size": box_size,
        "Lattice Constant": lat_const,
        "Position Initialisation": position_init_method,
        "Velocity Initialisation": velocity_init_method,
        "Simulator Type": simulator_type,
        "Target Temperature": temperature,
        "Average Temperature (after equilibrium)": f"{avg_temp:.02f}",
        "Temperature Error": f"{temp_error:.02f}%",
        "Equilibrium achieved at Timestep": eq_timestep,
    }

    # generate plots from the results
    if "energies" in outputs:
        sim_plots.plot_energy(
            kinetic,
            potential,
            step_size,
            timesteps,
            eq_timestep,
            file_name=f"energies_{simulator_type}.png",
        )
    if "distances" in outputs:
        print("distance output is deprecated, please remove it from your config")
    if "animation" in outputs:
        sim_plots.create_animation(
            pos,
            timesteps,
            step_size,
            eq_timestep,
            box_size,
            file_name=f"particles_{simulator_type}.mp4",
        )

    if eq_timestep < 0:
        print(f"simulation did not converge, not computing observables")
        sim_plots.print_outputs(variables, simulator_type, export_csv)
        continue

    if "pair_correlation" in outputs:
        # Compute pair correlation
        r_values, g_r = observables.compute_pair_correlation(
            box_size, histograms, amount_of_particles, bin_size
        )
        sim_plots.plot_pair_correlation(
            r_values, g_r, file_name=f"pair_correlation_{simulator_type}.png"
        )

    if "MSD" in outputs:
        msd = observables.compute_msd(pos, eq_timestep)
        time = np.arange(eq_timestep, timesteps + 1, 1) * step_size
        eq_time = eq_timestep * step_size
        state, exponent, amp, r2_pow = sim_plots.best_fit(msd, time)
        sim_plots.plot_MSD(
            msd, time, amp, exponent, file_name=f"MSD_{simulator_type}.png"
        )
        variables["State of Matter"] = state

    if "compressibility" in outputs:
        compressibility = observables.compute_compressibility_factor(
            temperature, amount_of_particles, virials, eq_timestep
        )
        variables["Compressibility Factor"] = compressibility.round(4)

    if "specific_heat" in outputs:
        specific_heat = observables.compute_specific_heat(
            kinetic, amount_of_particles, eq_timestep
        )
        variables["Specific Heat"] = specific_heat.round(4)
        cV_mean, cV_error = observables.bootstrap_specific_heat(
            kinetic[eq_timestep + 1 :], amount_of_particles
        )
        variables["Specific Heat (Bootstrap)"] = (
            f"{cV_mean.round(4)} ± {cV_error.round(4)}"
        )

    # Program Outputs
    sim_plots.print_outputs(variables, simulator_type, export_csv)
