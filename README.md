# Project 1: Molecular dynamics simulation of Argon atoms

A numpy based simulator for the molecular dynamics of Argon, has multiple integrators which
can be chosen using the config file `config.json`.

Created by: Mankrit, Noa and Ragi.

## Config file layout

The configuration file uses the json format. The options are listed below with their default value and meaning.

| Option | Default | Description |
| ------ | ------- | ----------- |
| particles | `3` | The amount of particles to simulate |
| box | `[5,5,5]` | The size of the periodic boundary conditions in dimensionless units, uses the format `[x,y,z]` |
| lattice_constant | `1.5` | The lattice constant of the fcc lattice, only used when using `"fcc"` initial positions |
| corner_offset | `[0,0,0]` | The amount to offset the corner of the fcc lattice relative to origin, only used then using `"fcc"` initial positions |
| time_steps | `1000` | The amount of timesteps to simulate |
| step_size | `0.01` | The size of a single timestep in dimensionless units |
| temperature | `0.3` | The dimensionless temperature of the system |
| position_method | `"uniform"` | The way the starting positions of the particles get generated, can be `"uniform"`, `"fcc"` or `"static"`.
| velocity_method | `"zero"` | The way the starting velocities of the particles get generated, can be `"zero"` for zero velocities or `"mbdist"` for a maxwell-boltzmann distribution.
| seed | random | The seed to use for the initial random generation.
| simulator_type | `["verlet", "euler"]` | What integrators to use on the starting position, the simulator executes every one in the list, can be `"verlet"`, `"euler"` or `"leapfrog"`.
| plots | `["energies", "distances", "animation"]` | The outputs the simulator should create from the simulation results.
| do_caching | `true` | The simulator will cache simulations in `./cache/` if this is `true`, changing the code does **not** wipe this cache, changing relevant configuration sections should create their own entries.


