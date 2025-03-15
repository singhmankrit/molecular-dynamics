import numpy as np
from utilities import dprint
from alive_progress import alive_bar


def simulate(init_pos, init_vel, num_tsteps, timestep, box_dim, equilibrium_steps, target_temperature, temperature_tolerance, equilibrium_stable_check, integrator):
    """
    Molecular dynamics simulation using the selected integrator to
    integrate the equations of motion.

    Parameters
    ----------
    init_pos : np.ndarray
        The initial positions of the atoms in Cartesian space
    init_vel : np.ndarray
        The initial velocities of the atoms in Cartesian space
    num_tsteps : int
        The total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : np.ndarray(float)
        Dimensions of the simulation box
    equilibrium_steps : int
        Number of steps after which we apply velocity rescaling (if applicable)
    target_temperature : float
        The target temperature of the system
    temperature_tolerance : float
        The tolerated error in temperature np.abs(actual_temp/target_temp - 1)
    equilibrium_stable_check : int
        Number of stable steps after which we stop rescaling
    integrator : str
        The integrator to use (leapfrog/verlet/euler)

    Returns
    -------
    Any quantities or observables that you wish to study.
    """
    amount_of_particles = len(init_pos)
    positions = np.zeros((num_tsteps+1, amount_of_particles, 3))
    positions[0, :, :] = init_pos
    velocities = np.zeros((num_tsteps+1, amount_of_particles, 3))
    velocities[0, :, :] = init_vel

    current_positions = init_pos
    current_velocities = init_vel
    dprint(
        f"starting positions are {current_positions} and starting velocities are {current_velocities}"
    )

    # we calculate these so we can calculate the "next step" only from now on
    relative_positions, distances = atomic_distances(current_positions, box_dim)
    current_forces = lj_force(relative_positions, distances)

    kinetic_energies = np.zeros((num_tsteps+1))
    kinetic_energies[0] = kinetic_energy(init_vel)
    potential_energies = np.zeros((num_tsteps+1))
    potential_energies[0] = potential_energy(distances)
    distance_list = np.zeros((num_tsteps+1, amount_of_particles, amount_of_particles))
    distance_list[0, :, :] = distances

    # Counter for equilibrium stability
    stable_counter = 0  
    apply_rescale = True
    # Timestep when rescaling is stopped
    equilibrium_timestep = -1
    temperature_list = np.array([])

    with alive_bar(num_tsteps) as bar:
        for step in np.arange(1,num_tsteps+1):
            dprint(
                f"""
                at step {step} the particles are at {current_positions}
                the particles have velocities {current_velocities}
                """
            )

            if integrator == "verlet":
                return verlet_step(init_pos, init_vel, num_tsteps, timestep, box_dim, equilibrium_steps, target_temperature, temperature_tolerance, equilibrium_stable_check)
            if integrator == "euler":
                return euler_step(init_pos, init_vel, num_tsteps, timestep, box_dim, equilibrium_steps, target_temperature, temperature_tolerance, equilibrium_stable_check)
            if integrator == "leapfrog":
                return leapfrog_step(init_pos, init_vel, num_tsteps, timestep, box_dim, equilibrium_steps, target_temperature, temperature_tolerance, equilibrium_stable_check)
            
            current_temperature = compute_temperature(current_kinetic_energy, amount_of_particles)
            # rescale velocities if applicable
            if apply_rescale == False:
                temperature_list[step-equilibrium_timestep-1] = current_temperature
            if step % equilibrium_steps == 0 and apply_rescale:
                if abs(current_temperature/target_temperature - 1) > temperature_tolerance:
                    current_velocities *= compute_rescale_factor(amount_of_particles, target_temperature, current_kinetic_energy)
                    stable_counter = 0
                else:
                    # Increase stability count if temperature is stable
                    stable_counter += 1
                    # Exit rescaling if stable for longer than equilibrium_stable_check and store timestep
                    if (stable_counter > equilibrium_stable_check):
                        apply_rescale = False
                        equilibrium_timestep = step
                        temperature_list = np.zeros((num_tsteps - step))
            bar()
    average_temperature = np.mean(temperature_list)
    return positions, velocities, kinetic_energies, potential_energies, distance_list, equilibrium_timestep, average_temperature

def verlet(init_pos, init_vel, num_tsteps, timestep, box_dim, equilibrium_steps, target_temperature, temperature_tolerance, equilibrium_stable_check):
    """
    Molecular dynamics simulation using the Velocity Verlet's algorithm
    to integrate the equations of motion. Calculates energies and other
    observables at each timestep.

    Parameters
    ----------
    init_pos : np.ndarray
        The initial positions of the atoms in Cartesian space
    init_vel : np.ndarray
        The initial velocities of the atoms in Cartesian space
    num_tsteps : int
        The total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : np.ndarray(float)
        Dimensions of the simulation box
    equilibrium_steps : int
        Number of steps after which we apply velocity rescaling (if applicable)
    target_temperature : float
        The target temperature of the system
    temperature_tolerance : float
        The tolerated error in temperature np.abs(actual_temp/target_temp - 1)
    equilibrium_stable_check : int
        Number of stable steps after which we stop rescaling

    Returns
    -------
    Any quantities or observables that you wish to study.
    """
            # use velocity-verlet to calculate the new positions and velocities
            current_positions = (
                current_positions
                + current_velocities * timestep
                + current_forces * timestep * timestep / 2
            ) % box_dim
            relative_positions, distances = atomic_distances(current_positions, box_dim)
            new_forces = lj_force(relative_positions, distances)
            current_velocities += (current_forces + new_forces) * timestep / 2

            current_kinetic_energy = kinetic_energy(current_velocities)

            # add the current statistics to the logs
            kinetic_energies[step] = current_kinetic_energy
            potential_energies[step] = potential_energy(distances)
            distance_list[step,:,:] = distances

            # keep track of the positions and velocities
            positions[step,:,:] = current_positions
            velocities[step,:,:] = current_velocities

            # update the forces so n -> n+1
            current_forces = new_forces


def euler(init_pos, init_vel, num_tsteps, timestep, box_dim, equilibrium_steps, target_temperature, temperature_tolerance, equilibrium_stable_check):
    """
    Molecular dynamics simulation using the Euler algorithm
    to integrate the equations of motion. Calculates energies and other
    observables at each timestep.

    Parameters
    ----------
    init_pos : np.ndarray
        The initial positions of the atoms in Cartesian space
    init_vel : np.ndarray
        The initial velocities of the atoms in Cartesian space
    num_tsteps : int
        The total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : np.ndarray(float)
        Dimensions of the simulation box
    equilibrium_steps : int
        Number of steps after which we apply velocity rescaling (if applicable)
    target_temperature : float
        The target temperature of the system
    temperature_tolerance : float
        The tolerated error in temperature np.abs(actual_temp/target_temp - 1)
    equilibrium_stable_check : int
        Number of stable steps after which we stop rescaling

    Returns
    -------
    Any quantities or observables that you wish to study.
    """
            # get the n-by-3 matrix of all the total forces on the particles
            forces = lj_force(relative_positions, distances)

            current_kinetic_energy = kinetic_energy(current_velocities)

            # Euler integration step, we'll have to rewrite this to improve energy conservation
            current_positions = (
                current_positions + current_velocities * timestep
            ) % box_dim
            current_velocities += forces * timestep

            # create the n-by-n matrix of all the distances and the n-by-n-by-3 matrix of the relative positions
            relative_positions, distances = atomic_distances(current_positions, box_dim)

            # get current energies and distances and append
            kinetic_energies[step]=current_kinetic_energy
            potential_energies[step]=(potential_energy(distances))
            distance_list[step,:,:] =(distances)

            # append the new positions and velocities to the arrays
            positions[step,:,:] = current_positions
            velocities[step,:,:] = current_velocities

def leapfrog(init_pos, init_vel, num_tsteps, timestep, box_dim, equilibrium_steps, target_temperature, temperature_tolerance, equilibrium_stable_check):
    """
    Molecular dynamics simulation using the Leapfrog algorithm
    to integrate the equations of motion.

    Parameters
    ----------
    init_pos : np.ndarray
        The initial positions of the atoms in Cartesian space
    init_vel : np.ndarray
        The initial velocities of the atoms in Cartesian space
    num_tsteps : int
        The total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : np.ndarray(float)
        Dimensions of the simulation box
    equilibrium_steps : int
        Number of steps after which we apply velocity rescaling (if applicable)
    target_temperature : float
        The target temperature of the system
    temperature_tolerance : float
        The tolerated error in temperature np.abs(actual_temp/target_temp - 1)
    equilibrium_stable_check : int
        Number of stable steps after which we stop rescaling

    Returns
    -------
    Any quantities or observables that you wish to study.
    """

    # Leapfrog starts with a half-step
    half_velocities = current_velocities + (current_forces*timestep/2)

            # Position update using half-step velocities
            current_positions = (current_positions + half_velocities * timestep) % box_dim

            # Update the positions and forces
            relative_positions, distances = atomic_distances(current_positions, box_dim)
            current_forces = lj_force(relative_positions, distances)

            # Full-step velocity update
            current_velocities = half_velocities + (current_forces * timestep / 2)

            # Half-step velocity update for the next step
            half_velocities += current_forces * timestep

            # add the current statistics to the logs
            kinetic_energies[step] = kinetic_energy(current_velocities)
            potential_energies[step]=potential_energy(distances)
            distance_list[step,:,:] = distances

            # keep track of the positions and velocities
            positions[step,:,:] = current_positions
            velocities[step,:,:] = current_velocities

            # we apply rescaling at half timesteps for leapfrog
            half_kinetic_energy = kinetic_energy(half_velocities)
            half_temperature = compute_temperature(half_kinetic_energy, amount_of_particles)
            if apply_rescale == False:
                temperature_list[step-equilibrium_timestep-1] = half_temperature
            if step % equilibrium_steps == 0 and apply_rescale:
                if abs(half_temperature/target_temperature - 1) > temperature_tolerance:
                    half_velocities *= compute_rescale_factor(amount_of_particles, target_temperature, half_kinetic_energy)
                    stable_counter = 0
                else:
                    # Increase stability count if temperature is stable
                    stable_counter += 1
                    # Exit rescaling if stable for longer than equilibrium_stable_check and store timestep
                    if (stable_counter > equilibrium_stable_check):
                        apply_rescale = False
                        equilibrium_timestep = step
                        temperature_list = np.zeros((num_tsteps - step))
            bar()
    average_temperature = np.mean(temperature_list)
    return positions, velocities, kinetic_energies, potential_energies, distance_list, equilibrium_timestep, average_temperature


def compute_temperature(kinetic_energy, amount_of_particles):
    """
    Computes the temperature of the system using the kinetic energy.

    Parameters
    ----------
    kinetic_energy : float
        The total kinetic energy of the system.
    amount_of_particles : int
        The number of particles in the system.

    Returns
    -------
    float
        The computed temperature of the system.
    """
    # Equipartition theorem
    return (2 * kinetic_energy) / (3 * (amount_of_particles - 1))


def compute_rescale_factor(amount_of_particles, target_temperature, kinetic_energy):
    """
    Returns the rescale factor for the velocities of the particles to reach a target temperature.

    Parameters
    ----------
    amount_of_particles : int
        The amount of particles in the system
    target_temperature : float
        The target temperature of the system
    kinetic_energy : float
        The total kinetic energy of the system.

    Returns
    -------
    scaled_velocities : float
        The rescale factor
    """
    return np.sqrt(3*(amount_of_particles-1)*target_temperature/(2*kinetic_energy))


def atomic_distances(
    pos: np.typing.NDArray[np.float64], box_dim: np.typing.NDArray[np.float64]
) -> np.typing.NDArray[np.float64]:
    """
    Calculates relative positions and distances between particles.

    parameters
    ----------
    pos : np.ndarray
        The positions of the particles in cartesian space
    box_dim : float
        The dimension of the simulation box

    returns
    -------
    rel_pos : np.ndarray
        Relative positions of particles
    rel_dist : np.ndarray
        The distance between particles
    """
    xpos, ypos, zpos = np.unstack(pos, axis=-1)

    # create meshgrids to make an n-by-n matrix of distances
    central_x, other_x = np.meshgrid(xpos, xpos, sparse=True, copy=False)
    central_y, other_y = np.meshgrid(ypos, ypos, sparse=True, copy=False)
    central_z, other_z = np.meshgrid(zpos, zpos, sparse=True, copy=False)

    # below (3 lines) is the part that takes long but it's already optimised
    x_dist = (central_x - other_x + box_dim[0] * 0.5) % box_dim[0] - box_dim[0] * 0.5
    y_dist = (central_y - other_y + box_dim[1] * 0.5) % box_dim[1] - box_dim[1] * 0.5
    z_dist = (central_z - other_z + box_dim[2] * 0.5) % box_dim[2] - box_dim[2] * 0.5

    relative_positions = np.stack((x_dist, y_dist, z_dist))
    distances = np.ma.masked_values(
        np.sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist),
        0.0,
        rtol=1e-60,
        atol=1e-60,
    )
    dprint(f"there are {np.ma.count_masked(distances)} masked distance values")

    # return the full distance matrix of shape n-by-n
    return (relative_positions, distances)


def lj_force(rel_pos, rel_dist):  # units of epsilon/sigma
    """
    Calculates the net forces on each atom from the matrices containing the positions and distances.

    Parameters
    ----------
    rel_pos : np.ndarray
        Relative particle positions as obtained from atomic_distances
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    np.ndarray
        nx3 array having the net vector force acting on particle i due to all other particles
    """
    # Compute force magnitude using the Lennard-Jones force formula
    force_magnitude = 24 * (1 / rel_dist) ** 7 - 48 * (1 / rel_dist) ** 13
    dprint(
        f"the minimal force magnitude is {np.min(force_magnitude)} and the maximum force is {np.max(force_magnitude)}"
    )

    # Compute force matrix
    force_direction = rel_pos / rel_dist

    force_matrix = force_magnitude * force_direction

    # Sum forces acting on each particle
    net_force = -np.sum(force_matrix, axis=1)

    return net_force.T


def kinetic_energy(vel):  # units of epsilon
    """
    Computes the kinetic energy of an atomic system.

    Parameters
    ----------
    vel: np.ndarray
        Velocity of particle

    Returns
    -------
    float
        The total kinetic energy of the system.
    """
    # Kinetic energy for each particle
    ke_individual = 1 / 2 * np.sum(vel**2, axis=1)
    # Total Kinetic Energy
    ke = np.sum(ke_individual)
    return ke


def lj_potential(distance):  # units of epsilon
    return 4 * ((1 / distance) ** 12 - (1 / distance) ** 6)


def potential_energy(rel_dist):
    """
    Computes the potential energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    float
        The total potential energy of the system.
    """

    return 1 / 2 * np.sum(lj_potential(rel_dist))
