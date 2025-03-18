import numpy as np
from utilities import dprint
from alive_progress import alive_bar
from scipy.integrate import solve_ivp


def simulate(
    init_pos,
    init_vel,
    num_tsteps,
    timestep,
    box_dim,
    equilibrium_steps,
    target_temperature,
    temperature_tolerance,
    equilibrium_stable_check,
    integrator,
    delta_r,
):
    """
    Molecular dynamics simulation based on the specified integrator

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
        The integrator based on which we want to simulate
    delta_r : float
        The bin size for the pair_correlation graph

    Returns
    -------
    Any quantities or observables that you wish to study.
    """
    amount_of_particles = len(init_pos)
    positions_list = np.zeros((num_tsteps + 1, amount_of_particles, 3))
    positions_list[0, :, :] = init_pos
    velocities_list = np.zeros((num_tsteps + 1, amount_of_particles, 3))
    velocities_list[0, :, :] = init_vel

    current_positions = init_pos
    current_velocities = init_vel
    dprint(
        f"starting positions are {current_positions} and starting velocities are {current_velocities}"
    )

    relative_positions, distances = atomic_distances(current_positions, box_dim)
    force_magnitudes, current_forces = lj_force(
        relative_positions, distances
    )  # used in verlet and leapfrog
    half_velocities = current_velocities + (
        current_forces * timestep / 2
    )  # used in leapfrog

    kinetic_energies_list = np.zeros((num_tsteps + 1))
    kinetic_energies_list[0] = kinetic_energy(init_vel)
    potential_energies_list = np.zeros((num_tsteps + 1))
    potential_energies_list[0] = potential_energy(distances)

    # For pair_correlation
    r_max = np.sqrt(box_dim[0] ** 2 + box_dim[1] ** 2 + box_dim[2] ** 2) /2
    bins = np.arange(0, r_max, delta_r)
    histograms = np.zeros((num_tsteps + 1, len(bins)-1))
    pairwise_dist = distances[np.triu_indices(amount_of_particles, k=1)]
    histograms[0], _ = np.histogram(pairwise_dist, bins=bins)


    # Counter for equilibrium stability
    stable_counter = 0
    apply_rescale = True
    # Timestep when rescaling is stopped
    equilibrium_timestep = -1
    temperature_list = np.array([])
    is_leapfrog = False

    virials = np.zeros((num_tsteps + 1))
    virials[0] = 0.5 * np.sum(distances * force_magnitudes, axis=(0, 1))

    with alive_bar(num_tsteps) as bar:
        for step in np.arange(1, num_tsteps + 1):
            dprint(
                f"""
                at step {step} the particles are at {current_positions}
                the particles have velocities {current_velocities}
                """
            )

            if integrator == "verlet":
                (
                    current_positions,
                    current_velocities,
                    current_forces,
                    distances,
                    force_magnitudes,
                ) = verlet_step(
                    current_positions,
                    current_velocities,
                    current_forces,
                    timestep,
                    box_dim,
                )
            elif integrator == "euler":
                force_magnitudes, current_forces = lj_force(
                    relative_positions, distances
                )
                current_positions, current_velocities = euler_step(
                    current_positions,
                    current_velocities,
                    current_forces,
                    timestep,
                    box_dim,
                )
                relative_positions, distances = atomic_distances(
                    current_positions, box_dim
                )
            elif integrator == "leapfrog":
                is_leapfrog = True
                (
                    current_positions,
                    half_velocities,
                    current_velocities,
                    distances,
                    force_magnitudes,
                ) = leapfrog_step(current_positions, half_velocities, timestep, box_dim)
            elif integrator == "scipy_rk45":
                (
                    current_positions,
                    current_velocities,
                    current_forces,
                    distances,
                    force_magnitudes,
                ) = scipy_rk45_step(
                    current_positions, current_velocities, timestep, box_dim
                )
            else:
                print(
                    f"Please select a valid integrator ('leapfrog', 'verlet', 'euler', 'scipy_rk45'), currently: {integrator}"
                )
                exit(4)

            current_kinetic_energy = kinetic_energy(current_velocities)

            # add the current statistics to the logs
            kinetic_energies_list[step] = current_kinetic_energy
            potential_energies_list[step] = potential_energy(distances)
            virials[step] = 0.5 * np.sum(distances * force_magnitudes, axis=(0, 1))

            # keep track of the positions and velocities
            positions_list[step, :, :] = current_positions
            velocities_list[step, :, :] = current_velocities

            # For pairwise correlation
            pairwise_dist = distances[np.triu_indices(amount_of_particles, k=1)]
            histograms[step], _ = np.histogram(pairwise_dist, bins=bins)

            current_temperature = compute_temperature(
                current_kinetic_energy, amount_of_particles
            )
            if is_leapfrog:
                # refer to half step temperature and velocities for scaling
                current_temperature = compute_temperature(kinetic_energy(half_velocities), amount_of_particles)
            # rescale velocities if applicable
            if apply_rescale == False:
                temperature_list[step-equilibrium_timestep-1] = current_temperature
            if step % equilibrium_steps == 0 and apply_rescale:
                if abs(current_temperature/target_temperature - 1) > temperature_tolerance:
                    rescale_factor = compute_rescale_factor(amount_of_particles, target_temperature, current_kinetic_energy)
                    if is_leapfrog:
                        half_velocities *= rescale_factor
                    else:
                        current_velocities *= rescale_factor
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
    return (
        positions_list,
        velocities_list,
        kinetic_energies_list,
        potential_energies_list,
        virials,
        histograms,
        equilibrium_timestep,
        average_temperature,
    )


def verlet_step(positions, velocities, forces, timestep, box_dim):
    """
    Uses the velocity-Verlet algorithm to integrate the equations of motion.

    Parameters
    ----------
    positions : np.ndarray
        The current positions of the particles in Cartesian space
    velocities : np.ndarray
        The current velocities of the particles in Cartesian space
    forces : np.ndarray
        The current forces on the particles in Cartesian space
    timestep : float
        The duration of a single simulation step
    box_dim : np.ndarray(float)
        Dimensions of the simulation box

    Returns
    -------
    new_positions : np.ndarray
        The new positions of the particles in Cartesian space
    velocities : np.ndarray
        The updated velocities of the particles in Cartesian space
    new_forces : np.ndarray
        The new forces on the particles in Cartesian space
    new_distances : np.ndarray
        The new distances between the particles
    new_magnitudes : np.ndarray
        The new forces between the particles
    """
    new_positions = (
        positions + velocities * timestep + forces * timestep * timestep / 2
    ) % box_dim
    relative_positions, new_distances = atomic_distances(new_positions, box_dim)
    new_magnitudes, new_forces = lj_force(relative_positions, new_distances)
    velocities += (forces + new_forces) * timestep / 2
    return new_positions, velocities, new_forces, new_distances, new_magnitudes


def euler_step(positions, velocities, forces, timestep, box_dim):
    """
    Uses the Euler method to integrate the equations of motion.

    Parameters
    ----------
    positions : np.ndarray
        The current positions of the particles in Cartesian space
    velocities : np.ndarray
        The current velocities of the particles in Cartesian space
    forces : np.ndarray
        The current forces acting on the particles
    timestep : float
        The duration of a single simulation step
    box_dim : np.ndarray(float)
        Dimensions of the simulation box

    Returns
    -------
    new_positions : np.ndarray
        The updated positions of the particles in Cartesian space
    velocities : np.ndarray
        The updated velocities of the particles in Cartesian space
    """
    new_positions = (positions + velocities * timestep) % box_dim
    velocities += forces * timestep
    return new_positions, velocities


def leapfrog_step(positions, half_velocities, timestep, box_dim):
    """
    Uses the Leapfrog method to integrate the equations of motion.

    Parameters
    ----------
    positions : np.ndarray
        The current positions of the particles in Cartesian space
    half_velocities : np.ndarray
        The current velocities of the particles in Cartesian space, advanced by half a time step
    timestep : float
        The duration of a single simulation step
    box_dim : np.ndarray(float)
        Dimensions of the simulation box

    Returns
    -------
    new_positions : np.ndarray
        The updated positions of the particles in Cartesian space
    half_velocities : np.ndarray
        The updated velocities of the particles in Cartesian space, advanced by half a time step
    current_velocities : np.ndarray
        The full-step velocities of the particles in Cartesian space
    distances : np.ndarray
        The distances between all particles in Cartesian space
    """
    # Position update using half-step velocities
    new_positions = (positions + half_velocities * timestep) % box_dim

    # Update the positions and forces
    relative_positions, distances = atomic_distances(new_positions, box_dim)
    forces = lj_force(relative_positions, distances)

    # Full-step velocity update
    current_velocities = half_velocities + (forces * timestep / 2)

    # Half-step velocity update for the next step
    half_velocities += forces * timestep

    return new_positions, half_velocities, current_velocities, distances


def scipy_rk45_step(positions, velocities, timestep, box_dim):
    """
    Uses the scipy.integrate.solve_ivp method to integrate the equations of motion.

    Parameters
    ----------
    positions : np.ndarray
        The current positions of the particles in Cartesian space
    velocities : np.ndarray
        The current velocities of the particles in Cartesian space
    timestep : float
        The duration of a single simulation step
    box_dim : np.ndarray(float)
        Dimensions of the simulation box

    Returns
    -------
    new_positions : np.ndarray
        The updated positions of the particles in Cartesian space
    new_velocities : np.ndarray
        The updated velocities of the particles in Cartesian space
    """
    amount_of_particles = len(positions)
    y0 = np.concatenate([positions.flatten(), velocities.flatten()])

    t_span = (0, timestep)
    t_eval = np.array([timestep])

    sol = solve_ivp(
        fun=molecular_dynamics_rhs,
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        args=(box_dim,),
        method="RK45"  # Runge-Kutta 4th/5th order
    )

    # Reshape solution: extract positions and velocities
    current_positions = sol.y[:3 * amount_of_particles].reshape(
        amount_of_particles, 3, -1).transpose(2, 0, 1).squeeze()
    current_velocities = sol.y[3 * amount_of_particles:].reshape(
        amount_of_particles, 3, -1).transpose(2, 0, 1).squeeze()
    relative_positions, distances = atomic_distances(
        current_positions, box_dim
    )
    current_forces = lj_force(relative_positions, distances)
    return current_positions, current_velocities, current_forces, distances


def molecular_dynamics_rhs(t, y, box_dim):
    """
    Compute derivatives for the molecular dynamics ODE system.

    Parameters
    ----------
    t : float
        Current time (not used explicitly, but required by solve_ivp).
    y : np.ndarray
        Flattened array containing both positions and velocities.
    box_dim : np.ndarray
        Dimensions of the simulation box.

    Returns
    -------
    dydt : np.ndarray
        Flattened derivative array containing velocity and acceleration.
    """
    amount_of_particles = len(y) // 6  # 3 for positions, 3 for velocities
    positions = y[:3 * amount_of_particles].reshape(amount_of_particles, 3)
    velocities = y[3 * amount_of_particles:].reshape(amount_of_particles, 3)

    # Compute distances and forces
    relative_positions, distances = atomic_distances(positions, box_dim)
    forces = lj_force(relative_positions, distances)

    # First derivatives (velocities)
    dpdt = velocities.flatten()

    # Second derivatives (accelerations)
    dvdt = forces.flatten()  # Assuming unit mass

    return np.concatenate([dpdt, dvdt])


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
        nxn array having the magnitudes of the forces acting on particle i due to all other particles
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

    return (force_magnitude, net_force.T)


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
