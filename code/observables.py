import numpy as np


def compute_pair_correlation(box_size, histograms, num_particles, delta_r):
    """
    Computes the time-averaged pair correlation function g(r) given a distance matrix.

    Parameters:
        box_size (tuple): The size of the simulation box.
        histograms (np.ndarray): A list of distance histograms for each timestep.
        num_particles (int): The total number of particles in the system.
        delta_r (float): The bin size for the histogram.

    Returns:
        r_values (numpy array): The midpoints of the histogram bins.
        g_r (numpy array): The time-averaged pair correlation function values.
    """
    volume = box_size[0] * box_size[1] * box_size[2]
    r_max = np.sqrt(box_size[0]**2 * box_size[1]**2 * box_size[2]**2)/2
    bins = np.arange(0, r_max, delta_r)

    # Time-averaged histogram
    avg_hist = np.mean(histograms, axis=0)[:, np.newaxis]
    # Compute g(r)
    r_values = (bins[:-1] + delta_r / 2)[:, np.newaxis]
    shell_volumes = 4 * np.pi * r_values**2 * delta_r
    norm_factor = (2 * volume) / (num_particles * (num_particles - 1))
    g_r = norm_factor * avg_hist / shell_volumes
    return r_values, g_r


def compute_msd(pos, t_eq):
    """
    Computes the mean squared displacement (MSD) given a 3D array of
    particle positions, starting from a given equilibrium time step.

    Parameters:
        pos (numpy array): A 3D array of shape (time_steps, num_particles, 3)
            containing the positions of particles at each time step.
        t_eq (int): The time step at which to start computing the MSD.

    Returns:
        msd (numpy array): A 1D array of shape (time_steps - t_eq,)
            containing the mean squared displacement at each time step.
    """

    pos_eq = pos[t_eq:]  # shape: (time_steps - t_eq, num_particles, 3)
    pos_0 = pos_eq[0]    # reference positions at t_eq

    # Compute displacements
    displacements = pos_eq - pos_0  # (time_steps - t_eq, num_particles, 3)

    # Square displacements and sum over x, y, z
    # (time_steps - t_eq, num_particles)
    squared_displacements = np.sum(displacements**2, axis=2)

    # Average over particles
    msd = np.mean(squared_displacements, axis=1)  # (time_steps - t_eq,)
    return msd


def compute_compressibility_factor(temperature, amount_of_particles, virials, eq_timestep):
    """
    Computes the compressibility factor for the given temperature and amount of particles.

    Parameters:
        temperature (float): The temperature of the system.
        amount_of_particles (int): The number of particles in the system.
        virials (numpy array): A 1D array of the virial at each timestep.

    Returns:
        float: The compressibility factor.
    """

    beta = 1.0/temperature
    average_virial = np.mean(virials[eq_timestep:])
    compressibility_factor = (1-beta/(3*amount_of_particles)*average_virial)
    return compressibility_factor


def compute_specific_heat(kinetic_energy, amount_of_particles, eq_timestep):
    """
    Computes the specific heat for the given kinetic energy and amount of particles.

    Parameters:
        kinetic_energy (float): The total kinetic energy of the system.
        amount_of_particles (int): The number of particles in the system.

    Returns:
        float: The specific heat.
    """

    # Consider only timesteps after equilibrium
    ke_eq = np.array(kinetic_energy[eq_timestep:])

    # Calculate mean and variance of kinetic energy
    K_mean = np.mean(ke_eq)
    K_sq_mean = np.mean(ke_eq**2)
    delta_K_sq = K_sq_mean - K_mean**2

    # Compute specific heat
    c_V = 1/(2/(3*amount_of_particles) - (delta_K_sq/(K_mean**2)))

    return c_V/amount_of_particles

def bootstrap_specific_heat(kinetic_energy, N, num_samples=10000, block_size=100):
    """
    Computes the error in the specific heat for the given kinetic energy and amount of particles.

    Parameters:
        kinetic_energy (float): The total kinetic energy of the system.
        amount_of_particles (int): The number of particles in the system.

    Returns:
        float: The specific heat.
    """
    
    kinetic_energy = np.array(kinetic_energy)
    num_blocks = len(kinetic_energy) // block_size
    block_averages_1 = np.array([np.mean(
        kinetic_energy[i * block_size: (i + 1) * block_size]) for i in range(num_blocks)])

    block_averages_2 = np.array([np.mean(
        kinetic_energy[i * block_size: (i + 1) * block_size]**2) for i in range(num_blocks)])

    bootstrap_values = []
    
    for _ in range(num_samples):
        # check seed stuff
        bootstrap_indices = np.random.randint(0, len(block_averages_1), size=num_blocks)
        
        # Compute bootstrap averages
        K_mean_bs = np.mean(block_averages_1[bootstrap_indices])
        K2_mean_bs = np.mean(block_averages_2[bootstrap_indices])
        delta_K_sq = K2_mean_bs - K_mean_bs**2

        # Compute specific heat
        c_V = 3/(2*(1-(3*N*delta_K_sq)/(2*K_mean_bs**2)))
        bootstrap_values.append(c_V)

    cV_mean = np.mean(bootstrap_values)
    cV_error = np.sqrt(
        np.mean(np.square(bootstrap_values)) - cV_mean**2)

    return cV_mean, cV_error


