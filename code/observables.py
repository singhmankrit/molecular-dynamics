import numpy as np


def compute_pair_correlation(distances, volume, num_particles, r_max, delta_r):
    """
    Computes the time-averaged pair correlation function g(r) given a distance matrix.

    Parameters:
        distances (numpy array): A (num_particles, num_timesteps) array of pairwise distances.
        volume (float): The total volume of the simulation cell.
        num_particles (int): The total number of particles in the system.
        r_max (float): The maximum distance to consider for g(r).
        delta_r (float): The bin size for the histogram.

    Returns:
        r_values (numpy array): The midpoints of the histogram bins.
        g_r (numpy array): The time-averaged pair correlation function values.
    """
    time_steps = distances.shape[0]
    bins = np.arange(0, r_max, delta_r)
    histograms = np.zeros((time_steps, len(bins) - 1))
    # Compute histograms for all time steps
    for t in range(time_steps):
        # Get upper triangle without diagonal (pairwise distances)
        pairwise_distances = distances[t][np.triu_indices(num_particles, k=1)]
        histograms[t], _ = np.histogram(pairwise_distances, bins=bins)

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