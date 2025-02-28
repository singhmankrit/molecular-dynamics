import json
import numpy as np


def parse_config(file_path):
    """
    Parse configuration located at `file_path`,
    return the parsed result

    Parameters
    ----------
    file_path: string
        the path to the config file

    Returns
    -------
    A tuple containing:
    - amount of particles: integer
    - time step size: float
    - time step amount: integer
    - temperature: float
    - the box size: np.ndarray
    """
    with open(file_path) as file:
        config = json.load(file)
        amount_of_particles = config["particles"]
        step_size = config["step_size"]
        time_steps = config["time_steps"]
        temperature = config["temperature"]
        box_x = config["box"]["x"]
        box_y = config["box"]["y"]
        box_z = config["box"]["z"]
        random_seed = config.get("seed", None)
        return (
            amount_of_particles,
            step_size,
            time_steps,
            temperature,
            np.array([box_x, box_y, box_z]),
            random_seed,
        )
