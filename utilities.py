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
        amount_of_particles = config.get("particles", 2)
        step_size = config.get("step_size", 0.01)
        time_steps = config.get("time_steps", 1000)
        temperature = config.get("temperature", 1)
        box = config.get("box")
        if box is not None:
            box_x = config["box"].get("x", 5.0)
            box_y = config["box"].get("y", 5.0)
            box_z = config["box"].get("z", 5.0)
        else:
            box_x, box_y, box_z = 5.0, 5.0, 5.0
        random_seed = config.get("seed", None)
        pos_method = config.get("position_method", "uniform")
        vel_method = config.get("velocity_method", "mbdist")
        return (
            amount_of_particles,
            step_size,
            time_steps,
            temperature,
            np.array([box_x, box_y, box_z]),
            random_seed,
            pos_method,
            vel_method,
        )
