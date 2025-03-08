import json
from typing import Any, TypeIs
import numpy as np
import os

from numpy.typing import NDArray

debug = True if os.environ.get("DEBUG") is not None else False


def dprint(str: object):
    """
    Prints the passed string only if the debug environment variable is set.

    Parameters
    ----------
    str : string
        The string to print
    """
    if debug:
        print(str)


def parse_config(file_path: str):
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
        config: dict[str, Any] = json.load(file)
        amount_of_particles: int = config.get("particles", 2)
        step_size: float = config.get("step_size", 0.01)
        time_steps: int = config.get("time_steps", 1000)
        equilibrium_steps = config.get("equilibrium_steps", 10)
        temperature: float = config.get("temperature", 1)
        temperature_tolerance: float = config.get("temperature_tolerance", 0.01)
        equilibrium_stable_check: int = config.get("equilibrium_stable_check", 3)
        box: dict[str, float] | None = config.get("box")
        if box is not None:
            box_x: float = config["box"].get("x", 5.0)
            box_y: float = config["box"].get("y", 5.0)
            box_z: float = config["box"].get("z", 5.0)
        else:
            box_x, box_y, box_z = 5.0, 5.0, 5.0
        random_seed: int = config.get("seed", np.random.randint(0, 1000000000))
        pos_method: str = config.get("position_method", "uniform")
        vel_method: str = config.get("velocity_method", "mbdist")
        simulator_type: list[str] = config.get("simulator_type", ["verlet"])
        outputs: list[str] = config.get("outputs", ["energies", "distances", "animation","pair_correlation","MSD"])
        enable_cache: bool = config.get("do_caching", True)
        lat_const: float = config.get("lattice_const", 1.5)
        corner_offset: list[float] = config.get("corner_offset", [0, 0, 0])
        bin_size: float = config.get("bin_size", 0.1)
        return (
            amount_of_particles,
            step_size,
            time_steps,
            equilibrium_steps,
            temperature,
            temperature_tolerance,
            equilibrium_stable_check,
            np.array([box_x, box_y, box_z]),
            random_seed,
            pos_method,
            vel_method,
            simulator_type,
            outputs,
            enable_cache,
            lat_const,
            corner_offset,
            bin_size,
        )


def is1d(
    arr: NDArray[np.float64],
) -> TypeIs[np.ndarray[tuple[int], np.dtype[np.float64]]]:
    return len(arr.shape) == 1


def is2d(
    arr: NDArray[np.float64],
) -> TypeIs[np.ndarray[tuple[int, int], np.dtype[np.float64]]]:
    return len(arr.shape) == 2


def is3d(
    arr: NDArray[np.float64],
) -> TypeIs[np.ndarray[tuple[int, int, int], np.dtype[np.float64]]]:
    return len(arr.shape) == 3
