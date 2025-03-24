import json
from typing import Any, TypeIs
import numpy as np
import os

from numpy.typing import NDArray

debug = True if os.environ.get("DEBUG") is not None else False


def fcc_for_box(
    box_size: NDArray[np.float64], num_particles: int
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Calculates the lattice constant such that the crystal repeats exactly in the box
    """
    lattice_const = box_size / np.power(num_particles / 4, 1 / 3)
    corner_offset = lattice_const / 2
    return lattice_const, corner_offset


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
    amount_of_particles: int
        The amount of particles to simulate
    step_size: float
        The timestep size
    time_steps: int
        How many timesteps to simulate
    equilibrium_steps: int
        How often to rescale towards the target temperature
    temperature: float
        The target temperature
    temperature_tolerance: float
        A relative tolerance on how close to be to the target temperature
    equilibrium_stable_check:
        How many equilibrium checks the temperature needs to be within the tolerance before
        declaring equilibrium
    box_size: np.array([box_x, box_y, box_z])
        The dimensions of the box
    random_seed: int | None
        The random seed to use for the simulation, or random if None
    pos_method: str
        The initialisation method for the particle positions
    vel_method: str
        The initialisation method for the particle velocities
    simulator_type: list[str]
        What simulators to simulate the initial state with
    outputs: list[str]
        What outputs to produce for each simulation
    enable_cache: bool
        Whether to use the cache or not
    lat_const: float
        The lattice constant for fcc lattice
    corner_offset: np.array([offset_x, offset_y, offset_z])
        How far to offset the corner of the fcc lattice from the box corner
    bin_size: int
        How big to make the bins for pair correlation
    export_csv: bool
        Whether to export the observable stats to csv or only print them
    """
    with open(file_path) as file:
        config: dict[str, Any] = json.load(file)
        amount_of_particles: int = config.get("particles", 108)
        step_size: float = config.get("step_size", 0.01)
        time_steps: int = config.get("time_steps", 1000)
        equilibrium_steps = config.get("equilibrium_steps", 25)
        temperature: float = config.get("temperature", 1)
        temperature_tolerance: float = config.get("temperature_tolerance", 0.2)
        equilibrium_stable_check: int = config.get("equilibrium_stable_check", 2)
        box: dict[str, float] | None = config.get("box")
        if box is not None:
            box_x: float = config["box"].get("x", 5.0)
            box_y: float = config["box"].get("y", 5.0)
            box_z: float = config["box"].get("z", 5.0)
        else:
            box_x, box_y, box_z = 5.0, 5.0, 5.0
        random_seed: int = config.get("seed", np.random.randint(0, 1000000000))
        pos_method: str = config.get("position_method", "fcc")
        vel_method: str = config.get("velocity_method", "mbdist")
        simulator_type: list[str] = config.get("simulator_type", ["verlet"])
        outputs: list[str] = config.get(
            "outputs",
            [
                "energies",
                "distances",
                "animation",
                "pair_correlation",
                "MSD",
                "compressibility",
                "specific_heat",
            ],
        )
        enable_cache: bool = config.get("do_caching", True)
        lat_const: float | None = config.get("lattice_const")
        corner_offset: list[float] | None = config.get("corner_offset")
        if lat_const is None or corner_offset is None:
            lat_consts, corner_offsets = fcc_for_box(
                np.array([box_x, box_y, box_z]), amount_of_particles
            )
            corner_offset = list(corner_offsets)
            lat_const = float(np.mean(lat_consts))
        bin_size: float = config.get("bin_size", 0.1)
        export_csv: bool = config.get("export_csv", False)
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
            export_csv,
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
