import os
import numpy as np
import pandas as pd
from typing import List, Dict
from pydantic import BaseModel
from pandas.core.frame import DataFrame
from sqlalchemy.future import select
from sqlalchemy.ext.asyncio import AsyncSession
from ccfatigue.models.database import Experiment, Test
from ccfatigue.config import settings

DATA_DIRECTORY: str = os.path.join(settings.data_path, "preprocessed")  # type: ignore
INTERVAL: int = 10
LOOP_SPACING: int = 1000
MAGNITUDE: int = -3


class HysteresisLoop(BaseModel):
    n_cycles: List[float]
    strain: List[float]
    stress: List[float]


class TestResult(BaseModel):
    specimen_id: int
    total_dissipated_energy: int
    hysteresis_loops: List[HysteresisLoop]
    n_cycles: List[float]
    creep: List[float]
    hysteresis_area: List[float]
    stiffness: List[float]


def get_dataframe(
    data_in: str,
    exp: Dict[str, str],
    specimen_id: int,
) -> DataFrame:
    """
    return extracted DataFrame related to that test from CSV
    """
    # FIXME researcher_name from a column value
    researcher_name = exp["researcher"].split(" ")[-1]
    if data_in == "HYS":
        filepath = os.path.join(
            DATA_DIRECTORY,
            f"TST_{researcher_name}_{exp['date']}_{exp['experiment_type']}",
            f"{data_in}_measure_{specimen_id:03d}.csv",
        )
    else:
        filepath = os.path.join(
            DATA_DIRECTORY,
            f"{data_in}_{researcher_name}_{exp['date']}_{exp['experiment_type']}",
            f"measure_{specimen_id:03d}.csv",
        )
    abspath = os.path.abspath(filepath)
    print(abspath)
    return pd.read_csv(abspath)


def get_total_dissipated_energy(hyst_df: DataFrame) -> int:
    """
    return calculated Total Dissipated Energy (TDE)
    """
    return np.sum(hyst_df["hysteresis_area"])


def compute_sub_indexes(df: DataFrame) -> List[int]:
    """
    Compute subset of index used for plotting curves
    """
    n_cycles_max = np.max(df.n_cycles)

    sub_index_intermediate = np.geomspace(
        start=LOOP_SPACING, stop=n_cycles_max, num=INTERVAL
    )
    sub_index = np.round(sub_index_intermediate, MAGNITUDE)
    sub_index[0] = 1

    return sub_index.tolist()


def create_sub_hystloops(df: DataFrame, sub_indexes: List[int]) -> List[HysteresisLoop]:
    """
    Conditional plotting of hysteresis loops
    we only plot the loops specified by sub_index
    """
    sub_hystloops: List[HysteresisLoop] = []
    for sub_index in sub_indexes:
        sub_hystloops_strain = []
        sub_hystloops_stress = []
        sub_hystloops_ncycles = []
        for i in range(len(df)):
            if df.Machine_N_cycles[i] == sub_index:
                sub_hystloops_strain.append(df.Machine_Displacement[i])
                sub_hystloops_stress.append(df.Machine_Load[i])
                sub_hystloops_ncycles.append(df.Machine_N_cycles[i])
        # make curve closed // Optional
        if len(sub_hystloops_ncycles) != 0:
            sub_hystloops_strain.append(sub_hystloops_strain[0])
            sub_hystloops_stress.append(sub_hystloops_stress[0])
            sub_hystloops_ncycles.append(sub_hystloops_ncycles[0])
        if (
            len(sub_hystloops_ncycles) > 0
            and len(sub_hystloops_strain) > 0
            and len(sub_hystloops_stress) > 0
        ):
            sub_hystloops.append(
                HysteresisLoop(
                    n_cycles=sub_hystloops_ncycles,
                    strain=sub_hystloops_strain,
                    stress=sub_hystloops_stress,
                )
            )
    return sub_hystloops


async def get_result(
    session: AsyncSession,
    experiment_id: int,
    test_id: int,
) -> TestResult:
    experiment: Dict[str, str] = (
        (
            await session.execute(
                select(
                    Experiment.laboratory,
                    Experiment.researcher,
                    Experiment.experiment_type,
                    Experiment.date,
                ).where(Experiment.id == experiment_id)
            )
        )
        .one()  # type: ignore
        ._asdict()
    )
    specimen_id: int = await session.scalar(
        select(Test.specimen_number)
        .where(Test.experiment_id == experiment_id)
        .where(Test.id == test_id)
    )  # type: ignore
    std_df = get_dataframe("TST", experiment, specimen_id)
    hyst_df = get_dataframe("HYS", experiment, specimen_id).fillna(value=0)
    hysteresis_loops = create_sub_hystloops(std_df, compute_sub_indexes(hyst_df))
    return TestResult(
        specimen_id=specimen_id,
        total_dissipated_energy=get_total_dissipated_energy(hyst_df),
        hysteresis_loops=hysteresis_loops,
        n_cycles=hyst_df["n_cycles"].to_list(),
        creep=hyst_df["creep"].to_list(),
        hysteresis_area=hyst_df["hysteresis_area"].to_list(),
        stiffness=hyst_df["stiffness"].to_list(),
    )
