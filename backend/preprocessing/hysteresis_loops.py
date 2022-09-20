#!/usr/bin/env python

"""
Directly copied
 from 30dd3efd7d1565fbd80d6da1744cb1e4a601c909
 Fatigue_test_dashboard/Hysteresis_loops.py (written by Scott)
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
from tst_data_lib import (
    Experiment,
    Logger,
    EXPERIMENTS_FOLDER,
    RAW_EXPERIMENT_FP_FOLDERS,
)


# Constantes
# DATA_DIRECTORY = "/Volumes/GoogleDrive/.shortcut-targets-by-id/306/FatigueDataPlatform files & data/Data Description/File directory example/"
# DATE = "2021-04-20"
# TEST_TYPE = "FA"
# TEST_NUMBER = "002"
# LAB = "CCLAB"
# RESEARCHER = "Vahid"
# DATA_STEP_IN = "STD"
# DATA_STEP_OUT = "HYS"


### Isolate Stress/Strain values for each cycles
def get_stress_strain(df, n_cycles):
    """
    Arguments:
        df: df is the dataframe generated by the Standardizing_data.py method with information about Stress and Strain arranged in standard format
        n_cycles:

    Returns: Void

    Description:
        This function takes the Machine_Load and Machine_Displacement columns from the standard format dataframe and isolates all the stress strain values for each cycles
        into a list of lists
    """
    Stress_N = []
    Strain_N = []
    for k in range(len(n_cycles)):
        Stress_N.append(list(df[df.Machine_N_cycles == n_cycles[k]].Machine_Load))
        Strain_N.append(
            list(df[df.Machine_N_cycles == n_cycles[k]].Machine_Displacement)
        )

    return Stress_N, Strain_N


def fill_hyst_df(df, meta_df, hyst_df):

    """
    Arguments:
        df: dataframe with raw data in standard format
        meta_df: metadata information contained in JSON file
        hyst_df: initialized dataframe with four columns generated with create_hyst_df()

    Returns: Void

    Description:
        This function does all the computations on the hysteresis loops allowing to find values for hysteresis area, stiffness, creep and the number of cycles
        In the first part we take the column Machine_N_cycles from the standard dataframe and make it more compact by allowing one row per cycles
        We then initialize 3 empty lists for the computations of hysteresis area, creep and stiffness
        The function get_stress_strain is then called, creating a list of len(n_cycles) lists each containing stress/strain information for a single hysteresis loop
        Then using a shoelace algorithm, we compute the area of every polygon described by hysteresis loops and store them in the Hysteresis_Area list
        Using the same points as for area computations, we evaluate the value for stiffness by performing a linear regression on each loop and extracting the slope
        We then evaluate the value for creep by taking the arithmetic mean the maximum and minimum value for strain at each cycle

    """

    # NB CYCLES
    ### Extract number of cycles without repeating values in other table, and store number of measurements per cycle
    n_cycles = sorted(set(df.Machine_N_cycles))
    hyst_df.n_cycles = n_cycles
    # n_measurements = Counter(df.Machine_N_cycles)
    # n_cycles_df = pd.DataFrame(n_cycles)

    # HYSTERESIS & STIFFNESS
    ### Definition of polyarea function
    def PolyArea(x, y):
        """
        Arguments:
            x: values measured along the x axis (strain)
            y: values measured along the y axis (stress)

        Returns:
            Single value for the area of the loop

        Description:
            The PolyArea function computes the area of each hysteresis loops using a shoelace algorithm
        """
        return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

    Hysteresis_Area = []
    Stiffness = []
    creep = []
    Stress_N, Strain_N = get_stress_strain(df, n_cycles)
    for j in range(len(n_cycles)):
        x = Stress_N[j]
        y = Strain_N[j]
        if j < (len(n_cycles) - 1):
            Hysteresis_Area.append(PolyArea(x, y) * (n_cycles[j + 1] - n_cycles[j]))
        else:
            Hysteresis_Area.append(PolyArea(x, y) * (int(meta_df.N_fail) - n_cycles[j]))
        if j > 0:
            slope, _, _, _, _ = stats.linregress(y, x)
            Stiffness.append(slope)
    hyst_df.hysteresis_area = Hysteresis_Area
    Stiffness.insert(0, np.nan)

    # CREEP
    ### Creep computation
    # Stress_max = [np.max(Stress_N[i]) for i in range(1,len(n_cycles))]
    # Stress_min = [np.min(Stress_N[i]) for i in range(1,len(n_cycles))]
    Strain_max = [np.max(Strain_N[i]) for i in range(1, len(n_cycles))]
    Strain_min = [np.min(Strain_N[i]) for i in range(1, len(n_cycles))]

    for i in range(len(n_cycles) - 1):
        # dStress = Stress_max[i]-Stress_min[i]
        # dStrain = Strain_max[i]-Strain_min[i]
        # Stiffness.append(dStress/dStrain)
        creep.append((Strain_max[i] + Strain_min[i]) / 2)
    creep.insert(0, np.nan)
    hyst_df.stiffness = Stiffness
    hyst_df.creep = creep


def main():
    with Logger(None) as logger:
        for experiment_raw_fp_folder in RAW_EXPERIMENT_FP_FOLDERS:
            experiment_metadata = Experiment(experiment_raw_fp_folder, logger)
            print(f"parsing {experiment_raw_fp_folder}")
            for measures in experiment_metadata.exp_meta_meta["measures"]:
                try:
                    logger.info(
                        f"Read measures {os.path.basename(measures['preprocessed_fp'])}"
                    )
                    with logger.indent:
                        df = pd.read_csv(measures["preprocessed_fp"])

                    hyst_df = pd.DataFrame(
                        columns=["n_cycles", "hysteresis_area", "stiffness", "creep"]
                    )

                    experiment_metadata.N_fail = 100
                    fill_hyst_df(df, experiment_metadata, hyst_df)
                    hys_fp = "HYS_" + os.path.basename(measures["preprocessed_fp"])
                    hyst_df.to_csv(
                        os.path.join(
                            experiment_metadata.exp_meta_meta["preprocessed_folder"],
                            hys_fp,
                        ),
                        index=False,
                    )
                except:
                    pass  # missing data to produce HYS -> skip for now


if __name__ == "__main__":
    main()
