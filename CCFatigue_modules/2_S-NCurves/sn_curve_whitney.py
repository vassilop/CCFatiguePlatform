#!/usr/bin/env python
""" CCFatigue - Module 2 - sn_curve_whitney.py
This code takes in AGG data and output SNC data
"""

import os
import math
import numpy as np
import pandas as pd
from scipy import stats
from itertools import chain


SRC_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(SRC_DIR, "..", "..", "Data")


INPUT_AGG_FILENAME = "AGG_input.csv"
INPUT_AGG_FILE = os.path.join(DATA_DIR, INPUT_AGG_FILENAME)
OUTPUT_SNC_JSON_FILENAME = "SNC_LogLog.json"
OUTPUT_SNC_JSON_FILE = os.path.join(DATA_DIR, OUTPUT_SNC_JSON_FILENAME)
OUTPUT_SNC_CSV_FILENAME = "SNC_LogLog.csv"
OUTPUT_SNC_CSV_FILE = os.path.join(DATA_DIR, OUTPUT_SNC_CSV_FILENAME)

CYCLES_TO_RUN_OUT = 5e7


def calculate_stress_ratio(
    so: float, pn: float, pow: float, faf: float, nn: float, f: float
) -> float:
    """
    Calculate stress ratio (SPN)
    https://github.com/EPFL-ENAC/CCFatiguePlatform/blob/develop/CCFatigue_modules/2_S-NCurves/S-N-Curve-Whitney.for#L329
    Inputs:
    - so: float
    - pn: float
    - pow: float
    - faf: float
    - nn: float
    - f: float
    Output:
    - stress ratio
    """
    spn = so * ((-math.log(pn(f) / 100)) ** (pow / faf)) * (nn ** (-pow))
    return spn


if __name__ == "__main__":

    # Import input file (AGG format)
    agg_df = pd.read_csv(INPUT_AGG_FILE)

    # Data are grouped by stress_ratio but one experiment
    # can have two separate groups with same stress_ratio so we need to identify
    agg_df["stress_ratio_id"] = (
        agg_df.stress_ratio != agg_df.stress_ratio.shift()
    ).cumsum()

    # Identify groups of stress level
    agg_df["stress_level_id"] = (
        agg_df.stress_level != agg_df.stress_level.shift()
    ).cumsum()

    # Average by stress_ratio
    avg_stress_ratios_df = (
        agg_df[
            [
                "stress_ratio_id",
                "stress_ratio",
                "cycles_to_failure",
            ]
        ]
        .groupby(["stress_ratio_id"])
        .mean()
    )

    avg_stress = (
        agg_df[
            [
                "stress_ratio_id",
                "stress_level_id",
                "stress_parameter",
                "cycles_to_failure",
            ]
        ]
        .groupby(["stress_ratio_id", "stress_level_id"])
        .mean()
    )

    # calculate average number of cycles at each stress level
    # ##########################################################################
    avg_stress["rout"] = CYCLES_TO_RUN_OUT / avg_stress.cycles_to_failure

    # # Number of failed data at each stress level
    # failed_data = (
    #     agg_df[["stress_ratio_id", "stress_level_id", "cycles_to_failure"]]
    #     .groupby(["stress_ratio_id", "stress_level_id"])
    #     .apply(lambda x: (x.cycles_to_failure < CYCLES_TO_RUN_OUT).sum())
    #     .reset_index(name="count")
    # )

    cycles_to_failure = pd.DataFrame(
        chain(
            range(1, 1000, 50),
            range(1000, 1001),
            range(10000, 2000000, 10000),
            range(3000000, 20000000, 1000000),
            range(30000000, 1400000000, 100000000),
        ),
        columns=["cycles_to_failure"],
    )

    snc_output_csv_df = pd.DataFrame(
        columns=[
            "stress_ratio",
            "cycles_to_failure",
            "stress_parameter",
            "stress_lowerbound",
            "stress_upperbound",
        ]
    )

    # Calculate A (Eq 4) and B (Eq 5)
    # linregress = agg_df.groupby("stress_ratio_id").apply(
    #     lambda x: stats.linregress(x.log10_stress_parameter, x.log10_number_of_cycles)
    # )

    # # Eq 4 and 5
    # # Add A and B in avg_by_stress_ratio
    # # Slope = A, Intercept = B
    # stress_ratios_df["slope"] = linregress.apply(lambda x: x[1])
    # stress_ratios_df["intercept"] = linregress.apply(lambda x: x[0])

    # Prepare SNC output list of cycles to failure

    snc_output_json_df = pd.DataFrame(
        columns=[
            "stress_ratio",
            "RSQL",
            "A",
            "B",
            "Fp",
        ],
    )

    snc_output_csv_df = pd.DataFrame(
        columns=[
            "stress_ratio",
            "cycles_to_failure",
            "stress_parameter",
            "stress_lowerbound",
            "stress_upperbound",
        ]
    )

    # For each group of stress_ratio
    for (stress_ratio_id, stress_ratio_df) in avg_stress_ratios_df.iterrows():

        # Prepare CSV data for each cycles to failure
        output = cycles_to_failure.copy()

        # Perform MLE procedure for each stress level
        # Calculate the value of shape parameter at each stress level
        ##########################################################################

        # Perform the second step of MLE procedure
        # Estimate characteristic number of cycles for each stress level
        ##########################################################################

        # Normalization of the No. of cycles at each level by the obtained char. No of cycles
        ##########################################################################

        # Perform MLE procedure
        # Calculate the value of shape parameter
        ##########################################################################

        # Perform the second step of MLE procedure
        ##########################################################################

        # Calculating the slope and the intercept of the S/N curve based on Linear Regression data fitting
        ##########################################################################

        a1 = 0
        a2 = 0
        a3 = 0
        a4 = 0
        a5 = 0
        for i in range(Nostrlev):

            noi = math.log10(xo * rn(i))
            log_s = math.log10(s)
            a1 = a1 + noi * log_s
            a2 = a2 + log_s
            a3 = a3 + noi
            a4 = a4 + log_s**2
            a5 = (a2) ** 2

        slope = (Nostrlev * a1 - a2 * a3) / (Nostrlev * a4 - a5)
        intercept = (a3 - slope * a2) / Nostrlev

        # Calculate So = 10**(-intercept / slope)

        output["stress_parameter"] = output.apply(
            lambda x: calculate_stress_ratio(
                nn=x,
            )
        )

        # stress_bounds = stress_parameter.apply(
        #     lambda x: stress_at_failure_bounds(
        #         stress_ratio_df.sample_count,
        #         stress_ratio_df.q,
        #         stress_ratio_df.slope,
        #         stress_ratio_df.intercept,
        #         x.cycles_to_failure,
        #         stress_ratio_df.pp,
        #         stress_ratio_df.xb,
        #     ),
        #     axis=1,
        # )
        # stress_parameter["stress_lowerbound"] = stress_bounds.apply(lambda x: x[0])
        # stress_parameter["stress_upperbound"] = stress_bounds.apply(lambda x: x[1])

        # stress_parameter["stress_ratio"] = stress_ratio_df.stress_ratio

        # snc_output_csv_df = pd.concat([snc_output_csv_df, stress_parameter])

        # # Prepare JSON
        # json_df = pd.DataFrame(
        #     {
        #         "stress_ratio": stress_ratio_df.stress_ratio,
        #         "RSQL": RELIABILITY_LEVEL,
        #         "A": stress_ratio_df.stress_parameter_aa,
        #         "B": stress_ratio_df.stress_parameter_bb,
        #         "Fp": stress_ratio_df.fp,
        #     },
        #     index=[0],
        # )
        # snc_output_json_df = pd.concat([snc_output_json_df, json_df], ignore_index=True)

    # Export dataframes to files
    snc_output_json_df.to_json(OUTPUT_SNC_JSON_FILE, orient="records")

    snc_output_csv_df.to_csv(
        OUTPUT_SNC_CSV_FILE,
        index=False,
    )
