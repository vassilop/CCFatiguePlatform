#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import os
from numpy.core.numeric import NaN
import pandas as pd
import numpy as np
import bokeh as bk
from collections import Counter
from scipy import stats

### Importing data from csv file
def read_std(data_directory, data_type, res, date, test_type, test_number, lab, researcher, loading):
        '''
        Arguments:
            hyst_df: DataFrame containing the following computed values: n_cycles, hysteresis_area, stiffness, creep_strain
            data_directory, data_type, res, date, test_type, test_number, lab, researcher, loading: Information relative to file structure

        Returns: pandas dataframe with values in standard format

        Description: This function reads the csv file created with the Standardizing_data.py method
        '''

    filename = data_type+'_'+res+'_'+date+'_'+test_type+'_'+test_number+'.txt'
    filepath = os.path.join(data_directory, lab, researcher, loading, date, data_type, filename)
    return pd.read_csv(filepath, header = 0)


### Importing meta data from csv
def read_met(data_directory, data_type, res, date, test_type, test_number, lab, researcher, loading):
        '''
        Arguments:
            hyst_df: DataFrame containing the following computed values: n_cycles, hysteresis_area, stiffness, creep_strain
            data_directory, data_type, res, date, test_type, test_number, lab, researcher, loading: Information relative to file structure

        Returns: pandas dataframe with metadata

        Description: This function reads the JSON file created with the Standardizing_data.py method
        '''

    filename = data_type+'_'+res+'_'+date+'_'+test_type+'_'+test_number+'.txt'
    filepath = os.path.join(data_directory, lab, researcher, loading, date, data_type, filename)
    return pd.read_csv(filepath, header = 0)

### Writing output file (HYS)
def write_hys(hyst_df, data_directory, data_type, res, date, test_type, test_number, lab, researcher, loading):
    '''
    Arguments:
        hyst_df: DataFrame containing the following computed values: n_cycles, hysteresis_area, stiffness, creep_strain
        data_directory, data_type, res, date, test_type, test_number, lab, researcher, loading: Information relative to file structure

    Returns: Void

    Description: This function takes all the data computed within Hysteresis_loops.py and writes them as CSV to a specified location
    '''
    filename = data_type+'_'+res+'_'+date+'_'+test_type+'_'+test_number+'.txt'
    filepath = os.path.join(data_directory, lab, researcher, loading, date, data_type, filename)
    hyst_df.to_csv(path_or_buf=filepath, index=False)
    return



### Create dataframe function
def create_hyst_df():
    '''
    Arguments: No Arguments

    Returns: Initializes a dataframe with 4 empty ColumnDataSource

    Description: In this function we use the columns described by hyst_col and initialize an empty list that will contain all the computed information

    '''
    hyst_col = ['n_cycles', 'hysteresis_area', 'stiffness', 'creep']
    return pd.DataFrame(columns = hyst_col)


### Isolate Stress/Strain values for each cycles
def get_stress_strain(df, n_cycles):
    '''
    Arguments:
        df: df is the dataframe generated by the Standardizing_data.py method with information about Stress and Strain arranged in standard format
        n_cycles:

    Returns: Void

    Description:
        This function takes the Machine_Load and Machine_Displacement columns from the standard format dataframe and isolates all the stress strain values for each cycles
        into a list of lists
    '''
    Stress_N = []
    Strain_N = []
    df0=np.array(df)
    for k in range(len(n_cycles)):
        Stress_N.append([df0[i,1] for i in range(df0.shape[0]) if df0[i,0]==n_cycles[k]])
        Strain_N.append([df0[i,2] for i in range(df0.shape[0]) if df0[i,0]==n_cycles[k]])

    return Stress_N, Strain_N

def fill_hyst_df(df, meta_df, hyst_df):

    '''
    Arguments:
        df: dataframe with raw data in standard format
        meta_df: metadata information contained in JSON file
        hyst_df: initialized dataframe with four columns generated with create_hyst_df()

    Returns: Void

    Description: This function does all the computations on the hysteresis loops allowing to find values for hysteresis area, stiffness, creep and the number of cycles
        In the first part we take the column Machine_N_cycles from the standard dataframe and make it more compact by allowing one row per cycles
        We then initialize 3 empty lists for the computations of hysteresis area, creep and stiffness
        The function get_stress_strain is then called, creating a list of len(n_cycles) lists each containing stress/strain information for a single hysteresis loop
        Then using a shoelace algorithm, we compute the area of every polygon described by hysteresis loops and store them in the Hysteresis_Area list
        Using the same points as for area computations, we evaluate the value for stiffness by performing a linear regression on each loop and extracting the slope
        We then evaluate the value for creep by taking the arithmetic mean the maximum and minimum value for strain at each cycle

    '''

    # NB CYCLES
    ### Extract number of cycles without repeating values in other table, and store number of measurements per cycle
    n_cycles = sorted(set(df.Machine_N_cycles))
    hyst_df.n_cycles = n_cycles
    #n_measurements = Counter(df.Machine_N_cycles)
    #n_cycles_df = pd.DataFrame(n_cycles)


    # HYSTERESIS & STIFFNESS
    ### Definition of polyarea function
    def PolyArea(x,y):
        return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

    Hysteresis_Area = []
    Stiffness = []
    creep = []
    n = 0
    Stress_N, Strain_N = get_stress_strain(df, n_cycles)
    for j in range(len(n_cycles)):
        x=Stress_N[j]
        y=Strain_N[j]
        points=(x, y)
        if j<(len(n_cycles)-1):
            Hysteresis_Area.append(PolyArea(x, y)*(n_cycles[j+1]-n_cycles[j]))
        else: Hysteresis_Area.append(PolyArea(x, y)*(int(meta_df.N_fail[0])-n_cycles[j]))
        if j > 0:
            slope, _, _, _, _ = stats.linregress(y, x)
            Stiffness.append(slope)
    hyst_df.hysteresis_area = Hysteresis_Area
    Stiffness.insert(0, np.nan)

    # CREEP
    ### Creep computation
    #Stress_max = [np.max(Stress_N[i]) for i in range(1,len(n_cycles))]
    #Stress_min = [np.min(Stress_N[i]) for i in range(1,len(n_cycles))]
    Strain_max = [np.max(Strain_N[i]) for i in range(1,len(n_cycles))]
    Strain_min = [np.min(Strain_N[i]) for i in range(1,len(n_cycles))]

    for i in range(len(n_cycles)-1):
        #dStress = Stress_max[i]-Stress_min[i]
        #dStrain = Strain_max[i]-Strain_min[i]
        #Stiffness.append(dStress/dStrain)
        creep.append((Strain_max[i]+Strain_min[i])/2)
    creep.insert(0, np.nan)
    hyst_df.stiffness = Stiffness
    hyst_df.creep = creep

    return


def main():

    res = 'VAH'
    date = '210420'
    test_type = 'FA'
    test_number = '002'
    lab = "CCLab"
    researcher = 'Vahid'
    loading = 'Fatigue'
    data_directory = '/Volumes/GoogleDrive/.shortcut-targets-by-id/306/FatigueDataPlatform files & data/Data Description/File directory example'

    df = read_std(data_directory, 'STD', res, date, test_type, test_number, lab, researcher, loading)
    meta = read_met(data_directory, 'MET',  res, date, test_type, test_number, lab, researcher, loading)
    hyst = create_hyst_df()
    fill_hyst_df(df, meta, hyst)

    write_hys(hyst, data_directory, 'HYS', res, date, test_type, test_number, lab, researcher, loading)

    return


if __name__ == "__main__":
    main()
