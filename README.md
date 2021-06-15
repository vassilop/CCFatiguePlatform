# CCFatigue Platform

CCFatiguePlatform is an initiative from CCLab that aims to develop a web application to faciltate manipulation and harmonized storage of composite materials testing datasets. 

# Data
## Standard data format
[Example of data in standard format](https://drive.google.com/file/d/1-SuUHPbW-xFr65yqVIbrqHl1vb4ejMzo/view?usp=sharing "Shayan's data in standard format")




![Metadata Scheme](/Volumes/GoogleDrive/.shortcut-targets-by-id/306/FatigueDataPlatform files & data/Data Description/Metadata/Metadata Scheme - visual.png)



|Directory naming convention|
|---------------------------|
|\Laboratory\Researcher\Test Type\Date\ *filename
|Laboratory: Full accronym (i.e. CCLab, RESSLab, IBeton,...)
|Researcher: Last name
|Test type: Standard Fatigue, Combined Fatigue Fracture, Standard Quasi-Static, Combined Quasi-Static Fracture
|Date: YYYY-MM-DD


|Filenaming conventions|
|-----------------------|
|{Researcher} _ {Date} _ {Test type} _ {status code} _ {###}|
|Researcher: 3 first letters of last name {i.e. VAS, WEI, MAT,...}
|Date: YYMMDD
|Test type: {FA} = standard fatigue, {QS} = standard quasi-static, {FF} = combined fatigue/fracture, {SF} = combined quasi-static/fracture
|Status code: {RAW} = raw data, {TRE} = treated
|###: numerical identifier
|Example: SHA _ 210420 _ FF _ RAW _ 001

# Contents

### Standardizing data
Code Standardizing_data.py converts raw data from use cases to the standard format
Raw data format supported: The raw data is read in CSV format

Usage: python run Standardizing_data.py

### Hysteresis analysis
The code section 'Hysteresis_loops.py' uses the data in it's standard form and performs computations to extract:

Hysteresis Area:
>Area within each of the loading/unloading cycles

Stiffness:
>Slope of linear regression for each loading/unloading cycles

Creep:
>Residual value of strain for each loading/unloading cycle

Output is a csv file with 4 columns: Number of cycles (n_cycles), Hysteresis area (hysteresis_area), Stiffness (stiffness), Creep (creep)



### Plots
Plotting.py creates the following plots:

Stress / strain curve:

![Stress - Strain](/Volumes/GoogleDrive/.shortcut-targets-by-id/306/FatigueDataPlatform files & data/Data Description/Captures d'ecran - plots/StressStrain.png "Stress - Strain")
>On this graph we use the raw inputs for Stress and Strain, specific loops have been selected in order for the graph not to be too busy

![Stiffness Evolution](/Volumes/GoogleDrive/.shortcut-targets-by-id/306/FatigueDataPlatform files & data/Data Description/Captures d'ecran - plots/Stiffness.png "Stiffness evolution")
>This graph shows the evolution of the sample's stiffness. Stiffness corresponds to the slope of hysteresis loops and is closely linked to hooke's law. E is comparable to the constant k in the context of springs. The analysis is made from the stress - strain raw data, for each hysteresis loop, we evaluate the slope with a linear regression for each hysteresis loops for stress and strain.

![Hysteresis area evolution](/Volumes/GoogleDrive/.shortcut-targets-by-id/306/FatigueDataPlatform files & data/Data Description/Captures d'ecran - plots/hystarea.png "Hysteresis area evolution")
>Here we see how the area of each hysteresis loop evolves as the fatigue test goes on. The hysteresis loop area is closely related to the amount of energy dissipated through deformations, the sum of hysteresis areas is denoted as TDE, or total dissipated energy

![Creep evolution](/Volumes/GoogleDrive/.shortcut-targets-by-id/306/FatigueDataPlatform files & data/Data Description/Captures d'ecran - plots/creep.png "Creep evolution")
>Creep corresponds to the residual deformations during a test or life cycle of a sample

%>>>> 

### Aggregating data
aggregate.py takes all files in an experiment in standard format, analyses them and creates an aggregated set with the right format for analysis in CCFatigue. Output is a csv file with 6 columns: Stress ratio, Reliability level, Stress level no., Stress parameter, Number of cycles to failure, Residual strength 

### Modules
Using the notebooks in the modules folder, one can compute the right parameters for plotting i.e. the notebook 'Hysteresis loops.ipynb' uses stress/strain info to compute the TDE and evolution of stiffness


# Reproduce

[Environment file](environment.yml "Environment file")

# License
TBD.

# Authors
Charlotte Weil, Scott M. Salmon, Samuel Bancal.
