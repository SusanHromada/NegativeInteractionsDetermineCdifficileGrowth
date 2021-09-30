# Overview
These scripts were used to estimate parameters of the generalized Lotka-Volterra model based on sets of experimental training data. They were used to generate the Preliminary Model and Full Model in Hromada et al. 2021.

## Authorship
The scripts were written by Dr. Ryan Clark and are similar to those published in Clark, R.L., Connors, B.M., Stevenson, D.M. et al. Design of synthetic human gut microbiome assembly and butyrate production. Nat Commun 12, 3254 (2021).

## Information 
All MATLAB scripts were written for MATLAB R2018b and run on a private computational server (Linux).

## Input data
The parameter estimation scripts read in data from a single folder containing a set of .csv files, with one .csv file for each initial condition of each community. These folders can be found at Parameter_estimation/input_data. The format of the .csv files is explained in Parameter_estimation/input_data/Experiment_csv_file_format.png
SEH22_PreliminaryModel is used to infer the Preliminary Model
SEH25_FullModel is used to infer the FullModel
SEH24_FullModel_validation_set_removed leaves out 24 communities that can be used as validation test data. 
