% Hromada et al "Negative interactions determine C. difficile growth in
% synthetic human gut communities". Submitted for publication in Molecular
% Systems Biology. Created by Yili Qian, Venturelli Lab, Jan 2021.

function [T,Y_mean,Y_std] = LoadFileFcn(dataFileName)
% load experimental data file and extract mean and std

load(['data/' dataFileName])
T = time_vec;
Y_mean = data_mean_vec;
Y_std = data_std_vec;