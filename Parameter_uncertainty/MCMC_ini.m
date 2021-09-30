% Hromada et al "Negative interactions determine C. difficile growth in
% synthetic human gut communities". Submitted for publication in Molecular
% Systems Biology. Created by Yili Qian, Venturelli Lab, Jan 2021.

clear all; close all; clc; warning off
%% MCMC setup

dyn_sys = @gLV_ParaVec;

num_chains = 1;

num_states = 14; % total number of species considered (i.e., set from which subcommunities are picked)

num_iti = 50000;   % number of samples generated

k_adapt = 0.5;  % adaptive stepping size. std of markov jump = k_adapt*std of previous jumps

burn_in = 100;

display(['num_iti=' num2str(num_iti)])
%% initialization parameters
% initial step size for the burn-in period. Use adaptive stepping afterwards
delta_r = 0.0001*ones(1,num_states);
delta_A = 0.0001*ones(num_states,num_states);

% minimum step size
step_min = 0.0001;

step_ini = [delta_r reshape(delta_A',1,num_states^2)];

% initial parameter guess
load('para_seed')

%% set up experimental data
data_file = 'ExpData';    % experimental data extracted from the same file
load(['data/' data_file])

% measured species in each experiment
exp_species = species_idx_vec;

%% run MCMC 
for q = 1:1:num_chains
    [para_vec,acceptance_vec,LogPosterior_vec] = AdaptiveMCMC_DynPara_v2(para0,dyn_sys,exp_species,...
        data_file,num_iti,burn_in,step_ini,k_adapt,step_min,num_states);
end