% Hromada et al "Negative interactions determine C. difficile growth in
% synthetic human gut communities". Submitted for publication in Molecular
% Systems Biology. Created by Yili Qian, Venturelli Lab, Jan 2021.

function [para_vec,acceptance_vec,LogPosterior_vec] = AdaptiveMCMC_DynPara_v2(para0,dyn_sys,exp_species,...
    data_files,num_iti,burn_in,step_ini,k_adapt,step_min,num_states)

% Adaptive random walkl MCMC to identify dynamical system parameters.
%%%%%%%%%%%%%
% Inputs:
% para0 - initial guess of the parameters [row vector 1 x p, p = # of parameters]
% dyn_sys - dynamical system to simulate [use function handle to specify]
% data_file - file name of experimental data [string, saved in subfolder named '/data']
% num_iti - total number of iterations [positive integer]
% burn_in - # of iterations to discard as burn-in [poisitve integer << num_iti]
% step_ini - initial std of proposal (normal) distribution before
% adaptation begines [row vector 1xp]
% k_adapt - adaptation factor that determines step size online [scalar]
% exp_species - a 1xM cell with each cell indicating the species measured
% in that experiment, M is the total number of experiments
%%%%%%%%%%%%
% Outputs:
% para_vec - (parameter) samples drawn from posterior [matrix num_iti x p]
% acceptance_vec - whether a (parameter) sample is accepted (1) or rejected (0) [vector num_itix1]
% LogPosterior_vec - posterior probability density of each sample

%% set up
% get the number of parameters
num_parameters = length(para0);

% records log likelihood (start with index 0 to always accept the first sample)
LogPosterior_vec = -inf*ones(num_iti+1,1);

% records if a sample is accepted or not
acceptance_vec = zeros(num_iti,1);

% acceptance probability
a_vec = zeros(num_iti+1,1);

% sampled r vector
para_vec = zeros(num_iti+1,num_parameters);

% generate uniform sample from [0,1] to compare with acceptance ratio take
% log because the likelihoods are also in log
seed = log((rand(num_iti,1)));

% generate normal sample from [0,1] to perform Markovian jump
para_seed_jump = randn(num_iti,num_parameters);

%% compute loglikelihood of para0
% loglikelihood associated with dynamic measurements
LogLikelihood_dyn = Dyn_LogLikelihood_ExpPool(exp_species,data_files,para0,dyn_sys,num_states);

% compute log posterior
load('det_fit') % use diterministic fit to set the priors
para_det = [r_vec A_vec];
% the second the third arguments in LogPrior are deterministic fits
LogPosterior = LogLikelihood_dyn + LogPrior_v2(para0,para_det);

LogPosterior_vec(1) = LogPosterior;

%% run Markov chain
para_vec(1,:) = para0;

% current state of parameters
para_current = para0;

% set jump std during burn-in to be the initial step size
delta_para = step_ini;

for j = 2:1:num_iti+1
    
    % display progress every 100 jumps
%     if mod(j,10) == 0
        display(['iteration ' num2str(j)])
%     end
    
    % adaptive MCMC scheme from Robserts and Rosenthal 2008
    if j>burn_in
%         para_std = 2.38^2*std(para_vec(1:j-1,:))/num_parameters;
        para_std = k_adapt*std(para_vec(1:j-1,:))/num_parameters;
        
        % use standard deviation from existing samples to set new "jump size"
        delta_para = max(para_std,step_min);
    end
    
    % propose next jump state
    para_current = para_current + para_seed_jump(j-1,:).*delta_para;
    
    % compute log likelihood
    LogLikelihood_dyn = Dyn_LogLikelihood_ExpPool(exp_species,data_files,para_current,dyn_sys,num_states);
    
    % compute log posterior
    LogPosterior = LogLikelihood_dyn + LogPrior_v2(para_current,para_det);
    
    % log acceptance ratio
    a = min(0,LogPosterior-LogPosterior_vec(j-1));
    a_vec(j-1) = a;
    
    if a > seed(j-1)
        % if greater than seed, accept, record parameter sample and likelihood
        acceptance_vec(j-1) = 1;
        
        para_vec(j,:) = para_current;

        LogPosterior_vec(j) = LogPosterior;
    else
        % if less than seed, reject, do not record parameter, replace
        % likelihood with previous likelihood
        acceptance_vec(j-1) = 0;
        
        para_vec(j,:) = para_vec(j-1,:);
        
        para_current = para_vec(j-1,:);
        
        LogPosterior_vec(j) = LogPosterior_vec(j-1);
    end
end