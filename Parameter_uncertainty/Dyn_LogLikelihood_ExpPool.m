% Hromada et al "Negative interactions determine C. difficile growth in
% synthetic human gut communities". Submitted for publication in Molecular
% Systems Biology. Created by Yili Qian, Venturelli Lab, Jan 2021.

function LogLikelihood = Dyn_LogLikelihood_ExpPool(exp_species,data_files,para_current,dyn_sys,num_states)
% compute loglikelihood for a set of experiments

% number of experimental data sets to pool together
L = length(exp_species);

% put parameters in matrix format
r = para_current(1:num_states);
A = reshape(para_current(num_states+1:end),num_states,num_states)';

% cumulative likelihood
LogLikelihood_cum = 0;

% this loads the data from all experimental conditions
[T,Y_mean,Y_std] = LoadFileFcn(data_files);

for i = 1:1:L
    species_idx = exp_species{i};
    
    % parameters needed for this experiment
    r_exp = r(species_idx);
    A_exp = A(species_idx,species_idx);
    
    [nn,~] = size(A_exp);
    
    A_exp_vec = reshape(A_exp',1,nn^2);
    
    para_exp_vec = [r_exp,A_exp_vec];
    
    y_mean = Y_mean{i};
    y_std = Y_std{i};
    y_std = max(y_std,0.025);
    
    % use the first measurement as initial condition
    X0 = y_mean(1,:);
    
    LogLikelihood_cum = LogLikelihood_cum + Dyn_LogLikelihood_SingleExp(para_exp_vec,dyn_sys,T{i},X0,y_mean,y_std);
end
LogLikelihood = LogLikelihood_cum;