% Hromada et al "Negative interactions determine C. difficile growth in
% synthetic human gut communities". Submitted for publication in Molecular
% Systems Biology. Created by Yili Qian, Venturelli Lab, Jan 2021.

function LogLikelihood = Dyn_LogLikelihood_SingleExp(para_current,dyn_sys,time,X0,Y_mean,Y_std)
% compute loglikelihood for a single experiment

L = length(time);

[~,X] = ode23s(@(t,x) dyn_sys(t,x,para_current),time,X0);

% % discard states not included in species_idx (e.g., when these species are not included in experiments,
% % so mean and std are both 0).
% X = X(:,species_idx);

% if system unstable set trajectory to infity
[L_sim,W_sim] = size(X);

if L_sim < L
    X = inf*ones(L,W_sim);
end

p_y_theta = zeros(L-1,1);   % compute likelihood of each temporal data point, remove initial point since it is deterministic

for i = 1:1:L-1
    for j = 1:1:W_sim
        % evaluate likelihood: assume channels are independent
        p_y_theta(i) = p_y_theta(i)+lognormalpdf(X(i+1,j),Y_mean(i+1,j),Y_std(i+1,j));
    end
end

LogLikelihood = sum(p_y_theta);

function y = lognormalpdf(x,y_mean,y_std)

if sum(isnan(y_mean)) == 0
    y = log(normpdf(x,y_mean,y_std));
else
    y = 0;
end