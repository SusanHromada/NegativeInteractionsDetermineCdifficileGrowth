% Hromada et al "Negative interactions determine C. difficile growth in
% synthetic human gut communities". Submitted for publication in Molecular
% Systems Biology. Created by Yili Qian, Venturelli Lab, Jan 2021.

function dx = gLV_ParaVec(t,x,para)

% gLV taking vector parameters

% find the number of species
N = (-1+sqrt(1+4*length(para)))*0.5;

% first N values are growth rates
R_model = para(1:N);

% rest are interaction coefficients -> matrix form
A_model_vec = para(N+1:end);
A_model = reshape(A_model_vec',N,N)';

dx = zeros(N,1);

for i = 1:1:N
    
    Interactions = 0;
    for j = 1:1:N
        Interactions = Interactions + A_model(i,j)*x(j);
    end
    
    dx(i) = x(i)*(R_model(i)+Interactions);
end