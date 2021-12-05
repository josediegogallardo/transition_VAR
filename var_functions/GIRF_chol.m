function [GIRF] = GIRF_chol(var_results, IRper, shockv, GIRF_opt)
%% Reference
% Written by J. Diego Gallardo (2021) to get Generalized Impulse Response
% Funtions(*) with a Cholesky decomposition
%
% (*) Koop, G., Pesaran, H., y Potter, S. (1996). Impulse respones analysis
% in nonlinear multivariate models. Journal of Econometrics.


%% Usage
% [GIRF] = GIRF_chol(var_results, IRper, shockv, sim)
% var_results:  Results form estimatevar function
% IRper: Number of periods for 
% GIRF_opt.sim: Number of simulations
% GIRF_opt.shock: Size and sign of shock
% GIRF_opt.hist.var: variable to follow history
% GIRF_opt.hist.thresh: Threshold to switch history

% Paths to store results
path_shock = [];
path = [];

% Reorder error matrix
random_u = var_results.u(randperm(size(var_results.u, 1)), :);
% variables for bootstrap
x_boot = var_results.X;
y_boot = var_results.X*var_results.A + random_u;
A_boot = x_boot\y_boot;
u_boot = y_boot - x_boot*A_boot;
Omega_boot = cov(u_boot);
S_boot = chol(Omega_boot)';

% See Martinez & Guevara and Bigio to structural shocks by size and sign

%% Sentence with no history
% GIRF_{y}(n,v_{t}) = E[Y_{t+n}|v_{t}] - E[Y_{t+n}]
if ~isempty(GIRF_opt.hist)

    
%% Sentence with history    
% GIRF_{y}(n,v_{t},w_{t-i}) = E[Y_{t+n}|v_{t},w_{t-i}] - E[Y_{t+n}|w_{t-i}]
else 
    
end
