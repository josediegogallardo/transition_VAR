function [GIRF] = GIRF_chol(var_results, IR, shockv, GIRF_opt)
%% Reference
% Written by J. Diego Gallardo (2021) to get Generalized Impulse Response
% Funtions(*) with a Cholesky decomposition
%
% (*) Koop, G., Pesaran, H., y Potter, S. (1996). Impulse respones analysis
% in nonlinear multivariate models. Journal of Econometrics.


%% Usage
% [GIRF] = GIRF_chol(var_results, IRper, shockv, sim)
% var_results:  Results form estimatevar function
% IR.per: Number of periods for
% IR.shocktype: Type of shock (unit vs s.d.)
% GIRF_opt.sim: Number of simulations
% GIRF_opt.shock: Size and sign of shock
% GIRF_opt.hist.var: variable to follow history
% GIRF_opt.hist.thresh: Threshold to switch history

%% Initial set up
% Paths to store results
path_shock = [];
path = [];


%% Bootstrapíng
% Reorder error matrix
random_u = var_results.u(randperm(size(var_results.u, 1)), :);

% variables for bootstrap
x_boot = var_results.X;
y_boot = var_results.X*var_results.A + random_u;
A_boot = x_boot\y_boot;
u_boot = y_boot - x_boot*A_boot;
Omega_boot = cov(u_boot);
S_boot = chol(Omega_boot)';

% Shocks
switch(IR.shocktype)
    case 0
        shocks = eye(var_results.n).*(diag(S_boot).^(-1)); % one unit shock
    case 1
        shocks = eye(var_results.n); % one std deviation shock
end

% Structural errors
v = inv(S_boot)*u_boot';
v = v';
vb = v; % non shocked errors
vb_shock = v; % shocked errors
%% Sentence with no history
% GIRF_{y}(n,v_{t}) = E[Y_{t+n}|v_{t}] - E[Y_{t+n}]
if isempty(GIRF_opt.hist)
    % Random point of history
    junk=randi([1 var_results.t-IR.per]);
    % Errors in history
    vb(junk,shockv) = 0;
    vb_shock(junk,shockv) = shocks(shockv,shockv);
    % Reduced form errors
    eb   = (S_boot*vb')';
    eb_shock = (S_boot*vb_shock')';
    % Y paths
%     yb(1,:) = x_boot(junk,:)*A_boot + eb(junk,:);
%     yb_shock(1,:) = x_boot(junk,:)*A_boot + eb_shock(junk,:);
%     % Path of Y
%     for i = 1:var_results.p
%         yb(i+1,:) =  x_boot(junk+i,:)*A_boot + eb(junk+i,:);
%     end  

%% Sentence with history    
% GIRF_{y}(n,v_{t},w_{t-i}) = E[Y_{t+n}|v_{t},w_{t-i}] - E[Y_{t+n}|w_{t-i}]
else 
    disp('ohana significa familia')
end

disp('hola')
