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
%
p = var_results.p;
n = var_results.n;
c_case = var_results.c_case;
Xex = var_results.Xex;

% Paths to store results
path_shock = [];

%% Bootstrapíng
fprintf(1,'IRF collection process: \t\t\t\t');
for ii = 1:GIRF_opt.sim
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
            shocks = (eye(n).*(diag(S_boot).^(-1)))*GIRF_opt.shock; % custom shock
        case 1
            shocks = eye(n); % one std deviation shock
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
        junk=randi([1 var_results.t-IR.per-n]);
        % Errors in history
        vb(junk,shockv) = 0;
        vb_shock(junk,shockv) = shocks(shockv,shockv);
        % Reduced form errors
        eb   = (S_boot*vb')';
        eb_shock = (S_boot*vb_shock')';
        % First obs
        yb(1,:) = x_boot(junk,:)*A_boot + eb(junk,:);
        yb_shock(1,:) = x_boot(junk,:)*A_boot + eb_shock(junk,:);


        % Collecting IRF for VAR(1)
        if p == 1

            for i = 1:IR.per-1
                if isempty(Xex)
                    yb(i+1,:) = [x_boot(junk+i,1:c_case) yb(i,:)]*A_boot + eb(junk+i,:);
                    yb_shock(i+1,:) = [x_boot(junk+i,1:c_case) yb_shock(i,:)]*A_boot + eb_shock(junk+i,:);
                else
                    yb(i+1,:) = [x_boot(junk+i,1:c_case) yb(i,:) Xex(junk+i,:)]*A_boot + eb(junk+i,:);
                    yb_shock(i+1,:) = [x_boot(junk+i,1:c_case) yb_shock(i,:) Xex(junk+i,:)]*A_boot + eb_shock(junk+i,:);
                end
            end
            
        irf = yb_shock - yb;
            
        elseif p == 2       
        % Second obs 
            if isempty(Xex)
                yb(2,:) = [x_boot(junk+1,1:c_case) yb(1,:) x_boot(junk,c_case+1:n+c_case)]*A_boot + eb(junk+1,:);
                yb_shock(2,:) = [x_boot(junk+1,1:c_case) yb_shock(1,:) x_boot(junk,c_case+1:n+c_case)]*A_boot + eb_shock(junk+1,:);
                else
                yb(2,:) = [x_boot(junk+1,1:c_case) yb(1,:) x_boot(junk,c_case+1:n+c_case) Xex(junk+1,:)]*A_boot + eb(junk+1,:);
                yb_shock(2,:) = [x_boot(junk+1,1:c_case) yb_shock(1,:) x_boot(junk,c_case+1:n+c_case) Xex(junk+1,:)]*A_boot + eb_shock(junk+1,:);
            end
            
            for i = 1:IR.per-2
                if isempty(Xex)
                    yb(i+2,:) = [x_boot(junk+i+1,1:c_case) yb(i+1,:) yb(i,:)]*A_boot + eb(junk+i+1,:);
                    yb_shock(i+2,:) = [x_boot(junk+i+1,1:c_case) yb_shock(i+1,:) yb_shock(i,:)]*A_boot + eb_shock(junk+i+1,:);
                else
                    yb(i+2,:) = [x_boot(junk+i+1,1:c_case) yb(i+1,:) yb(i,:) Xex(junk+i,:)]*A_boot + eb(junk+i+1,:);
                    yb_shock(i+2,:) = [x_boot(junk+i+1,1:c_case) yb_shock(i+1,:) yb_shock(i,:) Xex(junk+i+1,:)]*A_boot + eb_shock(junk+i+1,:);
                end
            end
            
        
        irf = yb_shock - yb;
           
        else
                disp('pera tantito')

        end

    %% Sentence with history    
    % GIRF_{y}(n,v_{t},w_{t-i}) = E[Y_{t+n}|v_{t},w_{t-i}] - E[Y_{t+n}|w_{t-i}]
    else
    end

    
    % Colecting IRF
    path_shock(:,ii,:) = irf;
    % Print progress
    fprintf(1,'\b\b\b\b%3.0f%%',100*(ii/GIRF_opt.sim));
end


for i = 1:n
    for j = 1:IR.per
    irf_mean(j,i) = mean(path_shock(j,:,i));
    irf_sup(j,i)  = prctile(path_shock(j,:,i),95);
    irf_inf(j,i)  = prctile(path_shock(j,:,i),5);
    end
end

GIRF.mean = irf_mean;
GIRF.sup = irf_sup;
GIRF.inf = irf_inf;

