%% Clear enviroment and add paths
clear
clc
addpath(genpath('.\var_functions'))

%% Importing data
data = xlsread('data_base.xlsx');
dates   = data(:,1); dates = x2mdate(dates, 0); % Import excel dates
ti      = data(:,2);  % Terms of trade (Annual percentage variance)
pbi     = data(:,3);  % GDP (Annual percentage variance)
pr      = data(:,4);  % IPC (Annual percentage variance)
ir      = data(:,5);  % Monetary Policy Reference Rate
bm      = data(:,6);  % Monetary Base (Annual percentage variance)
er      = data(:,7);  % Bilateral Exchange Rate (Annual percentage variance)
cp      = data(:,8);  % Commodities Price Index (Annual percentage variance)
%% VAR Estimation
y = [ti, pbi, pr, ir, bm, er];
p = 2;
c = 1;
exo_data = [];
VAR_lineal = estimateVAR(y, p, c, exo_data);

%% GIRF
var_results          = VAR_lineal; % Results from linear VAR for bootstrap
IRper                = 24; % Horizonts for IRF's
shockv               = 4; % Position of Endogenous variable to be shocked
GIRF_opt.sim         = 500; % Number of simulations
GIRF_opt.shock       = 1; % Shock size and sign
GIRF_opt.hist        = []; % No history
%GIRF_opt.hist.var    = 2; % Posittion of variable in estimateVAR.X for history dependence
%GIRF_opt.hist.thresh = 0.0; % Threshold to mark history or anothe

GIRF = GIRF_chol(var_results, IRper, shockv, GIRF_opt);