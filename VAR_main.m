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
exo_data = [cp];
VAR_lineal = estimateVAR(y, p, c, exo_data);

%% GIRF
var_results          = VAR_lineal; % Results from linear VAR for bootstrap
IR.per               = 24; % Horizonts for IRF's
IR.shocktype         = 0; % 0: custom shock in GIRF_opt.shock, 1: one s.d.
shockv               = 4; % Position of Endogenous variable to be shocked
GIRF_opt.sim         = 2000; % Number of simulations
GIRF_opt.shock       = 0.25; % Shock size and sign
GIRF_opt.hist        = []; % No history
%GIRF_opt.hist.var    = 2; % Posittion of variable in estimateVAR.X for history dependence
%GIRF_opt.hist.thresh = 0.0; % Threshold to mark history or anothe

GIRF = GIRF_chol(var_results, IR, shockv, GIRF_opt);

vnames =  { 'Terms of trade',
            'Output ',    
            'Inflation',    
            'Interest rate',
            'Monetary agregate',
            'Exchange rate'};
            
periods=1:IR.per;
zero_line = zeros(IR.per,1);
   figure(7)
for j=1:var_results.n
    % Subplot 1
    subplot(var_results.n/2,var_results.n/2,j);
    plot(periods,GIRF.sup(:,j),'b:',...
         periods,GIRF.mean(:,j),'b-',...
         periods,GIRF.inf(:,j),'b:',...
         periods,zero_line,'k'); 
    xlim([1 IR.per]);
    xticks(0:IR.per/6:IR.per+1);
    grid on;
    hold on;
    
    xlabel('Periods');
    ylabel(vnames{j});
end