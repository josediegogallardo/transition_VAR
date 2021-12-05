% Clean data and enviroment
clear;
clc
format short

% Add to path the fucntions for this code, note that we depend from variuos
% function from the Spatial Toolbox by James P. LeSage
addpath(genpath('.\functions'))

% Importing data
data = xlsread('data_base.xlsx');
dates   = data(:,1); dates = x2mdate(dates, 0); % Import excel dates
ti      = data(:,2);  % Terms of trade (Annual percentage variance)
pbi     = data(:,3);  % GDP (Annual percentage variance)
pr      = data(:,4);  % IPC (Annual percentage variance)
ir      = data(:,5);  % Monetary Policy Reference Rate
bm      = data(:,6);  % Monetary Base (Annual percentage variance)
er      = data(:,7);  % Bilateral Exchange Rate (Annual percentage variance)
cp      = data(:,8);  % Commodities Price Index (Annual percentage variance)

y = [ti, pbi, pr, ir, bm, er];
trans_var = cp;
nlag = 1;
lags2check = 6;
test_type = 3;
results = non_linear_test(y,nlag,lags2check,test_type,trans_var);

lags2tab = 1:1:lags2check;
res = [lags2tab' results.LMtests, results.LRtests'];
disp(res)

shockv = 4; % variable position to shock
IRper = 24;
thresh = 5.275;
sim = 2000;
translag = 1;

irf_low = sevar_irf(y,nlag,translag,test_type,trans_var,shockv,IRper,thresh,0,sim);
irf_high = sevar_irf(y,nlag,translag,test_type,trans_var,shockv,IRper,thresh,1,sim);

neqs = 6;
irfh = irf_high.irfs;
suph = irf_high.irfs_sup;
infh = irf_high.irfs_inf;
irfl = irf_low.irfs;
supl = irf_low.irfs_sup;
infl = irf_low.irfs_inf;

vnames =  { 'Terms of trade',
            'Output ',    
            'Inflation',    
            'Interest rate',
            'Monetary agregate',
            'Exchange rate'};
            
   periods=1:IRper+1;
   figure(7)
for j=1:neqs
    subplot(neqs/2,neqs/2,j);
    plot(periods,suph(:,j),'b:',periods,infh(:,j),'b:',periods,irfh(:,j),'b-',periods,0,'k'); hold on
    plot(periods,supl(:,j),'r:',periods,infl(:,j),'r:',periods,irfl(:,j),'r-',periods,0,'k'); 
    xlabel('Periods');
    ylabel(vnames{j});
end
