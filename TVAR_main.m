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