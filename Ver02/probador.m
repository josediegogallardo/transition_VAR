clear;
clc;

load peru.data;
% a test data set containing
% data from the Peruvian Economy:
% Discount Rate (level), Real Exchange Rate (Log Annual Dif.)
% GDP (Log Annual Dif.), Price Index (Log Annual Dif.) 

results = linearity_test(peru,2,3,6,4,0);
% results = linearity_test(Z,p,thv,thp,n1,n2);
% Evaluar si los p-val mejoran quitando el control por exógenas
%     Z = vector of dependent and exogenous variables Z = [Y X].
%     p = Lag lenght of VAR model.
%   thv = Threshold variable position in Z. 
%   thp = Treshold Maximun lag to test.
%    n1 = Number of dependent variables.
%    n2 = Number of exogenous variables.
param = [2,3,2,4,0];
% param = [p,thv,thp,n1,n2]
%     p = Lag lenght of VAR model.
%   thv = Threshold variable position in Z. 
%   thp = Treshold lag to test.
%    n1 = Number of dependent variables.
%    n2 = Number of exogenous variables.
results2 = optsmooth(peru,param, step_smooth, points_smooth, tol_thresh, points_thresh);
%              Z = ya tu sa
%          param = previous important parameters
%    step_smooth = size of jumps in  
%  points_smooth = number of jumps
%     tol_thresh = range from highest to lowest point in transition variable
%  points_thresh = points in range 

           