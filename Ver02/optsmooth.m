function results = optsmooth(y, param, step_smooth, points_smooth, tol_thresh, points_thresh, x)

% PURPOSE: Performs a grid search in order to find optimal threshold and 
%---------------------------------------------------
% USAGE:  result = optsmooth(y, param, step_smooth, points_smooth, tol_thresh, points_thresh, x)
% where:    y    = an (nobs x neqs) matrix of y-vectors
%        
%           y should be fixed from most endogenous to most exogenous
%
%           Param is a 1x10 vector that includes the following 
%           information in order:
%
%           nlag = the lag length of the VAR structure
%         shockv = won't matter for function
%          trans = position of the transition variable in the y matrix
%         thresh = won't matter for function
%         smooth = won't matter for function                
%          shock = won't matter for function
%            sim = won't matter for function
%          IRper = won't matter for function
%        history = won't matter for function
%       translag = the lag of the transition variable  
%
%           (e.g.)-> param=[1 2 1 0 1 1 1000 24 1 1]
%
%    step_smooth = size of jumps in  
%  points_smooth = number of jumps
%     tol_thresh = range from highest to lowest point in transition variable
%  points_thresh = points in range 
%
%           (e.g.)-> param=[2 2 2 0 1 1 1000 24 1 2]
%
%           x    = optional matrix of variables (nobs x nx)
%                 (NOTE: constant vector automatically included)
%
% NOTE: Smooth transition function is not applied to exogenous
% variables
%
%---------------------------------------------------
% RETURNS a structure
% results.smooth        = optimal smoothness parameter ; 
% results.smooth_line   = grid point for the smoothnes ;
% results.threshold     = optimal threshold value      ;
% results.threshold_line= grid points for the threshold values;
% results.LRpval        = minimu level LRpval                       ;
% results.pval          = P-Value for the Optimal set of parameters;
% results.start_thresh  = the starting threshold value;
%
%---------------------------------------------------
% SEE ALSO: vare, varf, prt_var, prt_granger, prt_ftests (from LeSage's econometrics
% toolbox) and optlagtranslag, stvar, sstvar (from Saki's toolkit)
%
%---------------------------------------------------
% Programmed by Saki Bigio (c) 2006
% Econometric Modelling Unit
% Central Reserve Bank of Peru

[nobs neqs]=size(y)                                         ;
mini=min(y(:,param(3)))                                     ;
maxi=max(y(:,param(3)))                                     ;
start_thresh=mini+(maxi-mini)*tol_thresh/2                  ;
step_thresh=(maxi-mini)*(1-tol_thresh)/points_thresh        ;
temppval=1                                                  ;
smooth1 = 0                                                 ;
param(7)= 1                                                 ;
param(8)= 1                                                 ;
param(9)= 1+param(1)                                        ;

for j = 1:points_smooth
                smooth_line(j)=1+(j-1)*step_smooth          ;
end;    

for i = 1:points_thresh
            param(4)=start_thresh+i*step_thresh             ;
            threshold_line(i)=start_thresh+(i-1)*step_thresh;
            for j = 1:points_smooth
                param(5)=1+(j-1)*step_smooth                ;
                
                if nargin == 7 
                    stvar1 = stvar_tests_gfunc(y, param, x) ;
                else
                    stvar1 = stvar_tests_gfunc(y, param)    ;
                end
        
                LRpval(j,i)= stvar1.LRpval                  ;
                if LRpval(j,i) < temppval
                    temppval=LRpval(j,i)                    ;
                    smooth=j*step_smooth                    ;
                    threshold=start_thresh+i*step_thresh    ;
                end;
            end;
end;

results.smooth=smooth                       ;
results.smooth_line=smooth_line             ;
results.threshold=threshold                 ;
results.threshold_line=threshold_line       ;
results.LRpval=LRpval                       ;
results.pval=temppval                       ;
results.start_thresh=start_thresh           ;
results.min=mini                            ;
results.max=maxi                            ;

surf(threshold_line,smooth_line,LRpval)     ;
axis tight;
ylabel('smoothness parameter (\lambda)')    ;
xlabel('threshold value (\theta)')          ;    
zlabel('LR P-Value')                        ;

orient landscape;
print -dpsc2  -loose -r300 -f 'Parameter Grid Search'