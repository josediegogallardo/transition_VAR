function results = non_linear_test(y,nlag,lags2check,test_type,trans_var)

% PURPOSE: Peform non-linearity test in a VAR(p) model, by evaluating
% non-linearity in the whole VAR with LR tests and by every single AR
% equation with LM tests.
%
% This test can be performed in an endogenous transition variable,or in 
% an exogenous transition variable or an independent transition variable 
% (is not part of endogenous or exogenous variables, just marks transition 
% like a markov chain in MS-VAR)
% 
% This functions is based on the opt_lag_translag function written by 
% Saki Bigio (2005) and Winkelried (2003) assimetric pash trough paper.
%--------------------------------------------------------------------------

[nobs neqs] = size(y);
results.nobs = nobs; % # of observations
results.neqs = neqs; % # of equations

nlag = nlag;
results.nlag = nlag; % # of lags in VAR(p)

lags2check = lags2check;
results.lags2check = lags2check; % # of lags to test non-linearity

if test_type == 1

    
elseif test_type == 2
   performar test exogeno
elseif test_type == 3
   performar test independiente
else
   disp('Invalid type of test') 
end
    


