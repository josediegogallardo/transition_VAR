function non_linear_test()

% PURPOSE: Peform non-linearity test in a VAR(p) model, by evaluating
% non-linearity in the whole VAR with LR tests and by every single AR
% equation with LM tests.
%
% This test can be performed in an endogenous transition variable, 
% an exogenous transition variable and an independent transition variable 
% (is not part of endogenous or exogenous variables, just marks transition 
% like a markov chain in MS-VAR)
% 
% This functions is based on the opt_lag_translag function written by 
% Saki Bigio (2005)
%--------------------------------------------------------------------------



