function results = linearity_test(Z,p,thv,thp,n1,n2)

% PURPOSE: Evaluates Lagrange Multiplier non-linearity tests and Likelihood 
% Ration tests to every variable set in Y and at differnt lags.
%
%     Z = vector of dependent and exogenous variables Z = [Y X].
%     p = Lag lenght of VAR model.
%   thv = Threshold variable position in Z. 
%   thp = Treshold Maximun lag to test.
%    n1 = Number of dependent variables.
%    n2 = Number of exogenous variables.

% DEPENDENCES:
% stvar_testv2

% Parameters of endogenous VAR equations
[T,k] = size(Z(:,1:n1));

% Empty storage for tests
LRtest = zeros(thp,1);
LRtest_pval = zeros(thp,1);
Fstat=zeros(thp,k);
Fstat_pval=zeros(thp,k);

% Generating results
for i = 1:thp
    thp = i;
    tests=stvar_testsv2(Z,p,thv,thp,n1,n2);
    LRtest(i,1)         = tests.LR;
    LRtest_pval(i,1)    = tests.LRpval;
    Fstat(i,1:k)     = tests.fstat;
    Fstat_pval(i,1:k)= tests.fstatpval;
end
results.LRtest=LRtest;
results.LRtest_pval=LRtest_pval;
results.Fstat=Fstat; 
results.Fstat_pval=Fstat_pval; 












