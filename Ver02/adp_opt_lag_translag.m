function results=adp_opt_lag_translag(y,nlag,lags2check,x)

% PURPOSE: Evaluates Lagrange Multiplier non-linearity tests and Likelihood 
% Ration tests to every variable set in Y and at differnt lags.
%---------------------------------------------------
% USAGE:  result = opt_lag_translag(y,param,lags2check,x)
% where:    y    = an (nobs x neqs) matrix of y-vectors of dependent
%                  variables
%
%           nlag = the lag length of the VAR
%          
%          param = [nlag shockv trans thresh smooth shock sim IRper history transilag]
%          trans : imagine y and x as a block, and then choose the trans
%          variable position in that block
%
%           (e.g)-> param=[1 2 1 0 1 1 1000 24 30 1]
%
%           x    = optional matrix of variables (nobs x nx)
%                 (NOTE: constant vector automatically included)
%           
%     lags2check = is the largest lag of a transition variable to perform 
%                  the tests on.
%---------------------------------------------------
% RETURNS a structure
% results.LRtests=LRtests;
% results.fstat=fstat    ; 
%---------------------------------------------------
% SEE ALSO: 
% vare, varf, prt_var, prt_granger, prt_ftests (from LeSage's econometrics
% toolbox)
%
% sstvar, stvar_tests (from Saki's Toolkit)
%
%---------------------------------------------------

% written by:
% Saki Bigio (c) 2005
% Econometric Modelling Unit
% Central Bank of Peru
% sakibigio@hotmail.com

%==========================================================================

% Variable Names (to be filled in by users)

names(1).name=['TI'];
names(2).name=['Y'];
names(3).name=['PR'];
names(4).name=['I'];
names(5).name=['BM'];
names(6).name=['ER'];
names(7).name=['D12LQ_US'];

%===========================================================================
param=[nlag 1 4 1 1 1 1 1 1 1] 
[nobs neqs]=size(y) ;
vars=neqs           ;
transilag=lags2check;

LRtests=zeros(transilag,vars);
fstat=zeros(transilag,vars,vars);

if param(3)<=vars
    for i=1:vars
        for j=1:transilag
            param(3)=i  ;
            param(10)=j ;
            if nargin==4
                tests=adp_stvar_tests(y,param,x);
            else
                tests=adp_stvar_tests(y,param);
            end;
            LRtests(j,i)     = tests.LRpval    ;
            fstat(j,i,1:vars)= tests.fstat_pval;
        end;
    end;
else
    for j=1:transilag
        param(10)=j ;
        tests=adp_stvar_tests(y,param,x);
        LRtests(j,1)     = tests.LRpval    ;
        fstat(j,1,1:vars) = tests.fstat_pval;
    end;
end;

results.LRtests=LRtests;
results.fstat=fstat    ; 

disp('Modify this Function to Insert Corresponding variable names')
pause;
% Use the next code if you only use endogenous transition variables (this
% will perform Saki Bigio original tests):
% opt_lag_table;

disp('Hola')
disp(LRtests)

%For exogenous varaibles run:
tab_index = [1;2;3;4;5;6];
testLRvalues = [LRtests(:,1)];
resulttab = [tab_index,testLRvalues];
disp('   lag   LR-test')
disp('--------------')
disp(resulttab)

