function results=stvar_tests_gfunc(y,param,x)

% PURPOSE: performs a smooth transition vector autoregression
% only presenings tests
%---------------------------------------------------
% USAGE:  result = stvar_tests_gfunc(y, param, x)
% where:    y    = an (nobs x neqs) matrix of y-vectors
%        
%           y should be fixed from endog to exog
%
%           Param is a 1x9 vector that includes the following 
%           information in order:
%
%           nlag = the lag length
%         shockv = variable of y being shocked (column of y)
%          trans = position of the transition variable in the y matrix
%         thresh = the threshold value
%         smooth = the smoothnes parameter                
%          shock = the standard deviation shocks 
%            sim = the number of montecarlo simulations
%          IRper = the number of periods for the impulse response function
%        history = the period we want to atart with        
%
%          param = [nlag shockv trans thresh smooth shock sim IRper history]
%
%           (e.g)-> param=[1 2 1 0 1 1 1000 24 30]
%
%           x    = optional matrix of variables (nobs x nx)
%                 (NOTE: constant vector automatically included)
%
%       NOTE: Smooth transition function is not applied to exogenous
%             variables
%
%
%---------------------------------------------------
% RETURNS a structure
% results.meth = 'STVAR'
% results.nobs = nobs, # of observations
% results.neqs = neqs, # of equations
% results.nlag = nlag, # of lags
% results.nvar = nlag*neqs+nx+1, # of variables per equation
% --- the following are referenced by equation # --- 
% results(eq).beta  = bhat for equation eq
% results(eq).tstat = t-statistics 
% results(eq).tprob = t-probabilities
% results(eq).resid = residuals 
% results(eq).yhat  = predicted values 
% results(eq).y     = actual values 
% results(eq).sige  = e'e/(n-k)
% results(eq).rsqr  = r-squared
% results(eq).rbar  = r-squared adjusted
% results(eq).boxq  = Box Q-statistics
% results(eq).ftest = Granger F-tests
% results(eq).fprob = Granger marginal probabilities
%---------------------------------------------------
% SEE ALSO: 
%vare, varf, prt_var, prt_granger, prt_ftests (from LeSage's econometrics
% toolbox)
%
%stvar2, stvar_search (from Saki's Toolkit)
%
%---------------------------------------------------

% written by:
% Saki Bigio, Macroeconomic Analysis Dept.
% Banco Central de Reserva del Peru
% Paul de Beaudiez 530,
% Lima L27,  PERU
% sbigio@bcrp.gob.pe

results.meth='STVAR_tests'  ;

[nobs neqs] = size(y)       ;

results.nobs = nobs         ; % # of observations
results.neqs = neqs         ; % # of equations

%Setting Default Values
nlag   = param(1);
if nlag ==0
    nlag=1                  ;
end;

results.nlag = nlag         ; % # of lags

shockv = param(2)           ;
if shockv == 0
    shockv = neqs           ;
end;

trans  = param(3);
if trans == 0
    trans = 1;
end;

thresh = param(4);
if thresh == 0
    thresh = 0;
end;

smooth = param(5);
if smooth == 0
    smooth = 7;
end;

shock  = param(6);
if shock == 0
    shock = 1;
end;

sim    = param(7);
if sim == 0
    sim   = 1000;
end;

IRper  = param(8);
if IRper == 0
    IRper = 24;
end;

Hist   = param(9);
if Hist == 0
    Hist = nlag+1;
end;

translag = param(10);
if translag == 0
    translag = 1;
end;

if nargin == 3
    [nobsx nx] = size(x);
if (nobsx ~= nobs)
    error('var: nobs in x-matrix not the same as y-matrix');
end;
    results.nvar = 2*nlag*neqs+nx+2             ; % # of variables per equation
end;

results.nvar = 2*nlag*neqs+2*1              ; % # of variables per equation

%Building the Smooth Transition Function
transi   = y(:,trans)                                    ;
standard = std(transi)                                   ;
transi = lag(transi,translag)                            ;
gfunc    =(1+exp(-smooth*(transi-thresh)/standard)).^(-1)         ;  % Does not Include Standard Deviation
g        = gfunc                                         ;

%Generalizing nx for case when we have exog. var.
nx = 0                                                   ;

%Building the non-linear data

gy=[]            ; 
ylag=mlag(y,nlag);

for i=1:neqs*nlag
    gy=[gy g.*ylag(:,i) ];
end    

if nargin == 3
    gy=[gy x]            ;
end;  

%We should fix the first period cause of the lag and SMTV
varres=vare2(y,nlag,gy) ;

%Loading Variables for Impulse Response
n = varres(1).nobs      ;
k = varres(1).neqs      ;

%get coefficients %%%%%%Generalize
b=[]                        ;
for i=1:neqs          
    b  = [b varres(i).beta] ;
end

%Extract reduced form residuals and (using Cholesky decomposition, obtain
%functional form uncorrelated errores)
e = zeros(n,k)                  ;
for i=1:k
    e(:,i) = varres(i).resid    ;
end

%compute tests statistics for non-linearity
%start with the equation by equation test
varres_aux=vare2(y,nlag)        ;

res_r=zeros(n,k);
for i=1:k
    res_r(:,i) = varres_aux(i).resid;   %get restricted residuals
end;
res_ur=e;                                %get unrestricted residuls

%Calculating the Equation by Equation LM estimator
FStat_eq=zeros(neqs,1);
for i=1:neqs
    SSR0=sum(res_r(:,i).^2);
    SSR1=sum(res_ur(:,i).^2);
    FStat_eq(i)=((SSR0-SSR1)/(nlag*neqs))/((SSR1/((nobs-nlag)-(2*nlag*neqs+1))));%calculate small sample F-Statistic
end;

pval=1-fdis_prb(FStat_eq,nlag*neqs,(nobs-nlag)-(2*nlag*neqs+1));        %calculate corresponding p-value
results.fstat_pval=pval;

%the overall LR-test

%omega0= (1/((nobs-nlag)-neqs*nlag-1))*res_r'*res_r;                  %build cov. matrices
%omega1= (1/((nobs-nlag)-2*neqs*nlag-1))*res_ur'*res_ur;

%LR=((nobs-nlag)-nlag*neqs-1)*(log(det(omega0))-log(det(omega1)));
%LRpval=1-chis_prb(LR,nlag*neqs^2);

omega0= (1/((nobs-nlag)))*res_r'*res_r;                  %build cov. matrices
omega1= (1/((nobs-nlag)))*res_ur'*res_ur;

T=nobs-nlag;

% Adding Sim's Correction: 
k=1+nx+2*neqs*nlag;
LR=(T-k)*(log(det(omega0))-log(det(omega1)));
LRpval=1-chis_prb(LR,nlag*neqs^2);

results.omega0=omega0;
results.omega1=omega1;      
results.LR=LR        ;
results.LRpval=LRpval;
results.beta=b       ;
results.g=g          ;
