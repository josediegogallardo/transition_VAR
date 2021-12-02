function results=stvar_testsv2(Z,p,thv,thp,n1,n2)

%     Z = vector of dependent and exogenous variables Z = [Y X].
%     p = Lag lenght of VAR model.
%   thv = Threshold variable position in Z. 
%   thp = Treshold Maximun lag to test.
%    n1 = Number of dependent variables.
%    n2 = Number of exogenous variables.

nlag = p;
results.meth='STVAR_tests';
 
[nobs, neqs] = size(Z(:,1:n1)); % endogenous part values

results.nobs = nobs;         % # of observations
results.neqs = neqs;         % # of equations
results.nlag = nlag;         % # of lags
results.nvar = 2*nlag*n1+1;  % # of variables

%Building the Smooth Transition Function
transi = Z(:,thv);
g = lag(transi,thp);

%Building the non-linear data
gy=[];
ylag=mlag(Z(:,1:n1),nlag);

for i=1:neqs*p
    gy=[gy g.*ylag(:,i)];
end

%Ojo, la reg auxiliar de contraste no incluye exógenas (ver si se las
%incluye o no: Linea 59 aprox
gy=[gy Z(:,n1+1:end)];  % See what happens if we do not include exogenous, see how to
%activate and not

%We should fix the first period cause of the lag and SMTV
varres=vare2(Z(:,1:n1),nlag,gy);

%Loading Variables
n = varres(1).nobs;
k = varres(1).neqs;

%get coefficients
b=[];
for i=1:neqs          
    b  = [b varres(i).beta];
end

%Extract reduced form residuals and (using Cholesky decomposition, obtain
%functional form uncorrelated errores)
res_ur = zeros(n,k);
for i=1:k
    res_ur(:,i) = varres(i).resid;
end

%%compute tests statistics for non-linearity
% LR test for Overall-Non-Linearities
varres_aux=vare2(Z(:,1:n1),nlag,Z(:,n1+1:end)); %Ojo en la auxiliar no incluye las exógenas
res_r=zeros(n,k);
for i=1:k
    res_r(:,i) = varres_aux(i).resid;   %get restricted residuals
end;

T=nobs-nlag;

% VAR-COV of Residuals:
omega0= (1/((T)))*res_r'*res_r;                  %build cov. matrices
omega1= (1/((T)))*res_ur'*res_ur;

% Adding Sim's Correction: 
k=1+n2+2*neqs*nlag;
LR=(T-k)*(log(det(omega0))-log(det(omega1)));
LRpval=1-chis_prb(LR,nlag*neqs^2);

%% Calculating the Bootstraped Equation by Equation LM estimator
FStat_eq=zeros(neqs,1);
for i=1:neqs
    SSR0=sum(res_r(:,i).^2);
    SSR1=sum(res_ur(:,i).^2);
    FStat_eq(i)=((SSR0-SSR1)/(nlag*neqs))/((SSR1/((nobs-nlag)-(2*nlag*neqs+1))));%calculate small sample F-Statistic
end;

pval=1-fdis_prb(FStat_eq,nlag*neqs,(nobs-nlag)-(2*nlag*neqs+1));        %calculate corresponding p-value
results.fstat = FStat_eq;
results.fstatpval=pval;

%% Registering Results
results.omega0=omega0   ;
results.omega1=omega1   ;      
results.LR=LR           ;
results.LRpval=LRpval   ;
results.beta=b  ;
results.g=g;
