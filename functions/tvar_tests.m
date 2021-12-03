function results = tvar_tests(y,nlag,transilag,test_type,trans_var)
    
[nobs, neqs] = size(y);

nvar = 2*nlag*neqs+2*1; % # of variables per equation

nx = 0;

%Building the Smooth Transition Function
transi   = y(:,trans_var);
transi   = lag(transi,transilag);
g        = transi;

%Building the non-linear data
gy=[];
ylag=mlag(y,nlag);

for i=1:neqs*nlag
    gy=[gy g.*ylag(:,i)];
end

varres=vare2(y,nlag,gy);
    
    %Loading Variables for Impulse Response
n = varres(1).nobs;
k = varres(1).neqs;

%get coefficients %%%%%%Generalize
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
varres_aux=vare2(y,nlag);
res_r=zeros(n,k);
for i=1:k
    res_r(:,i) = varres_aux(i).resid;   %get restricted residuals
end;

T=nobs-nlag;

% VAR-COV of Residuals:
omega0= (1/((T)))*res_r'*res_r;                  %build cov. matrices
omega1= (1/((T)))*res_ur'*res_ur;

% Adding Sim's Correction: 
k=1+nx+2*neqs*nlag;
LR=(T-k)*(log(det(omega0))-log(det(omega1)));
LRpval=1-chis_prb(LR,nlag*neqs^2);

%% Calculating the Bootstraped Equation by Equation LM estimator
FStat_eq=zeros(neqs,1);
for i=1:neqs
    SSR0=sum(res_r(:,i).^2);
    SSR1=sum(res_ur(:,i).^2);
    FStat_eq(i)=((SSR0-SSR1)/(nlag*neqs))/((SSR1/((nobs-nlag)-(2*nlag*neqs+1))));%calculate small sample F-Statistic
end;
results.fstat=FStat_eq;

pval=1-fdis_prb(FStat_eq,nlag*neqs,(nobs-nlag)-(2*nlag*neqs+1));        %calculate corresponding p-value
results.fstat_pval=pval;

%% Registering Results
results.omega0=omega0   ;
results.omega1=omega1   ;      
results.LR=LR           ;
results.LRpval=LRpval   ;
results.beta=b  ;
results.g=g;
end