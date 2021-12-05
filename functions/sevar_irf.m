function results = sevar_irf(y,nlag,translag,test_type,trans_var,shockv,IRper,thresh,Hist,sim,x)

tic 

results.meth='Self Exciting VAR';

[nobs, neqs] = size(y);

results.nobs = nobs; % # of observations
results.neqs = neqs; % # of equations
results.nlag = nlag; % # of lags

shock = 1;

if nargin == 11
    [nobsx nx] = size(x);
    if (nobsx ~= nobs)
        error('var: nobs in x-matrix not the same as y-matrix');
    end;
    results.nvar = 2*nlag*neqs+nx+2 ; % # of variables per equation
end;

results.nvar = 2*nlag*neqs+2       ; % # of variables per equation

%Building the non-linear data
gy_=[];
ylag=mlag(y,nlag);
if test_type == 1
    % Exogenous # var
    nx = 0;
    % Building the Smooth Transition Function
    transi   = y(:,trans_var);
    transi   = lag(transi,translag);
    for t = 1:nobs
        gfunc(t) = transi(t) > thresh  ; %Weise uses the function with -1/2, note that variables aren't normalized
    end    
    gfunc=double(gfunc)';
    g = gfunc;
    results.g=g  
    
    % Non linear subsistem to test
    for i=1:neqs*nlag
        gy_=[gy_ g.*ylag(:,i)];
    end
    
    % Exogenous variables
    if nargin == 11
        gy_ = [gy_ x];
    end;
      
elseif test_type == 2
    % Exogenous # var
    nx = 1;
    % Building the Smooth Transition Function
    transi   = trans_var;
    transi   = lag(transi,translag);
    for t = 1:nobs
        gfunc(t) = transi(t) > thresh  ; %Weise uses the function with -1/2, note that variables aren't normalized
    end    
    gfunc=double(gfunc)';
    g = gfunc; 
    % Non linear subsistem to test
    for i=1:neqs*nlag
        gy_=[gy_ g.*ylag(:,i)];
    end
    % Exogenous variables
    if nargin == 11
        gy_ = [gy_ x];
    end;
    % Transtition variable as exogenous var
    gy_=[gy_ trans_var];
    
elseif test_type == 3
    % Building the Smooth Transition Function
    transi   = trans_var;
    transi   = lag(transi,translag);
    for t = 1:nobs
        gfunc(t) = transi(t) > thresh  ; %Weise uses the function with -1/2, note that variables aren't normalized
    end    
    gfunc=double(gfunc)';
    g = gfunc;
    % Non linear subsistem to test
    for i=1:neqs*nlag
        gy_=[gy_ g.*ylag(:,i)];
    end
    
    % Exogenous variables
    if nargin == 11
        gy_ = [gy_ x];
    end;
else
    disp('Invalid test type')
end 

%Running the Non-Linear VAR 
varres=vare2(y,nlag,gy_) ;

%Loading Variables for Impulse Response
n = varres(1).nobs      ;
k = varres(1).neqs      ;
gy_=gy_((1+nlag:nobs),:)  ;
g=g((1+nlag:nobs),:)    ;

%Estimation Coefficients 
b=[]                        ;
for i=1:neqs          
    b  = [b varres(i).beta] ;
end

%Estimation T-Stats
tstat=[];
for i=1:neqs          
    tstat  = [tstat varres(i).tstat];
end

%Estimation P-Values
tprob=[];
for i=1:neqs          
    tprob  = [tprob varres(i).tprob];
end

%Estimation Residuals
e = zeros(n,k);
for i=1:k
    e(:,i) = varres(i).resid;
end

% Recovering Predicted Values
yhat=[];
for i=1:neqs          
    yhat  = [yhat varres(i).yhat];
end

% Recovering Alligned Actual Values
yact=[];
for i=1:neqs          
    yact  = [yact varres(i).y]  ;
end

% Recovering R-squared
rsqr=[];
for i=1:neqs          
    rsqr  = [rsqr varres(i).rsqr];
end

% Recovering Alligned Actual Values
yact=[];
for i=1:neqs          
    yact  = [yact varres(i).y]  ;
end

% Recovering R-Squared Adjusted
rbar=[];
for i=1:neqs          
    rbar  = [rbar varres(i).rbar];
end

%registering relevant statistics of the estimation
results.beta  = b       ;
results.tstat = tstat   ;
results.tprob = tprob   ;
results.resid = e       ;
results.yhat  = yhat    ;
results.y     = yact    ;
results.rsqr  = rsqr    ;
results.rbar  = rbar    ;

%Preparing x
if nargin == 3
    x=x((nlag+1:nobs),:);    
end

% Main Loop Starts Here:

    %disp(jj)
    % Bootstrapping the Residuals and Reconstructing Data:
    e_ = bootstrp(1,'equal',e);
    e_ = reshape(e_,n,k);
    y_ = yhat+e_;
    
    gy_=[]            ; 
    ylag_=mlag(y_,nlag);
    gfunc=[];
    % Re-Estimating Parameters with built-in Data:
    %Building the Smooth Transition Function
    if test_type == 1
    % Building the Smooth Transition Function
        transi   = y_(:,trans_var);
        transi   = lag(transi,translag);
        for t = 1:nobs-nlag
            gfunc(t) = transi(t) >= thresh  ; %Weise uses the function with -1/2, note that variables aren't normalized
        end    
        gfunc=double(gfunc)';
        g_ = gfunc; 
    
    % Non linear subsistem to test
        for i=1:neqs*nlag
            gy_=[gy_ g_.*ylag_(:,i)];
        end
    
    % Exogenous variables
        if nargin == 11
            gy_ = [gy_ x];
        end
      
    elseif test_type == 2
        % Building the Smooth Transition Function
        transi   = trans_var;
        transi   = lag(transi,translag);
        for t = 1:nobs-nlag
            gfunc(t) = transi(t) >= thresh  ; %Weise uses the function with -1/2, note that variables aren't normalized
        end
        gfunc=double(gfunc)';
        g_=gfunc;  
        % Non linear subsistem to test
        for i=1:neqs*nlag
            gy_=[gy_ g_.*ylag_(:,i)];
        end
        % Exogenous variables
        if nargin == 11
            gy_ = [gy_ x];
        end;
        % Transtition variable as exogenous var
        gy_=[gy_ trans_var(1:end-nlag,:)];
    
    elseif test_type == 3
        % Building the Smooth Transition Function
        transi   = trans_var;
        transi   = lag(transi,translag);
            for t = 1:nobs-nlag
                gfunc(t) = transi(t) >= thresh  ; %Weise uses the function with -1/2, note that variables aren't normalized
            end    
        gfunc=double(gfunc)';
        g_=gfunc;      
        % Non linear subsistem to test
        for i=1:neqs*nlag
            gy_=[gy_ g_.*ylag_(:,i)];
        end

        % Exogenous variables
        if nargin == 11
            gy_ = [gy_ x];
        end;
    else
        disp('Invalid test type')
    end 

    %Running the Non-Linear VAR 
    varres_=vare2(y_,nlag,gy_);
    
    % Reconstructing rsiduals:
    %Estimation Residuals
    e_ = [];
    for i=1:k
        e_(:,i) = varres_(i).resid;
    end
    
    % Re-constructing values for Beta's
    b_=[];
    for i=1:neqs          
        b_  = [b_ varres_(i).beta];
    end

    % Extracting reduced form residuals and (using Cholesky decomposition, obtain
    % functional form uncorrelated errors from the)
    Omega = (1/(n-k))*e_'*e_;
    C = chol(Omega);
    v = inv(C)*e_';
    v=v';

    %Creating empty matrix to store differenced paths
    paths_dif=[];

    %Starting main simulation loop
for i = 1:sim
        disp(i)
        % Choosing a Random period that coincides with the relevant history
        junk=randi([(IRper+translag+nlag+1) (nobs-(IRper+translag+nlag))]);
        
        if test_type == 1
            if Hist ==1
                while y(junk-translag,trans_var) < thresh
                      junk=randi([(IRper+translag+nlag+1) (nobs-(IRper+translag+nlag))]);
                end
                History=junk;
            else
                while y(junk-translag,trans_var) >= thresh
                      junk=randi([(IRper+translag+nlag+1) (nobs-(IRper+translag+nlag))]);
                end;
                History=junk;
            end;    
        else    
            if Hist ==1
                while trans_var(junk-translag) < thresh
                      junk=randi([(IRper+translag+nlag+1) (nobs-(IRper+translag+nlag))]);
                end
                History=junk;
            else
                while trans_var(junk-translag) >= thresh
                      junk=randi([(IRper+translag+nlag+1) (nobs-(IRper+translag+nlag))]);
                end;
                History=junk;
            end;
        
        end
        %Reshufling uncorrelated errors sim times (montecarlo simulations)
        %and fixing size
        temp = bootstrp(1,'equal',v);
        vbs = reshape(temp,n-nlag,k);

        %Fix a given shock because we are using Cholesky decomposition, a direct addition can be
        %done that will be equal to fixing a "shock" sized standard deviation shock
        %(see Hamilton p. 400)

        vbss=vbs;
        vbss(History,shockv)= shock/(Omega(shockv,shockv)^(1/2)); % this is only the shock, should we add it to the residual...?  
        vbs(History,shockv)= 0;

        %returning dependence to shocks

        ebs_  = (C*vbs')';
        ebss_ = (C*vbss')';
        eb    = [ebs_ ebss_];

        %Setting an auxiliary variable and all complements to rebuild series
        paths=[];
        for i=0:neqs:neqs

            ebs=eb(:,i+1:i+neqs);
            path = y_;
            g_sim=gfunc(1+nlag:end,:);

            ylag = mlag(y_,nlag);
            ymat=ylag(nlag+1:nobs-nlag,:);

            gy_=[];
            for i=1:neqs*nlag   
                gy_=[gy_ g_sim.* ymat(:,i)]; 
            end   

            if nargin == 11
                gy_=[gy_ x(nobs-2*nlag,:)];    
            end; 
            
            if nargin == 11
                gy_=[gy_ x(nobs-2*nlag,:)];    
            end; 
            
            if test_type ==  2
                gy_ = [gy_ trans_var(1:nobs-2*nlag,:)];
            else
                gy_ = gy_;
            end

            gy_=[ymat gy_ ones(nobs-2*nlag,1)];

            %Rebuilding series for I-R period
            for t=History:History+IRper % (le quitamos 
                path(t,:)= gy_(t,:)*b_ + ebs(t,:); % ---> Beta from this loop.

                %updating gy matrix with the values just obtained
                for j = 0:k-1
                    ymat(t+1,j*nlag+1:(j+1)*nlag)=rot90(path((t-nlag+1):t,j+1),-1); %check this out (t or t+1)
                end;

                %Updating point t of the G function 
                if test_type == 1
                    g_sim(t) = double(path(t-translag,trans_var)>thresh);
                else
                    g_sim(t) = g_sim(t);
                end
                %Emptying GY matrix and updating it
                gy_=[];
                a = 1;
                
                for i=1:neqs*nlag
                    gy_=[gy_ g_sim.*ymat(:,i)];
                end;    
                    
                if nargin == 11
                    gy_=[gy_ x];    
                end;
                
                if test_type ==  2
                    gy_ = [gy_ trans_var(1:nobs-2*nlag,:)];
                else
                    gy_ = gy_;
                end
                
                gy_=[ymat gy_ ones(nobs-2*nlag,1)]; 
            end;

            path=path((History):(History+IRper),:); 

            % Store resulting path
            paths=[paths path];   
        end;

        paths_us  = paths(:,1:neqs)             ;
        paths_s   = paths(:,neqs+1:2*neqs)      ;
        paths_dif = [paths_dif paths_s-paths_us];
end;  
    
    
%Reordering main loop's output
for j=1:neqs
    for i=1:sim
    IRFS_aux(:,i+(j-1)*sim)=paths_dif(:,j+neqs*(i-1));
    end;
end;


%Building I-R and confidence bands over the 90%
inf(IRper+1,neqs)=0;
sup(IRper+1,neqs)=0;
mode(IRper+1,neqs)=0;

for i=1:neqs
    for t=1:IRper+1
        IRF(t,i)=mean(IRFS_aux(t,(1+(i-1)*sim):i*sim))  ;
        % IRF_mode(t,i)=mode(IRFS_aux(t,(1+(i-1)*sim):i*sim));
        [prob val]=ecdf(IRFS_aux(t,(1+(i-1)*sim):i*sim));
        inf(t,i)=val(min(find(prob>0.05)))              ;
        sup(t,i)=val(min(find(prob>=0.95)-1))           ;
    end;
end;
results.irfs=IRF        ;
%results.irfs_mode=IRF_mode;
results.irfs_sup=sup    ;
results.irfs_inf=inf    ;



