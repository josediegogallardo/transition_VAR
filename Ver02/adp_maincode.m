%% [0] Pseudo code
% 1. Test de no linealidad
% 2. Grid search
% 3. STVAR estimation
% 4. IRFS cumputation
clear
clc
%% [1] TEST DE NO LINEALIDAD
% Se busca performar el test de no linealidad incluso si queremos que la
% varaible de transición sea endógena o exógena.

% Para el desarrollo se va a acordar que la variable Discount Rate va a
% hacer de variable exógena.

load peru.data;% a test data set containing
               % data from the Peruvian Economy:
               % Discount Rate (level), Real Exchange Rate (Log Annual Dif.)
               % GDP (Log Annual Dif.), Price Index (Log Annual Dif.) 
               
% Non linearity test
% In adp_opt_laf_translag: you must set up variables names (line 46)
%                          you must set the params correctly (line 55)
    nlag        = 2;
    y           = peru(:,2:end);
    x           = peru(:,1);
    lags2check  = 6;
    adp_opt_lag_translag(y,nlag,lags2check,x)