function RP_5_vec = RP_5(parameters, n, e_vec)

%{
##########################################################################
This files computes the optimal risk premium by using the system of Euler 
Equations.
Input variables: 
- parameters:   exogenous model parameters
- n:            number of steps
- e_vec:        vector of optimal exchange rate in period1
Output variables:
- RP_5_vec:     vector of optimal risk premia
##########################################################################
%}

% Unpack Structure --------------------------------------------------------
beta = parameters.beta;     % discount factor
sigma = parameters.sigma;   % risk aversion coefficient
i_d = parameters.i_d;       % domestic interest rate
i_f = parameters.i_f;       % foreign interest rate
y1 = parameters.y1;         % income in t
p = parameters.p;           % price

mu_y = parameters.mu_y;     % mean of log-normal distribution (income)
mu_e = parameters.mu_e;     % mean of log-normal distribution (exchange rate)
var_y = parameters.var_y;   % variance of log-normal distribution (income)
var_e = parameters.var_e;   % variance of log-normal distribution (exchange rate)
cov_ye = parameters.cov_ye; % covariance of income and exchange rate
b_f = parameters.b_f;       % grid for foreign bonds
e1 = parameters.e1;         % grid for exchange rate

% Create Grids ------------------------------------------------------------
var_e_5 = linspace(0.0001,0.04,n);    % variance exchange variance
var_y_5 = var_e_5;                    % variance income
cov_ye_5 = linspace(-0.0012,0.001,n); % covariance income and exchange rate

% Set up Anonymous Functions ----------------------------------------------
f_ye = @(y, e, var_e, var_y, cov_ye, rho) (exp(-((((log(y)-mu_y)./(sqrt(var_y))).^2-2.*rho.* ...
    ((log(y)-mu_y)./(sqrt(var_y))).*((log(e)-mu_e)./(sqrt(var_e)))+ ...
    ((log(e)-mu_e)./(sqrt(var_e))).^2)./(2.*(1-rho.^2))))./ ...
    (2.*pi.*sqrt(var_y).*sqrt(var_e).*sqrt(1-rho.^2).*y.*e)); % bivariate distribution

rho1 = cov_ye/(sqrt(var_y)*sqrt(var_e));

% Set up empty vectors ----------------------------------------------------
RP_5_vec = zeros(n,4);

for i=1:n
    % Variation in sigma
    expect_e2 = integral2((@(y2,e2) e2./e_vec(i,1).* ...
        f_ye(y2,e2,var_e,var_y,cov_ye, rho1)),0,Inf,0,Inf); % expected e2
    RP_5_vec(i,1) = 1+i_d-(1+i_f)*expect_e2;
    % Variation in the exchange variance
    rho2 = cov_ye/(sqrt(var_e_5(i))*sqrt(var_e));   % rho_ye
    expect_e2 = integral2((@(y2,e2) e2./e_vec(i,2).* ...
        f_ye(y2,e2, var_e_5(i),  var_y, cov_ye, rho2)),0,Inf,0,Inf);% expected e2
    RP_5_vec(i,2) = 1+i_d-(1+i_f)*expect_e2;
    % Variation in the income variance
    rho3 = cov_ye/(sqrt(var_y)*sqrt(var_e));        % rho_ye
    expect_e2 = integral2((@(y2,e2) e2./e_vec(i,3).* ...
        f_ye(y2,e2, var_e, var_y_5(i), cov_ye, rho3)),0,Inf,0,Inf); % expected e2
    RP_5_vec(i,3) = 1+i_d-(1+i_f)*expect_e2;
    % Variation in the exchange rate and the income
    rho4 = cov_ye_5(i)/(sqrt(var_y)*sqrt(var_e));   % rho_ye
    expect_e2 = integral2((@(y2,e2) e2./e_vec(i,4).* ...
        f_ye(y2,e2, var_e, var_y, cov_ye_5(i), rho4)),0,Inf,0,Inf); % expected e2
    RP_5_vec(i,4) = 1+i_d-(1+i_f)*expect_e2;
end
end