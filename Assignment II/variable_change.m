function [b_f_opt_vec, e1_opt_vec] = variable_change(parameters, n, n_b_f, n_e)

%{
##########################################################################
This files computes the optimal level of foreign bonds and the optimal
exchange rate by using the system of Euler Equations for various values of
the underlying parameters i) risk aversion coefficient, ii) variance of the
exchange rate, iii) variance of the income, and iv) covariance of the 
exchange rate and the income.
Input variables: 
- parameters:   exogenous model parameters
- n:            number of steps
- n_b_f:        number of grid points for foreign bond
- n_e:          number of grid points for exchange rate
Output variables:
- b_f_opt_vec:  vector of optimal level of foreign bonds
- e1_opt_vec:   vector of optimal exchange rate in period1
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
sigma_5 = linspace(1.55,2.45,n);      % risk aversion coefficient
var_e_5 = linspace(0.0001,0.04,n);    % variance exchange variance
var_y_5 = var_e_5;                    % variance income
cov_ye_5 = linspace(-0.0012,0.001,n); % covariance income and exchange rate

% Set up empty vectors ----------------------------------------------------
b_f_opt_vec = zeros(n,4);   % set up grid
e1_opt_vec = zeros(n,4);    % set up grid


% Loops -------------------------------------------------------------------
% Change in the ...
% i) risk aversion coefficient
for i=1:n
    parameters.sigma = sigma_5(i);
    [b_f_opt_vec(i,1),e1_opt_vec(i,1)] = SolveEE_sigma(parameters, n_b_f, n_e);
end
parameters.sigma = 2;   % reset parameter to intial values

% ii) variance of the exchange rate
for i=1:n
    parameters.var_e = var_e_5(i);
    [b_f_opt_vec(i,2),e1_opt_vec(i,2)] = SolveEE(parameters, n_b_f);
end
parameters.var_e = 0.02; % reset parameter to intial values

% iii) variance of the income
for i=1:n
    parameters.var_y = var_y_5(i);
    [b_f_opt_vec(i,3),e1_opt_vec(i,3)] = SolveEE(parameters, n_b_f);
end
parameters.var_y = 0.02; % reset parameter to intial values

% iv) covariance of the exchange rate and the income
for i=1:n
    parameters.cov_ye = cov_ye_5(i);
    [b_f_opt_vec(i,4),e1_opt_vec(i,4)] = SolveEE(parameters, n_b_f);
end
parameters.cov_ye = -0.0001; % reset parameter to intial values

end