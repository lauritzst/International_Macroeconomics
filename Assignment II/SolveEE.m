function [b_f_opt, e1_opt] = SolveEE(parameters,n_b_f)

%{
##########################################################################
This files computes the optimal level of foreign bonds and the optimal
exchange rate by using the system of Euler Equations.
Input variables: 
- parameters:   exogenous model parameters
- n_b_f:        number of grid points for foreign bond
Output variables:
- b_f_opt:  optimal level of foreign bonds
- e1_opt:   optimal exchange rate in period1
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



% Set up Anonymous Functions ----------------------------------------------
u_prime = @(c) c.^(-sigma); % marginal utility

rho = cov_ye/(sqrt(var_y)*sqrt(var_e));
f_ye = @(y, e, var_e, var_y, cov_ye) (exp(-((((log(y)-mu_y)./(sqrt(var_y))).^2-2.*rho.* ...
    ((log(y)-mu_y)./(sqrt(var_y))).*((log(e)-mu_e)./(sqrt(var_e)))+ ...
    ((log(e)-mu_e)./(sqrt(var_e))).^2)./(2.*(1-rho.^2))))./ ...
    (2.*pi.*sqrt(var_y).*sqrt(var_e).*sqrt(1-rho.^2).*y.*e));   % bivariate distribution

fun_dom=@(y2,e2,b_f) u_prime(y2 + e2./p.*(1+i_f).*b_f).* ...
    f_ye(y2,e2, var_e, var_y, cov_ye);    % set up function for domestic EE
fun_for=@(y2,e2,b_f) e2.*u_prime(y2 + e2./p.*(1+i_f).*b_f).* ...
    f_ye(y2,e2, var_e, var_y, cov_ye);  % set up function for foreign EE

% Rearrange both EE for e1 ------------------------------------------------
syms e_d(bf,rhsd) e1_sol bf rhsf 
e_d = p./bf.*y1 - p./bf.*(rhsd)^(-1./sigma); 
e_d = matlabFunction(e_d);      % Anonymous Function for domestic optimal exchange rate
e_f = solve(e1_sol.*u_prime(y1 - e1_sol./p.*bf) == rhsf, e1_sol);
e_f = matlabFunction(e_f(1));       % Anonymous Function for foreign optimal exchange rate

% Calculate Integrals (Expectations) on the RHS of EE ---------------------
integral_d = zeros(n_b_f,1); % set up vector
integral_f = zeros(n_b_f,1); % set up vector
% Loop to calculate Integrals for varying foreign bond levels
parfor i=1:n_b_f 
    integral_d(i)=integral2(@(y2,e2) fun_dom(y2,e2,b_f(i)),0,Inf,0,Inf);
    integral_f(i)=integral2(@(y2,e2) fun_for(y2,e2,b_f(i)),0,Inf,0,Inf);
end
RHS_d = beta.*(1+i_d).*integral_d;  % RHS domestic EE
RHS_f = beta.*(1+i_f).*integral_f;  % RHS foreign EE


% Calulate optimal exchange rate for variable foreign bonds ---------------
e_d_opt = e_d(b_f, RHS_d);
e_f_opt = e_f(b_f, RHS_f);


% Find equilibrium foreign bond -------------------------------------------
diff_b = e_d_opt - e_f_opt;                 % difference of both optimal et vectors
falls=(gradient(diff_b,b_f(2)-b_f(1))<0);   % find position where slope changes from 
                                            % positive to negative

b_f_opt=interp1(diff_b(falls),b_f(falls),0,'linear'); % optimal foreign bond level


% Find equilibrium exchange rate ------------------------------------------
fun_dom_opt=@(y2,e2) u_prime(y2 + e2./p.*(1+i_f).*b_f_opt).* ...
    f_ye(y2,e2, var_e, var_y, cov_ye); % Anonymous function with equilibrium foreign bond
integral_d_opt=integral2(@(y2,e2) fun_dom_opt(y2,e2),0,Inf,0,Inf); % Solve integral with equilibrium foreign bond
RHS_d = beta.*(1+i_d).*integral_d_opt; % RHS of domestic EE with equilibrium foreign bond

e1_opt = p./b_f_opt.*y1 - p./b_f_opt.*(RHS_d)^(-1./sigma); % optimal exchange rate

end
