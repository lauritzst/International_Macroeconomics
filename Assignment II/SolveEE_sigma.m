function [b_f_opt, e1_opt] =SolveEE_sigma(parameters, n_b_f, n_e)

%{
##########################################################################
This files computes the optimal level of foreign bonds and the optimal
exchange rate by using the system of Euler Equations. It uses a similar
method as the solveEE function to calculate the variations of sigma in task
5i as the solveEE function fails to calculate the different optimal level
of foreign bonds and  exchange rates for varying sigma.
Input variables: 
- parameters:   exogenous model parameters
- n_b_f:        number of grid points for foreign bond
- n_e:          number of grid points for exchange rate
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


% Calculate right-hand side of both EEs -----------------------------------
integral_d = zeros(n_b_f,1); % set up vector
parfor i=1:n_b_f 
    integral_d(i)=integral2(@(y2,e2) fun_dom(y2,e2,b_f(i)),0,Inf,0,Inf);
end
RHS_d = beta.*(1+i_d).*integral_d;  % RHS domestic EE

integral_f = zeros(n_b_f,1); % set up vector
parfor i=1:n_b_f 
    integral_f(i)=integral2(@(y2,e2) fun_for(y2,e2,b_f(i)),0,Inf,0,Inf);
end
RHS_f = beta.*(1+i_f).*integral_f;  % RHS foreign EE


% Solve left-hand sight (LHS) of the domestic EE and optimal foreign bond
% level -------------------------------------------------------------------
b_d_opt = zeros(n_e,1);

for k=1:n_e
    LHS_d = zeros(1,n_b_f);
    LHS_d(1:n_b_f) = u_prime(y1 - e1(k)/p*b_f);
    ddiff = LHS_d - RHS_d';
    falls = (gradient(ddiff,b_f(2)-b_f(1))>0);
    b_d_opt_loop = interp1(ddiff(falls),b_f(falls),0,'linear');
    b_d_opt(k) = b_d_opt_loop;
end


% Solve left-hand sight (LHS) of the foreign EE and optimal foreign bond
% level -------------------------------------------------------------------
b_ff_opt = zeros(n_e,1);

for k=1:n_e
    LHS_f = zeros(1,n_b_f);
    LHS_f(1:n_b_f) = u_prime(y1 - e1(k)/p*b_f)*e1(k);
    fdiff = LHS_f - RHS_f';
    falls = (gradient(fdiff,b_f(2)-b_f(1))>0);
    b_f_opt_loop = interp1(fdiff(falls),b_f(falls),0,'linear');
    b_ff_opt(k) = b_f_opt_loop;
end


% Find equilibrium exchange rate and foreign bond level -------------------
diff_b = b_d_opt - b_ff_opt;        % difference of both optimal et vectors
falls=(gradient(diff_b,b_f(2)-b_f(1))>0);   % find position where slope 
                                            % changes from positive to negative

b_f_opt = interp1(diff_b(falls),b_d_opt(falls),0,'linear'); % optimal foreign bond level
e1_opt=interp1(diff_b(falls),e1(falls),0,'linear'); % optimal exchange rate

end
