function RP_i_vec = RP_change(parameters,n,n_b_f)

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

% Calculate Risp Premia ------------------------—------------------------—-

u_prime = @(c) c.^(-sigma); % marginal utility

var_e7 = linspace(0.02,0.001,n);
RP_i_vec = zeros(n,1); 
RP_ii_vec = zeros(n,1); 

for i=1:n
    % i)
    rho = cov_ye/(sqrt(var_y)*sqrt(var_e7(i)));
    f_ye = @(y, e, var_e, var_y, cov_ye) (exp(-((((log(y)-mu_y)./(sqrt(var_y))).^2-2.*rho.* ...
        ((log(y)-mu_y)./(sqrt(var_y))).*((log(e)-mu_e)./(sqrt(var_e)))+ ...
        ((log(e)-mu_e)./(sqrt(var_e))).^2)./(2.*(1-rho.^2))))./ ...
        (2.*pi.*sqrt(var_y).*sqrt(var_e).*sqrt(1-rho.^2).*y.*e)); % bivariate distribution
    parameters.var_e = var_e7(i);
    [~,e1_opt7] = SolveEE(parameters, n_b_f);
    % [b_f_opt7,e1_opt7] = SolveEE(parameters, n_b_f);
    expect_e2 = integral2((@(y2,e2) e2./e1_opt7.* ...
        f_ye(y2,e2, var_e7(i), var_y, cov_ye)),0,Inf,0,Inf);
    RP_i_vec(i,1) = 1+i_d-(1+i_f)*expect_e2;
    
    % Check values
    %{
    % ii)
    c1 = y1 - 1/p*e1_opt7*b_f_opt7;
    M = @(c2) beta.*c2./u_prime(c1);
    expect_M = integral2((@(y2,e2) M(u_prime(y2 + e2./p.*(1+i_f).*b_f_opt7)).* ...
        f_ye(y2,e2, var_e7(i), var_y, cov_ye)),0,Inf,0,Inf);
    cov_Me = integral2((@(y2,e2) M(u_prime(y2 + e2./p.*(1+i_f).*b_f_opt7)).* ...
        e2./e1_opt7.* f_ye(y2,e2, var_e7(i), var_y, cov_ye)),0,Inf,0,Inf) - expect_M.*expect_e2;
    RP_ii_vec(i,1) = (1+i_f).*cov_Me./expect_M;
    %}
end

end