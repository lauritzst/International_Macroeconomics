function [b_f_opt, e1] = Fx_solveEE(parameters,n_b_f)
%{
##########################################################################
This files computes the optimal level of foreign bonds by using the system
of Euler Equations.
Input variables: 
- parameters: exogenous model parameters
- n_b_f: number of grid points
Output variables:
- b_f_opt: optimal level of foreign bonds
- e1: exchange rate in period1
##########################################################################
%}

% unpack structure
%- - - - - - - - - - - - - - - - 
beta = parameters.beta; % discount factor
sigma = parameters.sigma; % risk aversion coefficient
i_d = parameters.i_d; % domestic interest rate
i_f = parameters.i_f; % foreign interest rate
e2 = parameters.e2; % exchange rate in t+1
y1 = parameters.y1; % income in t
p = parameters.p; % price
mu_y = parameters.mu_y; % mean of log-normal distribution
sd_y = parameters.sd_y; % standard deviation of log-normal distribution
b_f = parameters.b_f; % grid for foreign bonds

%- - - - - - - - - - - - - - - -  
% set up anonymous functions
%- - - - - - - - - - - - - - - - 

u_prime = @(c) c.^(-sigma); % marginal utility Fxn

%- - - - - - - - - - - - - - - -  
% equilibrium exchange rate
%- - - - - - - - - - - - - - - -  

e1 = (1+i_f)/(1+i_d) * e2;

%- - - - - - - - - - - - - - - -  
% euler equation for domestic bonds b_t+1
%- - - - - - - - - - - - - - - - 

% i) compute integral
fun=@(y2,b_f) u_prime(y2 + e2/p*(1+i_f).*b_f).*lognpdf(y2,mu_y,sd_y); % set up function 
% for which the integral should be computed
LHS = zeros(n_b_f,1); % set up matrix 
parfor i=1:n_b_f
    LHS(i)=integral(@(x) fun(x,b_f(i)),-Inf,Inf); % compute integral
end
LHS = beta*(1+i_d)*LHS;
%LHS = beta*(1+i_f)*e1/e2*LHS;
% ii) compute second term
RHS = (u_prime(y1 - e1./p.*b_f));
%- - - - - - - - - - - - - - - - 
% intersection point of both euler equations
%- - - - - - - - - - - - - - - - 
zero=LHS-RHS;
% Multiple intersections: choose depending on slope
falls=(gradient(zero,b_f(2)-b_f(1))<0);
b_f_opt=interp1(zero(falls),b_f(falls),0,'linear');
%- - - - - - - - - - - - - - - - 

end