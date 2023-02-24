

%{

        MATLAB Code for Assignment 2a
        Submitted by Lauritz Storch (21-607-015)

%}

clear all
close all
clc

% Please enter your working directory 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

pwd = ['/Users/lauritzstorch/MMM Dropbox/Lauritz Storch/Mac/Desktop/' ...
    'International Macro/Assignment 2a - Lauritz Storch'];

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% Please note that some plots are commented out


%% Pre-settings
rng(123); % set seed to ensure reproducibility of results
set(0,'DefaultLineLineWidth',2) % set linewidth to 2 for all plots

% create folders for tables and plots
mkdir '../Assignment 2a - Lauritz Storch' /plots
mkdir '../Assignment 2a - Lauritz Storch' /data

% LaTeX notation in labels
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 4a) Define Exogeneous model parameters

beta = 0.96;    % discount factor
sigma = 2;      % risk aversion coefficient
i_d = 0.02;     % domestic interest rate
i_f = 0.02;     % foreign interest rate
y1 = 3.3;       % income endowment in period t
p = 1;          % price (normalized)

% income and exchange rate parameters for log-normal distribution
mu_y = 1;       % income mean
mu_e = 1;       % exchange mean
var_y = 0.02;   % income variance
var_e = 0.02;   % exchange variance
cov_ye = -0.0001;% covariance of income and echange rate

%% 4b) Set up grids

% Foreign bond grid
b_f_min = 0.05; % lower bound of foreign bond grid
b_f_max = 1;    % upper bound of foreign bond grid
n_b_f = 5000;   % number of grid points (foreign bond)
b_f = linspace(b_f_min, b_f_max, n_b_f)'; % linearly spaced grid (foreign bond)

% Exchange rate grid
e_min = 0;      % lower bound of exchange rate grid
e_max = 5;      % upper bound of exchange rate grid
n_e = 10000;    % number of grid points (exchange rate)
e1 = linspace(e_min, e_max, n_e)'; % linearly spaced grid (exchange rate)

%% 4c) Anonymous functions

u_prime = @(c) c.^(-sigma); % marginal utility

rho = cov_ye/(sqrt(var_y)*sqrt(var_e));
f_ye = @(y, e, var_e, var_y, cov_ye) (exp(-((((log(y)-mu_y)./(sqrt(var_y))).^2-2.*rho.* ...
    ((log(y)-mu_y)./(sqrt(var_y))).*((log(e)-mu_e)./(sqrt(var_e)))+ ...
    ((log(e)-mu_e)./(sqrt(var_e))).^2)./(2.*(1-rho.^2))))./ ...
    (2.*pi.*sqrt(var_y).*sqrt(var_e).*sqrt(1-rho.^2).*y.*e)); % bivariate distribution

%% 4 d) Equilibrium levels of exchange rate and foreign bonds

% save exogenous model parameters and grid as structure variable
parameters = struct('beta', beta, ...
              'sigma', sigma, ...
              'i_d', i_d, ...
              'i_f', i_f, ...
              'y1', y1, ...
              'p', p, ...
              'mu_y', mu_y, ...
              'mu_e', mu_e, ...
              'var_y', var_y, ...
              'var_e', var_e, ...
              'cov_ye', cov_ye, ...
              'b_f', b_f, ...
              'e1', e1);

% Calculate the equilibrium exchange rate and the optimal level of foreign
% bonds
% Takes a while to load ...
[b_f_opt, e1_opt] = SolveEE(parameters,n_b_f);

%% 4 e) Equilibrium Risk Premium

% i)
expect_e2 = integral2((@(y2,e2) e2./e1_opt.* ...
    f_ye(y2,e2, var_e, var_y, cov_ye)),0,Inf,0,Inf);
RP_i = 1+i_d-(1+i_f)*expect_e2;


% ii)
c1 = y1 - 1/p*e1_opt*b_f_opt;
M = @(c2) beta.*c2./u_prime(c1);
expect_M = integral2((@(y2,e2) M(u_prime(y2 + e2./p.*(1+i_f).*b_f_opt)).* ...
    f_ye(y2,e2, var_e, var_y, cov_ye)),0,Inf,0,Inf);
cov_Me = integral2((@(y2,e2) M(u_prime(y2 + e2./p.*(1+i_f).*b_f_opt)).* ...
    e2./e1_opt.* f_ye(y2,e2, var_e, var_y, cov_ye)),0,Inf,0,Inf) - expect_M.*expect_e2;
RP_ii = (1+i_f).*cov_Me./expect_M;


%% 5 Variations in the parameters and its impact on the exchange rate 
% and the level of foreign bonds

% Columns correspond to the subtasks (e.g. column 2 = task 5.ii)

n = 11;
% Takes a while to load (couple of hours) ...
%[b_vec, e_vec] = variable_change(parameters, n, n_b_f, n_e);

%writematrix(b_vec,'/data/foreignbonds.csv') % Save Data as CSV
%writematrix(e_vec,'/data/exchangerate.csv') % Save Data as CSV

% Recommended -------------------------------------------------------------
b_vec = readmatrix('./data/foreignbonds.csv');
e_vec = readmatrix('./data/exchangerate.csv');
RP_vec = RP_5(parameters,n,e_vec);

%{
% i)
[linspace(1.55,2.45,n)',b_vec(:,1),e_vec(:,1),RP_vec(:,1)]'

% ii)
[linspace(0.0001,0.04,n)',b_vec(:,2),e_vec(:,2),RP_vec(:,2)]'

% iii)
[linspace(0.0001,0.04,n)',b_vec(:,3),e_vec(:,3),RP_vec(:,3)]'

% iv)
[linspace(-0.0012,0.001,n)',b_vec(:,4),e_vec(:,4),RP_vec(:,4)]'
%}

%% Plots task 5
% Plot all changes in the exchange rate
figure;
hold on
plot(e_vec(:,1),'LineWidth',3)
plot(e_vec(:,2),'LineWidth',3)
plot(e_vec(:,3),'LineWidth',3)
plot(e_vec(:,4),'LineWidth',3)
plot(6,e1_opt,'ro','LineWidth',3)
xlim([1 n])
xlabel('Index')
ylabel('Equilibrium Exchange Rate')
%title('Change in the equilibrium exchange rate e_t')
legend('$e_t(\sigma)$','$e_t(\sigma_e^2)$','$e_t(\sigma_y^2)$', ...
    '$e_t(\sigma_{ye})$','Initial point', 'FontSize',14)
legend('position',[.7 .15 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changee1_all.eps'], 'epsc')
close all

% Plot changes in the exchange rate due to changes in the risk aversion coeff. (i)
hold on
plot(e_vec(:,1),'LineWidth',3)
plot(6,e1_opt,'ro','LineWidth',3)
xlim([1 n])
set(gca,'XTickLabel',round(linspace(1.55,2.45,11),3));
xlabel('Risk aversion coefficient \sigma')
ylabel('Equilibrium Exchange Rate')
%title('Change in the equilibrium exchange rate e_t')
legend('$e_t(\sigma)$','Initial point','FontSize',14)
legend('position',[.73 .77 .12 .12])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changee1_riskcoeff.eps'], 'epsc')
close all

% Plot changes in the exchange rate due to changes in the variances (ii and iii)
hold on
plot(e_vec(:,2),'LineWidth',3)
plot(e_vec(:,3),'LineWidth',3)
plot(6,e1_opt,'ro','LineWidth',3)
xlim([1 n])
set(gca,'XTickLabel',round(linspace(0.0001,0.04,11),3));
xlabel('Variances \sigma_y^2 and \sigma_e^2')
ylabel('Equilibrium Exchange Rate')
%title('Change in the equilibrium exchange rate e_t')
legend('$e_t(\sigma_e^2)$','$e_t(\sigma_y^2)$','Initial point', ...
    'FontSize',14)
legend('position',[.73 .2 .12 .12])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changee1_var.eps'], 'epsc')
close all

% Plot changes in the exchange rate due to changes in the covariances (iv)
hold on
plot(e_vec(:,4),'LineWidth',3)
plot(6,e1_opt,'ro','LineWidth',3)
xlim([1 n])
set(gca,'XTickLabel',round(linspace(-0.0012,0.001,11),4));
xlabel('Covariance \sigma_{ye}')
ylabel('Equilibrium Exchange Rate')
%title('Change in the equilibrium exchange rate e_t')
legend('$e_t(\sigma_{ye})$','Initial point','FontSize',14)
legend('position',[.73 .77 .12 .12])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changee1_cov.eps'], 'epsc')
close all

% Plot all changes in the level of foreign bonds
figure;
hold on
plot(b_vec(:,1),'LineWidth',3)
plot(b_vec(:,2),'LineWidth',3)
plot(b_vec(:,3),'LineWidth',3)
plot(b_vec(:,4),'LineWidth',3)
plot(6,b_f_opt,'ro','LineWidth',3)
xlim([1 n])
xlabel('Index')
ylabel('Equilibrium Level of Foreign Bonds')
%title('Change in the equilibrium level of foreign bonds b^*_{t+1}')
legend('$b^*_{t+1}(\sigma)$','$b^*_{t+1}(\sigma_e^2)$', ...
    '$b^*_{t+1}(\sigma_y^2)$','$b^*_{t+1}(\sigma_{ye})$', ...
    'Initial point','FontSize',14)
legend('position',[.7 .2 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changebf_all.eps'], 'epsc')
close all

% Plot changes in the exchange rate due to changes in the risk aversion coeff. (i)
hold on
plot(b_vec(:,1),'LineWidth',3)
plot(6,b_f_opt,'ro','LineWidth',3)
xlim([1 n])
set(gca,'XTickLabel',round(linspace(1.55,2.45,11),3));
xlabel('Risk aversion coefficient \sigma')
ylabel('Equilibrium Level of Foreign Bonds')
%title('Change in the equilibrium level of foreign bonds b^*_{t+1}')
legend('$b^*_{t+1}(\sigma)$','Initial point','FontSize',14)
legend('position',[.73 .2 .12 .12])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changebf_riskcoeff.eps'], 'epsc')
close all

% Plot changes in the level of foreign bonds due to changes in the variances (ii and iii)
hold on
plot(b_vec(:,2),'LineWidth',3)
plot(b_vec(:,3),'LineWidth',3)
plot(6,b_f_opt,'ro','LineWidth',3)
xlim([1 n])
set(gca,'XTickLabel',round(linspace(0.0001,0.04,11),3));
xlabel('Variances \sigma_y^2 and \sigma_e^2')
ylabel('Equilibrium Level of Foreign Bonds')
%title('Change in the equilibrium level of foreign bonds b^*_{t+1}')
legend('$b^*_{t+1}(\sigma_e^2)$','$b^*_{t+1}(\sigma_y^2)$', ...
    'Initial point','FontSize',14)
legend('position',[.73 .2 .12 .12])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changebf_var.eps'], 'epsc')
close all

% Plot changes in the level of foreign bonds due to changes in the covariances (iv)
hold on
plot(b_vec(:,4),'LineWidth',3)
plot(6,b_f_opt,'ro','LineWidth',3)
xlim([1 n])
set(gca,'XTickLabel',round(linspace(-0.0012,0.001,11),4));
xlabel('Covariance \sigma_{ye}')
ylabel('Equilibrium Level of Foreign Bonds')
%title('Change in the equilibrium level of foreign bonds b^*_{t+1}')
legend('$b^*_{t+1}(\sigma_{ye})$','Initial point','FontSize',14)
legend('position',[.73 .16 .12 .12])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changebf_cov.eps'], 'epsc')
close all

% Plot all changes in the Risk Premium
figure;
hold on
plot(RP_vec(:,1),'LineWidth',3)
plot(RP_vec(:,2),'LineWidth',3)
plot(RP_vec(:,3),'LineWidth',3)
plot(RP_vec(:,4),'LineWidth',3)
plot(6,RP_i,'ro','LineWidth',3)
xlim([1 n])
xlabel('Index')
ylabel('Equilibrium risk premium')
%title('Change in the equilibrium risk premium RP_t')
legend('$RP_t(\sigma)$','$RP_t(\sigma_e^2)$','$RP_t(\sigma_y^2)$', ...
    '$RP_t(\sigma_{ye})$','Initial point', 'FontSize',14)
legend('position',[.7 .7 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changeRP_all.eps'], 'epsc')
close all


% Plot changes in the risk premium due to changes in the risk aversion coeff. (i)
hold on
plot(RP_vec(:,1),'LineWidth',3)
plot(6,RP_i,'ro','LineWidth',3)
xlim([1 n])
set(gca,'XTickLabel',round(linspace(1.55,2.45,11),3));
xlabel('Risk aversion coefficient \sigma')
ylabel('Equilibrium risk premium')
%title('Change in the equilibrium risk premium RP_t')
legend('$RP_t(\sigma)$','Initial point','FontSize',14)
legend('position',[.73 .77 .12 .12])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changeRP_riskcoeff.eps'], 'epsc')
close all

% Plot changes in the risk premium due to changes in the variances (ii and iii)
hold on
plot(RP_vec(:,2),'LineWidth',3)
plot(RP_vec(:,3),'LineWidth',3)
plot(6,RP_i,'ro','LineWidth',3)
xlim([1 n])
set(gca,'XTickLabel',round(linspace(0.0001,0.04,11),3));
xlabel('Variances \sigma_y^2 and \sigma_e^2')
ylabel('Equilibrium risk premium')
%title('Change in the equilibrium risk premium RP_t')
legend('$RP_t(\sigma_e^2)$','$RP_t(\sigma_y^2)$','Initial point', ...
    'FontSize',14)
legend('position',[.73 .77 .12 .12])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changeRP_var.eps'], 'epsc')
close all

% Plot changes in the risk premium due to changes in the covariances (iv)
hold on
plot(RP_vec(:,4),'LineWidth',3)
plot(6,RP_i,'ro','LineWidth',3)
xlim([1 n])
set(gca,'XTickLabel',round(linspace(-0.0012,0.001,11),4));
xlabel('Covariance \sigma_{ye}')
ylabel('Equilibrium risk premium')
%title('Change in the equilibrium risk premium RP_t')
legend('$RP_t(\sigma_{ye})$','Initial point','FontSize',14)
legend('position',[.73 .77 .12 .12])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changeRP_cov.eps'], 'epsc')
close all

%% 6 Compare results to a model with a perfectly anticipated exchange rate

n = 11;                         % Number of steps in linspace
sigmas = linspace(1.55,2.45,n); % Vector of different sigmas
var_ys = linspace(0.0001,0.04,n);% variance income
e2 = exp(1);                    % exchange rate in period t+1
parameters.e2 = e2;             % Add exchange rate in the second period
parameters.sd_y = sqrt(var_y);  % Add standard deviation of income

% set up matrix to store solution
b_fs_eq6 = zeros(n,4);
e1s_eq6 = zeros(n,4);


% risk aversion coefficient
for i=1:n
    parameters.sigma = sigmas(i);
    [b_fs_eq6(i,1),e1s_eq6(i,1)] = Fx_solveEE(parameters, n_b_f);
end
parameters.sigma = 2;

% standard deviation of income
for i=1:n
    parameters.sd_y = sqrt(var_ys(i));
    [b_fs_eq6(i,3),e1s_eq6(i,3)] = Fx_solveEE(parameters, n_b_f);
end
parameters.sd_y = sqrt(var_y);

% standard deviation future exchange rate and covariance of the exchange rate and the income
[b_fs_eq6(:,[2 4]),e1s_eq6(:,[2 4])] = Fx_solveEE(parameters, n_b_f);
%[b_fs_eq6(:,),e1s_eq6(:,4)] = Fx_solveEE(parameters, n_b_f);

%% Plots task 6
% Plot comparison sigma b_f
figure 
hold on
plot(sigmas,b_vec(:,1),'r','LineWidth',3)
plot(sigmas,b_fs_eq6(:,1),'b','LineWidth',3)
plot(sigmas(6),b_f_opt,'o','LineWidth',3)
plot(sigmas(6),b_fs_eq6(6,1),'o','LineWidth',3)
xlim([sigmas(1) sigmas(n)])
ylabel('Equilibrium Level of Foreign Bonds')
xlabel('Risk aversion coefficient \sigma')
%title('Change in the equilibrium level of foreign bonds b_{t+1}')
legend('uncertain $b^*_{t+1}(\sigma)$','perfectly anticipated $b^*_{t+1}(\sigma)$', ...
    'Initial point (uncertainty model)','Initial point (perfect model)', ...
    'FontSize',14)
legend('position',[.59 .2 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/comparison_b_sigma.eps'], 'epsc')
close all

% Plot comparison sigma e1
figure 
hold on
plot(sigmas,e_vec(:,1),'r','LineWidth',3)
plot(sigmas,e1s_eq6(:,1),'b','LineWidth',3)
plot(sigmas(6),e1_opt,'o','LineWidth',3)
plot(sigmas(6),e1s_eq6(6,1),'o','LineWidth',3)
xlim([sigmas(1) sigmas(n)])
ylabel('Equilibrium Exchange Rate')
xlabel('Risk aversion coefficient \sigma')
%title('Change in the equilibrium exchange rate e_t')
legend('uncertain $e_t(\sigma)$','perfectly anticipated $e_t(\sigma)$', ...
    'Initial point (uncertainty model)','Initial point (perfect model)', ...
    'FontSize',14)
legend('position',[.3 .4 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/comparison_e1_sigma.eps'], 'epsc')
close all

% Plot comparison variance_y bf
figure 
hold on
plot(var_ys,b_vec(:,3),'r','LineWidth',3)
plot(var_ys,b_fs_eq6(:,3),'b','LineWidth',3)
plot(var_ys(6),b_f_opt,'o','LineWidth',3)
plot(var_ys(6),b_fs_eq6(6,2),'o','LineWidth',3)
xlim([0 var_ys(n)])
ylabel('Equilibrium Level of Foreign Bonds')
xlabel('Income variance \sigma_y^2')
%title('Change in the equilibrium level of foreign bonds b_{t+1}')
legend('uncertain $b^*_{t+1}(\sigma^2_y)$','perfectly anticipated $b^*_{t+1}(\sigma^2_y)$', ...
    'Initial point (uncertainty model)','Initial point (perfect model)', ...
    'FontSize',14)
legend('position',[.59 .2 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/comparison_b_vary.eps'], 'epsc')
close all

% Plot comparison variance_y e1
figure 
hold on
plot(var_ys,e_vec(:,3),'r','LineWidth',3)
plot(var_ys,e1s_eq6(:,3),'b','LineWidth',3)
plot(var_ys(6),e1_opt,'o','LineWidth',3)
plot(var_ys(6),e1s_eq6(6,2),'o','LineWidth',3)
xlim([0 var_ys(n)])
ylabel('Equilibrium Exchange Rate')
xlabel('Income variance \sigma_y^2')
%title('Change in the equilibrium exchange rate e_t')
legend('uncertain $e_t(\sigma^2_y)$','perfectly anticipated $e_t(\sigma^2_y)$', ...
    'Initial point (uncertainty model)','Initial point (perfect model)', ...
    'FontSize',14)
legend('position',[.6 .4 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/comparison_e1_vary.eps'], 'epsc')
close all

% Plot comparison variance_e bf
figure 
hold on
plot(var_ys,b_vec(:,2),'r','LineWidth',3)
plot(var_ys,b_fs_eq6(:,2),'b','LineWidth',3)
plot(var_ys(6),b_f_opt,'o','LineWidth',3)
plot(0,b_fs_eq6(6,3),'o','LineWidth',3)
xlim([0 var_ys(n)])
ylabel('Equilibrium Level of Foreign Bonds')
xlabel('Exchange rate variance \sigma_e^2')
%title('Change in the equilibrium level of foreign bonds b_{t+1}')
legend('uncertain $b^*_{t+1}(\sigma^2_e)$','perfectly anticipated $b^*_{t+1}(\sigma^2_e)$', ...
    'Initial point (uncertainty model)','Initial point (perfect model)', ...
    'FontSize',14)
legend('position',[.3 .2 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/comparison_b_vare.eps'], 'epsc')
close all

% Plot comparison variance_e e1
figure 
hold on
plot(var_ys,e_vec(:,2),'r','LineWidth',3)
plot(var_ys,e1s_eq6(:,2),'b','LineWidth',3)
plot(var_ys(6),e1_opt,'o','LineWidth',3)
plot(0,e1s_eq6(6,3),'o','LineWidth',3)
xlim([0 var_ys(n)])
ylabel('Equilibrium Exchange Rate')
xlabel('Exchange rate variance \sigma_e^2')
%title('Change in the equilibrium exchange rate e_t')
legend('uncertain $e_t(\sigma^2_e)$','perfectly anticipated $e_t(\sigma^2_e)$', ...
    'Initial point (uncertainty model)','Initial point (perfect model)', ...
    'FontSize',14)
legend('position',[.3 .7 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/comparison_e1_vare.eps'], 'epsc')
close all

% Plot comparison covariance_ye bf
cov_yes = linspace(-0.0012,0.001,n);
figure 
hold on
plot(cov_yes,b_vec(:,4),'r','LineWidth',3)
plot(cov_yes,b_fs_eq6(:,4),'b','LineWidth',3)
plot(cov_yes(6),b_f_opt,'o','LineWidth',3)
plot(0,b_fs_eq6(6,4),'o','LineWidth',3)
xlim([cov_yes(1) cov_yes(n)])
ylabel('Equilibrium Level of Foreign Bonds')
xlabel('Covariance of the exchange rate and the income \sigma_{ye}')
%title('Change in the equilibrium level of foreign bonds b_{t+1}')
legend('uncertain $b^*_{t+1}(\sigma_{ye})$','perfectly anticipated $b^*_{t+1}(\sigma_{ye})$', ...
    'Initial point (uncertainty model)','Initial point (perfect model)', ...
    'FontSize',14)
legend('position',[.3 .5 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/comparison_b_cov.eps'], 'epsc')
close all

% Plot comparison covariance_ye e1
figure 
hold on
plot(cov_yes,e_vec(:,4),'r','LineWidth',3)
plot(cov_yes,e1s_eq6(:,4),'b','LineWidth',3)
plot(cov_yes(6),e1_opt,'o','LineWidth',3)
plot(0,e1s_eq6(6,3),'o','LineWidth',3)
xlim([cov_yes(1) cov_yes(n)])
ylabel('Equilibrium Exchange Rate')
xlabel('Covariance of the exchange rate and the income \sigma_{ye}')
%title('Change in the equilibrium exchange rate e_t')
legend('uncertain $e_t(\sigma_{ye})$','perfectly anticipated $e_t(\sigma_{ye})$', ...
    'Initial point (uncertainty model)','Initial point (perfect model)', ...
    'FontSize',14)
legend('position',[.3 .4 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/comparison_e1_cov.eps'], 'epsc')
close all

%% 7 Risk premium and variance of the exchange rate approaches zero

% Takes a while to load ...
% RP_i_vec = RP_change(parameters,11,n_b_f);  % Calculate risk premia

% writematrix(RP_i_vec,'./data/riskpremium.csv') % Save Data as CSV

% Recommended -------------------------------------------------------------
RP_i_vec = readmatrix('./data/riskpremium.csv');

%% Plot: Risk Premium with the variance of the exchange rate approaching zero
hold on
plot(RP_i_vec,'LineWidth',3)
plot(RP_i,'ro','LineWidth',3)
plot(linspace(0,0,11),'r--','LineWidth',3)
set(gca,'XTickLabel',round(linspace(0.02,0.001,11),3));
ylim([-0.0045 0.0001])
xlim([1 n])
xlabel('Changes in the variance \sigma^2_{e}')
ylabel('Risk Premium')
%title('Change in the risk premium due to a change in the variance \sigma^2_{e}')
legend('$\sigma^2_{e}$','Initial point','FontSize',14)
legend('position',[.7 .2 .18 .18])
set(gca,'FontSize',18)
hold off

saveas(gcf, [pwd '/plots/changeRP.eps'], 'epsc')
close all
