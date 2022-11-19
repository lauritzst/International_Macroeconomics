

%{

        MATLAB Code for Assignment 1a
        Submitted by Lauritz Storch

%}


% Please enter your working directory 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

pwd = ['/Users/lauritzstorch/MMM Dropbox/Lauritz Storch/' ...
    'Mac/Desktop/International Macro/Assignment 1a'];

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% Please note that some plots are commented out


%% Pre-settings
clear all
close all
clc

rng(123); % set seed to ensure reproduciability of results
set(0,'DefaultLineLineWidth',2) % set linewidth to 2 for all plots

% create folders for tables and plots
mkdir '../Assignment 1a'/tables
mkdir ../'Assignment 1a'/plots

% LaTeX notation in labels
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 4a) Define Exogeneous model parameters


beta = 0.96;    % discount factor
sigma = 2;      % risk aversion coefficient
r = 0.01;       % risk free interest rate
a0 = 0.4332;    % parameter of loss function
a1 = 0.6341;    % parameter of loss function
y1 = 1;         % income endownment of period 1
B1 = 1;         % initial debt
mu_y = .0175;   % exp(mu_y) is mean of income shocks
sd_y = .0685;   % standard deviation of of income shocks


%% 4b) linearly spaced grid for bonds B2

b_min = 0;
b_max = 1.5;
n_b = 1000;
B2 = linspace(b_min, b_max, n_b)';

%% 4c) Anonymous function for marginal utility

u_marg = @(c) c.^(- sigma);   % marginal utility function

%% 4d) Compute the bond price q1 

y_tilde=((B2.*(1+a1))./a0).^(1./(1+a1));  % calculation for threshold value of y2


prob = logncdf(y_tilde,mu_y,sd_y);  % default probability      
q1=(1 - prob)./(1+r).*ones(n_b,1);  % bond price

%% 4e) Compute the elasticity epsilon

% Anonymous function for epsilon depending on the level of debt
epsilon = @(b) (1./(1+r)).*lognpdf(y_tilde,mu_y,sd_y).*(y_tilde).^(-a1) .* (b./(a0.*q1)); 

% elasticity of q1 wrt B2
epsilon_q1B2=epsilon(B2); 

%{
% Plot
plot(epsilon_q1B2)
legend('$\textrm{elasticity } \varepsilon_{q_1,B_2}$')
xlabel('debt levels')
ylabel('elasiticity value')
saveas(gcf, [pwd '/plots/elasitictiy.eps'], 'epsc')
%}

%% 4f) Compute the lhs of the Euler Equation

close all

c1=y1 - B1 + q1.*B2;    % consumption in first period
euler_lhs = u_marg(c1).*q1.*(1 - epsilon_q1B2);     % left hand side (lhs) of the Euler Equation


%% 4g) Compute the rhs of the Euler Equation


% right-hand side (rhs) of the Euler Equation
euler_rhs = zeros(n_b,1);
for i=1:n_b
    euler_rhs(i) = beta.*integral(@(y) (u_marg(y - B2(i)).*lognpdf(y,mu_y,sd_y)), y_tilde(i) ,Inf);
end

%% 4h Compute intersection

% works if we have only one intersection
% look for sign swaps
diff = euler_lhs-euler_rhs
for i=2:length(diff)
    intersect(i) = diff(i).*diff(i-1);
end
index = find(intersect < 0);

% Check sum
% euler_rhs(index) - euler_lhs(index);

% Debt level and Bond price
fprintf(['\n\n\nThe government chooses   %d   values of debt and \n' ...
        'the according bond price is   %d   \n ' ...
        'the optimal probability is   %d   \n'], round(B2(index),5), round(q1(index),5), round(prob(index),5));


%{
% optimal level of debt
plot(B2), hold on
plot(index,B2(index),'r*')
%}


% Intersection plot
plot(euler_lhs), hold on
plot(euler_rhs), hold on
plot(index,euler_rhs(index),'r*','LineWidth',5)
legend ('$\textrm{euler}_{lhs}$', '$\textrm{euler}_{rhs}$', '$\textrm{intersection}$','FontSize',12)
xlabel('Levels of debt')
ylabel('Euler Equation dynamics')
set(gca,'XTick',linspace(0,1000,10));
set(gca,'XTickLabel',round(linspace(0,1.5,10),2));
saveas(gcf, [pwd '/plots/euler_equation.eps'], 'epsc')
hold off


%{
% Due to inaccuracy, one could increase the steps in linespace to 100000 for
% B2 (left out for simplicity):

n_b = 4000;
B2_new = linspace(b_min, b_max, n_b)';
y_tilde=((B2.*(1+a1))./a0).^(1./(1+a1)); 
prob = logncdf(y_tilde,mu_y,sd_y);
q1=(1 - prob)./(1+r).*ones(n_b,1);
epsilon = @(b) (1./(1+r)).*lognpdf(y_tilde,mu_y,sd_y).*(y_tilde).^(-a1) .* (b./(a0.*q1)); 
epsilon_q1B2 = epsilon(B2);
c1 = y1 - B1 + q1.*B2;
euler_lhs = u_marg(c1).*q1.*(1 - epsilon_q1B2);

% calculate lhs of the Euler equation
euler_rhs = zeros(n_b,1);
for i=1:n_b
    euler_rhs(i) = beta.*integral(@(y) (u_marg(y - B2(i)).*lognpdf(y,mu_y,sd_y)), y_tilde(i) ,Inf);
end

% look for sign swaps
diff = euler_lhs-euler_rhs
for i=2:length(diff)
    intersect(i) = diff(i).*diff(i-1);
end
index = find(intersect < 0);

% Check sum
% euler_rhs(index) - euler_lhs(index);

% Debt level and Bond price
fprintf(['\n\n\nThe government chooses   %d   values of debt and \n' ...
        'the according bond price is   %d   \n ' ...
        'the optimal probability is   %d   \n'], round(B2(index),5), round(q1(index),5), round(prob(index),5));

% intersection: index=604
% result: debt = 2.261800e-01, bond price = 9.434900e-01, optimal default
% probability = 4.708000e-02
%}

%% 5 Analyze graphically
close all
clc

% initial debt
% B1 = [1.9,1.6,1.3,1.2,1.1];
B1 = [0.1,0.2,0.3,0.5,0.7];

% income endownment of period 1
%y1 = [0.1,0.4,0.7,0.8,0.9];
y1 = [1.9,1.8,1.7,1.6,1.3];

% risk free interest rate
r = [0.1,0.5,2,10,20]/100;

%Variables = [B1,y1,r];


% -------------------------------------------------------------------
% -------------------------------------------------------------------
j = 1; 
%j = 2; 
%j = 3; 
% Selection index:          i)      j=1 changes in B1 
%                           ii      j=2 changes in y1
%                           iii)    j=3 changes in r
% -------------------------------------------------------------------
% -------------------------------------------------------------------


% Calculations ---
euler_lhs_loop = zeros(n_b,length(B1));
B2_opt_loop = zeros(1,length(B1));
prob_opt_loop = zeros(1,length(B1));
q1_equ_loop = zeros(1,length(B1));

% Loop to recalculate the values accordingly
for i = 1:length(B1)
    if j==1 
        % i) changes in B_1
        y1=1;
        r=0.01;
        c1_loop = y1 - B1(i) + q1.*B2;
        euler_lhs_loop(:,i) = u_marg(c1_loop).*q1.*(1 - epsilon_q1B2);

        diff_loop = euler_lhs_loop(:,i)-euler_rhs;
        for p=2:length(diff_loop)
            intersect_loop(p) = diff_loop(p).*diff_loop(p-1);
        end

        try
            index_loop(i) = find(intersect_loop < 0);
        catch 
            warning('No numerical solution. Graphs do not intersect.');
            index_loop(i) = 1000; % set index equal to previous numerical solution to detect 
                                   % if we have only an approximation but not intersection
        end
        
        B2_opt_loop(1,i) = B2(index_loop(i));
        prob_opt_loop(1,i) = prob(index_loop(i));
        q1_equ_loop(1,i) = q1(index_loop(i));
        
        lgd{i} = strcat('B1 = ',num2str(B1(i)));

        t=length(B1);
        
    elseif j==2
        % ii) changes in y_1
        B1=1;
        r=0.01;
        c1_loop = y1(i) - B1 + q1.*B2;
        euler_lhs_loop(:,i) = u_marg(c1_loop).*q1.*(1 - epsilon_q1B2);

        diff_loop = euler_lhs_loop(:,i)-euler_rhs;
        for p=2:length(diff_loop)
            intersect_loop(p) = diff_loop(p).*diff_loop(p-1);
        end

        try
            index_loop(i) = find(intersect_loop < 0);
        catch 
            warning('No numerical solution. Graphs do not intersect.');
            index_loop(i) = 1000; % set index equal to previous numerical solution to detect 
                                  % if we have only an approximation but not intersection
        end
        
        B2_opt_loop(1,i) = B2(index_loop(i));
        prob_opt_loop(1,i) = prob(index_loop(i));
        q1_equ_loop(1,i) = q1(index_loop(i));
        
        lgd{i} = strcat('y1 = ',num2str(y1(i)));

        t=length(y1);

    else
        % iii) changes in r
        y1=1;
        B1=1;
        q1_loop=(1 - prob)./(1+r(i)).*ones(n_b,1);  % bond price
        epsilon_loop = @(b) (1./(1+r(i))).*lognpdf(y_tilde,mu_y,sd_y).*(y_tilde).^(-a1) .* (b./(a0.*q1_loop)); 
        epsilon_q1B2_loop=epsilon_loop(B2);
        euler_lhs_loop(:,i) = u_marg(c1).*q1_loop.*(1 - epsilon_q1B2_loop);

        diff_loop = euler_lhs_loop(:,i)-euler_rhs;
        for p=2:length(diff_loop)
            intersect_loop(p) = diff_loop(p).*diff_loop(p-1);
        end

        try
            index_loop(i) = find(intersect_loop < 0);
        catch 
            warning('No numerical solution. Graphs do not intersect.');
            index_loop(i) = 1000; % set index equal to previous numerical solution to detect 
                                   % if we have only an approximation but not intersection
        end
        
        B2_opt_loop(1,i) = B2(index_loop(i));
        prob_opt_loop(1,i) = prob(index_loop(i));
        q1_equ_loop(1,i) = q1_loop(index_loop(i));

        % plot(euler_lhs_loop(:,i));
        lgd{i} = strcat('r = ',num2str(r(i)));

        t=length(r);

    end
end


% Plots ---
figure;
hold on;
plot(euler_lhs_loop);
if j == 1
    plot(euler_lhs, 'r--');
    lgd{length(B1)+1} = strcat('$B1 = ',num2str(1), ' (euler_{lhs})$');
    plot(euler_rhs, 'b--');
    lgd{length(B1)+2} = strcat('$euler_{rhs}$');
    legend(lgd);
    legend('position',[.67 .565 .2 .2])
    set(gca,'XTick',linspace(0,1000,10));
    set(gca,'XTickLabel',round(linspace(0,1.5,10),2));
    xlabel('Levels of debt (B_2)')
    ylabel('LHS of the Euler Equation')
elseif j == 2
    plot(euler_lhs, 'r--');
    lgd{length(y1)+1} = strcat('$y_1 = ',num2str(1), ' (euler_{lhs})$');
    plot(euler_rhs, 'b--');
    lgd{length(y1)+2} = strcat('$euler_{rhs}$');
    legend(lgd);
    legend('position',[.67 .565 .2 .2])
    set(gca,'XTick',linspace(0,1000,10));
    set(gca,'XTickLabel',round(linspace(0,1.5,10),2));
    xlabel('Levels of debt')
    ylabel('LHS of the Euler Equation')
else
    plot(euler_lhs, 'r--');
    lgd{length(r)+1} = strcat('$r = ',num2str(0.01), ' (euler_{lhs})$');
    plot(euler_rhs, 'b--');
    lgd{length(r)+2} = strcat('$euler_{rhs}$');
    legend(lgd);
    legend('position',[.67 .565 .2 .2])
    set(gca,'XTick',linspace(0,1000,10));
    set(gca,'XTickLabel',round(linspace(0,1.5,10),2));
    xlabel('Levels of debt B_2')
    ylabel('LHS of the Euler Equation')
end

% Zoom in figure ---
axes('position',[.53 .175 .34 .34]);
% put box around new pair of axes
box on 
plot(euler_lhs_loop(:,1:5)); hold on,
plot(euler_lhs, 'r--');
% zoom into intersection area 
ax = gca;
if j==1 || j==2
    ax.XLim = [90 210];
    ax.YLim = [-10 10];
elseif j==3
    ax.XLim = [151.55 151.75];
    ax.YLim = [1.43 1.47];
end
% add right-hand side of the euler equation
plot(euler_rhs, 'b--');
lgd{length(B1)+2} = strcat('$euler_{rhs}$');
% Add title
title('Zoomed in'); hold off,

% Table ---
% previous values
entry_nochange = [B2(index);q1(index);prob(index)];
var_name  = {}; 
entry = [];
% loop to fill newly calculated values into table format
for i=1:t
    if j == 1
        var_name = cat(2,var_name, strcat('B1 = ', num2str(B1(i))));
        no_change = 'B1=1';
    elseif j == 2
        var_name = cat(2,var_name, strcat('y1 = ', num2str(y1(i))));
        no_change='y1=1';
    else 
        var_name = cat(2,var_name, strcat('r = ', num2str(r(i))));
        no_change='r=0.01';
    end    
        
    eval([strcat('entry', num2str(i)), ' = [B2_opt_loop(1,i);q1_equ_loop(1,i);prob_opt_loop(1,i)]']);
    clc
end
var_name = cat(2,var_name,no_change);

% create table
T = table(entry1,entry2,entry3,entry4,entry5, ...
    entry_nochange,'VariableNames',var_name,'RowName',{'B2','q1','prob.(default)'}); 
disp(T)

% Save tables and plots
if j==1
    saveas(gcf, [pwd '/plots/sensitivity_B1.eps'], 'epsc')
    table2latex(T,[pwd '/tables/sensitivity_B1'])
elseif j==2
    saveas(gcf, [pwd '/plots/sensitivity_y1.eps'], 'epsc')
    table2latex(T,[pwd '/tables/sensitivity_y1'])
else
    saveas(gcf, [pwd '/plots/sensitivity_r.eps'], 'epsc')
    table2latex(T,[pwd '/tables/sensitivity_r'])
end


%% 6 Modified loss function

B1 = 1;         % initial debt
y1 = 1;         % income endownment of period 1
r = 0.01;       % risk free interest rate

% new loss parameter
alpha=0.8;

% euler_lhs derivation ---
% y_tilde
y_tilde_mod = B2./alpha;
%(alpha).^(-1).*B2;
% probability
prob_mod = logncdf(y_tilde_mod,mu_y,sd_y);
% bond price
q1_mod=(1 - prob_mod)./(1+r).*ones(n_b,1);
% elasticity
epsilon_mod = @(b) (1./(1+r)).*lognpdf(y_tilde_mod,mu_y,sd_y).*b./(alpha.*q1_mod);
epsilon_q1B2_mod = epsilon_mod(B2);
% euler_lhs
c1_mod = y1 - B1 + q1_mod.*B2;
euler_lhs_mod = u_marg(c1_mod).*q1_mod.*(1 - epsilon_q1B2_mod);


% euler_rhs derivation ---
euler_rhs_mod = zeros(n_b,1);
for i=1:n_b
    euler_rhs_mod(i) = beta.*integral(@(y) (u_marg(y - B2(i)).*lognpdf(y,mu_y,sd_y)), (y_tilde_mod(i)) ,Inf);
end

% new calculation for the intersection
intersect_mod = min(abs(euler_lhs_mod - euler_rhs_mod));
index_mod = find(intersect_mod == abs(euler_lhs_mod - euler_rhs_mod));


%{
% Check sum
euler_rhs_mod(index_mod) - euler_lhs_mod(index_mod)
%}

% old penalization values
fprintf(['\n\n\nThe government chooses   %d   values of debt and \n' ...
        'the according bond price is   %d   \n'...
        'the optimal probability is   %d   \n if L(y_2) = a_0 * (y_2^(1+a1))/(1+a1) \n'], ...
        round(B2(index),5), round(q1(index),5), round(prob(index),5));

% new penalization values
fprintf(['\n\n\nThe government chooses   %d   values of debt and \n' ...
        'the according bond price is   %d   \n' ...
        'the optimal probability is   %d   \n if L(y_2) = alpha*y_2  \n'], ...
        round(B2(index_mod),5), round(q1_mod(index_mod),5), round(prob_mod(index_mod),5));

%% 6. levels of debt (plot)
close all

plot(B2), hold on
plot(index,B2(index),'r*','LineWidth',5)
plot(index_mod,B2(index_mod),'b*','LineWidth',5)
legend ('$B_2$','$L(y_2) = a_0  \frac{y_2^{1+a1}}{1+a1}$', '$L(y_2) = \alpha y_2$')
xlabel('index')
ylabel('Amount of debt B_2')
legend('position',[.69 .72 .15 .15],'FontSize',14)
saveas(gcf, [pwd '/plots/B2_comparison.eps'], 'epsc')
hold off

%% 6. bond price (plot)
close all

hold on
plot(q1)
plot(q1_mod)
plot(index,q1(index),'r*','LineWidth',5)
plot(index_mod,q1_mod(index_mod),'b*','LineWidth',5)
ylabel('Bond price q_1(B_2)')
legend ('$q_{initial}$', '$q_{modified}$','$L(y_2) = a_0  \frac{y_2^{1+a1}}{1+a1}$', '$L(y_2) = \alpha y_2$')
legend('position',[.69 .72 .15 .15],'FontSize',14)
set(gca,'XTick',linspace(0,1000,10));
set(gca,'XTickLabel',round(linspace(0,1.5,10),2));
xlabel('Levels of debt B_2')
saveas(gcf, [pwd '/plots/q1_comparison.eps'], 'epsc')
hold off

%% 6. default probability (plot)
close all

hold on
plot(prob)
plot(prob_mod)
plot(index,prob(index),'r*','LineWidth',5)
plot(index_mod,prob_mod(index_mod),'b*','LineWidth',5)
ylabel('Default probability \delta')
legend ('$\delta_{initial}$', '$\delta_{modified}$','$L(y_2) = a_0  \frac{y_2^{1+a1}}{1+a1}$', '$L(y_2) = \alpha y_2$')
legend('position',[.69 .72 .15 .15],'FontSize',14)
saveas(gcf, [pwd '/plots/delta_comparison.eps'], 'epsc')
set(gca,'XTick',linspace(0,1000,10));
set(gca,'XTickLabel',round(linspace(0,1.5,10),2));
xlabel('Levels of debt B_2')
hold off

%% 7 Mulitple Equilibria

% written in word



