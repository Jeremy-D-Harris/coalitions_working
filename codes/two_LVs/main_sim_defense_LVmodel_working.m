% void = main_get_ss_LVmodel_022823(void)
% steady-state analysis of single LV-model

%%
clear all; close all; clc;

%% want to save?
save_ans = 1;
% 0: don't save
% 1: save


% filename = 'LV_defense_r_gamma_LAinductionHigh.mat';
filename = 'LV_defense_r_gamma_LAinductionLow.mat';


my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];

%% system parameters (units of hours, micrograms and mL).
% resident
r_fixed = 0.9995;
gam_fixed = 1e-5;
% invader
r_pt = 1.01;
gam_pt = 1e-3;

%% initial conditions
% VA0 = 10^0;
% LA0 = 10^0;
VB0 = 10^0;
LB0 = 10^0;


%% system parameters (units of hours, micrograms and mL).
% Assumes the system is a 500 mL flask running for ~ 24hr;
% conversion_efficiency = 5e-7; %ug/cell
% d_R = 0; % per hour
% mu_max = 1.2; % growth rate (per hour)
r_A = r_fixed; % growth rate (per hour)
gam_A = gam_fixed; % lysis rate (per hour)

r_B = r_pt; % growth rate (per hour)
gam_B = gam_pt; % lysis rate (per hour)

K = 2e8;
% R_in = 5; %ug/mL
% d_S = .2; % death rate susceptibles (per hour)
% d_E = .2; % death rate exposed (per hour)
d = 0.2; % death rate lysogens (per hour)
% d_I = .2; % death rate infected (per hour)
% lam = 2; % commitment rate (per hour)
bet = 5; % burst size
phi = 3.4e-10; %3.4e-10; % adsorption rate (mL/hr)
m = 1/24; % virus washout (per hour)
% alpha_s = 0; % selection coefficient: alpha_s>1 corresponds to advantage of lysogen over susceptible
% J = 0; %ug/mL-h

% rng(1);

%simulation parameters:
dt = 1; % hours
t_end = 200; % hours
t_span = transpose(0:dt:t_end); % time
% NRuns = 100;
% p = linspace(0,1,6); %dilution factor

% R0 = [10 logspace(2,8,6)]; %initial resource amount in ug/mL ( 500 mL flask)

% flask_volume = 500; %volume in mL


% set up parameter structure
% parameters
% params.conversion_efficiency = conversion_efficiency;
% params.d_R = d_R;
params.r_A = r_A;
params.gam_A = gam_A;
params.r_B = r_B;
params.gam_B = gam_B;
params.K = K;
% params.R_in = R_in;
% params.d_S = d_S;
% params.d_E = d_E;
params.d = d;
% params.d_I = d_I;
% params.lambda = lam;
params.bet = bet;
params.phi = phi;
params.m = m;
% params.alpha_s = alpha_s;
% params.J = J;
params.dt = dt;
params.t_span = t_span;
params.t_end = t_end;
% params.flask_volume = flask_volume; %%flask volume in mL
% params.L0 = L0;
% params.V0 = V0;




%% Simulate model
% options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events', @myEvent);
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% Opt    = odeset();
% [T, Y] = ode45(@YourFun, T, Y0, Opt);

% L0,V0 - 2 initial conditions
% black (10^-4,1.005), yellow-orange
Lysogen_equilibrium_nonzero = (K*(r_A-gam_A-d)/r_A);
phage_equilibrium_nonzero = bet*gam_A*Lysogen_equilibrium_nonzero/(phi*Lysogen_equilibrium_nonzero+m);
LA0 = Lysogen_equilibrium_nonzero;
VA0 = phage_equilibrium_nonzero;

init_conds = [LA0;VA0;LB0;VB0];

% 0, 10^-4
% gam_vals_range = [gam_pt,10^-4];

% init_fraction_LB = LB0/(LA0+LB0);
% init_fraction_VB = VB0/(VA0+VB0);
% init_conds_range = [2e4, 0.9e6;200,200];

% for count = 1:length(gam_vals_range)

%     this_init_conds = init_conds(:,count);

% this_L_traj = zeros(1,length(t_span));
% this_V_traj = zeros(1,length(t_span));
%     this_gam_A = gam_vals_range(count);
%     params.gam_A = this_gam_A;

[t_traj,y_traj] = ode45(@(t,y)simulate_2LVmodel(t,y,params), params.t_span, init_conds, options);

LA_traj = y_traj(:,1)';
VA_traj = y_traj(:,2)';
LB_traj = y_traj(:,3)';
VB_traj = y_traj(:,4)';


%% plot R0 wrt r vs. gamma
f1 = figure(1); set(f1, 'Position', [100 800 600 450]);
g(1) = semilogy(t_span,LA_traj,'k','linewidth',4); hold on;
g(2) = semilogy(t_span,VA_traj,'k--','linewidth',4); hold on;
g(3) = semilogy(t_span,VB_traj,'--','Color',my_rgb_colors(1,:),'linewidth',4); hold on;
g(4) = semilogy(t_span,LB_traj,'Color',my_rgb_colors(1,:),'linewidth',4); hold on;
xlim([0 t_end]);
ylim([10^0 10^10]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Time (hr)');
ylabel('Population Density');
% title('Basic Reproduction Number - lysogenic path, $\mathcal R_0$','interpreter','latex');
f2=gca;
% f2.CLim = [1 max(max(lysogen_reproduction_number_lysogenic))];
% f2.XScale = 'log';
% f1.YScale = 'log';
f2.LineWidth = 1;
f2.FontSize = 24;
f2.FontWeight = 'normal';
f2.FontName = 'Times New Roman';
f2.ColorScale = 'linear';

% legend_char1 = ['Resident Lysogen, $L_A$ ($\gamma_A > 0$)'];
% legend_char2 = ['Invading Lysogen, $L_B$ ($\gamma_A > 0$)'];
% legend_char3 = ['Resident Lysogen, $L_A$ ($\gamma_A = 0$)'];
% legend_char4 = ['Invading Lysogen, $L_B$ ($\gamma_A = 0$)'];

% legend(lg,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Location','SouthWest','FontSize',13);

%% collect results
results.y_traj = y_traj; 

%% save simulated data
if save_ans==1
    
    folder_location = './sim_data/';
    save(strcat(folder_location,filename),'params','results');
    
    fprintf('Saved to file: \n');
    fprintf(strcat(filename,'\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('data not saved.\n');
    
end

