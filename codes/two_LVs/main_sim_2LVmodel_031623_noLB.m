% void = main_single_LVmodel_022823(void)
% simulate single LV-model

%%
clear all; close all; clc;

%% want to save?
save_ans = 0;
% 0: don't save
% 1: save

filename = '2LV_phaseplane_overtime.mat';
% figure_name = 'LV_phaseplane_notraj.eps';

my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];

%% initial conditions
VB0 = 10^0;
LB0 = 10^0+0.01;

%% system parameters (units of hours, micrograms and mL).
% Assumes the system is a 500 mL flask running for ~ 24hr;
% conversion_efficiency = 5e-7; %ug/cell
% d_R = 0; % per hour
% mu_max = 1.2; % growth rate (per hour)
r_A = 1; % growth rate (per hour)
gam_A = 1e-4; % lysis rate (per hour)
r_B = 1; % growth rate (per hour)
gam_B = 0; % lysis rate (per hour)

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
t_end = 240; % hours
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



% VA0 = 1e6;
% LA0 = 1e8;

% L0,V0 - 4 initial conditions
% black, green, light blue, maroon
% V0 = [5e6,1e6,3e4,3e4];
% L0 = [5e8,1e6,5e4,2e9];
Lysogen_equilibrium_nonzero = (K*(r_A-gam_A-d)/r_A);
phage_equilibrium_nonzero = bet*gam_A*Lysogen_equilibrium_nonzero/(phi*Lysogen_equilibrium_nonzero+m);
LA0 = Lysogen_equilibrium_nonzero;
VA0 = phage_equilibrium_nonzero;



init_conds_range = [LA0;VA0;LB0;VB0];

init_fraction_LB = LB0/(LA0+LB0);
init_fraction_VB = VB0/(VA0+VB0)
% init_conds_range = [2e4, 0.9e6;200,200];

for count = 1%length(init_conds_range)
    
    this_init_conds = init_conds_range(:,count);
    
    % this_L_traj = zeros(1,length(t_span));
    % this_V_traj = zeros(1,length(t_span));
    
    [t_traj,y_traj] = ode45(@(t,y)simulate_2LVmodel(t,y,params), params.t_span, this_init_conds, options);
    
    this_LA_traj = y_traj(:,1)';
    this_VA_traj = y_traj(:,2)';
    this_LB_traj = y_traj(:,3)';
    this_VB_traj = y_traj(:,4)';
    
    this_fraction_LB = this_LB_traj./(this_LA_traj+this_LB_traj);
    this_fraction_VB = this_VB_traj./(this_VA_traj+this_VB_traj);
    
    LA_traj(count,:) = this_LA_traj;
    VA_traj(count,:) = this_VA_traj;
    LB_traj(count,:) = this_LB_traj;
    VB_traj(count,:) = this_VB_traj;
    
    fraction_LB(count,:) = this_fraction_LB;
    fraction_VB(count,:) = this_fraction_VB;
    
end



%% nullclines, equilibria, asymptotes
x_values = linspace(10^0,3e7,600);%10.^linspace(1,5,100);
x_values_pnull = linspace(10^0,1e10,100);%10.^linspace(1,5,100);
y_values = linspace(10^0,10^9,100);%10.^linspace(1,5,100);

Lysogen_equilibrium_zero = zeros(size(x_values));
% Lysogen_equilibrium_nonzero = (K*(r_A-gam_A-d)/r_A);
Lysogen_equilibrium_nonzero_vector = Lysogen_equilibrium_nonzero*ones(size(x_values_pnull));
phage_nullcline = m*x_values./(bet*gam_A-phi*x_values);%bet*eta*x_values./(phi*x_values+m);
% phage_equilibrium_nonzero = bet*gam_A*Lysogen_equilibrium_nonzero/(phi*Lysogen_equilibrium_nonzero+m);
phage_nullcline_asymptote = bet*gam_A/phi;



%% plot L,V over time
f1 = figure(1); set(f1, 'Position', [200 500 600 450]);
g(1)=semilogy(t_traj,LA_traj(1,:),'Color',default_rgb_colors(1,:),'linewidth',5); hold on;
g(2)=semilogy(t_traj,VA_traj(1,:),'--','Color',default_rgb_colors(1,:),'linewidth',5); hold on;
% g(3)=semilogy(t_traj,LB_traj(1,:),'Color',my_rgb_colors(1,:),'linewidth',2); hold on;
g(4)=semilogy(t_traj,VB_traj(1,:),'--','Color',[0.7 0.7 0.7],'linewidth',5); hold on;


% g(1)=semilogy(t_traj,LA_traj(1,:),'Color',default_rgb_colors(1,:),'linewidth',5); hold on;
% g(2)=semilogy(t_traj,VA_traj(1,:),'--','Color',default_rgb_colors(1,:),'linewidth',5); hold on;
% g(3)=semilogy(t_traj,LB_traj(1,:),'Color',my_rgb_colors(1,:),'linewidth',5); hold on;
% g(4)=semilogy(t_traj,VB_traj(1,:),'--','Color',my_rgb_colors(1,:),'linewidth',5); hold on;

xlim([0 t_end]);
ylim([10^0 10^10]);

xlabel('Time (hr)');
ylabel('Population Density');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 25;
f1.FontName = 'Times New Roman';
f1.FontWeight = 'bold';



% include legend
legend_char1_g = ['Resident Lysogen'];
legend_char2_g = ['Resident Virus'];
legend_char3_g = ['Invading Lysogen'];
legend_char4_g = ['Invading Virus'];
% legend(g,{legend_char1_g,legend_char2_g,legend_char3_g,legend_char4_g}, 'Interpreter','Latex','Location','SouthEast','FontSize',16);



%% collect results
results.Lysogen_equilibrium_nonzero_vector = Lysogen_equilibrium_nonzero_vector;
results.phage_nullcline = phage_nullcline;
results.phage_equilibrium_nonzero = phage_equilibrium_nonzero;
results.Lysogen_equilibrium_nonzero = Lysogen_equilibrium_nonzero;
results.L_traj = LA_traj;
results.V_traj = VA_traj;
results.x_values = x_values;
results.y_values = y_values;
results.t_span = t_span;
results.init_conds_range = init_conds_range;


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

