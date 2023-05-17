% void = main_single_LVmodel_022823(void)
% simulate single LV-model

%%
clear all; close all; clc;

%% want to save?
save_ans = 0;
% 0: don't save
% 1: save

filename = 'LV_phaseplane_overtime.mat';
% figure_name = 'LV_phaseplane_notraj.eps';

my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];

%% system parameters (units of hours, micrograms and mL).
% Assumes the system is a 500 mL flask running for ~ 24hr;
% conversion_efficiency = 5e-7; %ug/cell
% d_R = 0; % per hour
% mu_max = 1.2; % growth rate (per hour)
r = 1; % growth rate (per hour)
K = 2e8;
% R_in = 5; %ug/mL
% d_S = .2; % death rate susceptibles (per hour)
% d_E = .2; % death rate exposed (per hour)
d = 0.2; % death rate lysogens (per hour)
% d_I = .2; % death rate infected (per hour)
% lam = 2; % commitment rate (per hour)
gam = 1e-4; % lysis rate (per hour)
bet = 10; % burst size
phi = 3.4e-10; %3.4e-10; % adsorption rate (mL/hr)
m = 1/24; % virus washout (per hour)
% alpha_s = 0; % selection coefficient: alpha_s>1 corresponds to advantage of lysogen over susceptible
% J = 0; %ug/mL-h

% rng(1);

%simulation parameters:
dt = 0.05; % hours
t_end = 40; % hours
t_span = transpose(0:dt:t_end); % time
% NRuns = 100;
% p = linspace(0,1,6); %dilution factor

% R0 = [10 logspace(2,8,6)]; %initial resource amount in ug/mL ( 500 mL flask)

% flask_volume = 500; %volume in mL


% set up parameter structure
% parameters
% params.conversion_efficiency = conversion_efficiency;
% params.d_R = d_R;
params.r = r;
params.K = K;
% params.R_in = R_in;
% params.d_S = d_S;
% params.d_E = d_E;
params.d = d;
% params.d_I = d_I;
% params.lambda = lam;
params.gam = gam;
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
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% L0 = 1e4;
% V0 = 200;
% L0,V0 - 4 initial conditions
% black, green, light blue, maroon
V0 = [5e6,1e6,3e4,3e4];
L0 = [5e8,1e6,5e4,2e9];
init_conds_range = [L0;V0];
% init_conds_range = [2e4, 0.9e6;200,200];

for count = 1:length(init_conds_range)
    
    this_init_conds = init_conds_range(:,count);
    
    % this_L_traj = zeros(1,length(t_span));
    % this_V_traj = zeros(1,length(t_span));
    
    [t,y_traj] = ode45(@(t,y)simulate_LVmodel(t,y,params), params.t_span, this_init_conds, options);
    
    this_L_traj = y_traj(:,1)';
    this_V_traj = y_traj(:,2)';
    
    L_traj(count,:) = this_L_traj;
    V_traj(count,:) = this_V_traj;
    
end

%% nullclines, equilibria, asymptotes
x_values = linspace(10^3,8e7,600);%10.^linspace(1,5,100);
y_values = linspace(10^4,10^10,100);%10.^linspace(1,5,100);

Lysogen_equilibrium_zero = zeros(size(x_values));
Lysogen_equilibrium_nonzero = (K*(r-gam-d)/r);
Lysogen_equilibrium_nonzero_vector = Lysogen_equilibrium_nonzero*ones(size(x_values));
phage_nullcline = m*x_values./(bet*gam-phi*x_values);%bet*eta*x_values./(phi*x_values+m);
phage_equilibrium_nonzero = bet*gam*Lysogen_equilibrium_nonzero/(phi*Lysogen_equilibrium_nonzero+m);
phage_nullcline_asymptote = bet*gam/phi;


%% plot L-V phase plane
f1 = figure(1); set(f1, 'Position', [200 500 1200 450]);
subplot(1,2,1);
% semilogy(x_values,Lysogen_equilibrium_zero,'Color',my_rgb_colors(1,:)); hold on;
h(1)=loglog(x_values,Lysogen_equilibrium_nonzero_vector,'--','Color',my_rgb_colors(2,:),'linewidth',1); hold on;
h(2)=loglog(x_values,phage_nullcline,'--','Color',my_rgb_colors(1,:),'linewidth',1); hold on;
% loglog(phage_nullcline_asymptote*ones(size(y_values)),y_values,'k--','linewidth',1); hold on;
loglog(V_traj(1,:),L_traj(1,:),'k','linewidth',2); hold on;
loglog(V_traj(2,:),L_traj(2,:),'Color',default_rgb_colors(4,:),'linewidth',2); hold on;
loglog(V_traj(3,:),L_traj(3,:),'Color',default_rgb_colors(3,:),'linewidth',2); hold on;
loglog(V_traj(4,:),L_traj(4,:),'Color',default_rgb_colors(2,:),'linewidth',2); hold on;
h(3)=loglog(init_conds_range(2,:),init_conds_range(1,:),'ko','MarkerSize',9,'linewidth',2); hold on;
h(4)=loglog(phage_equilibrium_nonzero,Lysogen_equilibrium_nonzero,'k.','MarkerSize',30); hold on;

xlim([10^4 1e8]);
ylim([10^3 10^10]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Virus Population, $V$ (PFU/mL)','interpreter','latex');
ylabel('Lysogen Population, $L$ (CFU/mL)','interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 20;
f1.FontName = 'Times New Roman';
f1.FontWeight = 'light';


% include legend
legend_char1 = ['$L$-nullcline'];
legend_char2 = ['$V$-nullcline'];
legend_char3 = ['Initial condition'];
legend_char4 = ['$L$-$V$-equilibirum'];
legend(h,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Location','SouthEast','FontSize',16);

%% plot L,V over time
subplot(1,2,2);
g(1)=semilogy(t_span,L_traj(1,:),'Color',default_rgb_colors(1,:),'linewidth',2); hold on;
g(2)=semilogy(t_span,V_traj(1,:),'--','Color',default_rgb_colors(1,:),'linewidth',2); hold on;
semilogy(t_span,L_traj(2,:),'Color',default_rgb_colors(4,:),'linewidth',2); hold on;
semilogy(t_span,V_traj(2,:),'--','Color',default_rgb_colors(4,:),'linewidth',2); hold on;
semilogy(t_span,L_traj(3,:),'Color',default_rgb_colors(3,:),'linewidth',2); hold on;
semilogy(t_span,V_traj(3,:),'--','Color',default_rgb_colors(3,:),'linewidth',2); hold on;
semilogy(t_span,L_traj(4,:),'Color',default_rgb_colors(2,:),'linewidth',2); hold on;
semilogy(t_span,V_traj(4,:),'--','Color',default_rgb_colors(2,:),'linewidth',2); hold on;

xlim([0 t_end]);
ylim([10^3 10^10]);

xlabel('Time (hours)');
ylabel('Population Density');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 20;
f1.FontName = 'Times New Roman';
f1.FontWeight = 'light';



% include legend
legend_char1_g = ['Lysogen population'];
legend_char2_g = ['Virus population'];
legend(g,{legend_char1_g,legend_char2_g}, 'Interpreter','Latex','Location','NorthEast','FontSize',16);



%% collect results
results.Lysogen_equilibrium_nonzero_vector = Lysogen_equilibrium_nonzero_vector;
results.phage_nullcline = phage_nullcline;
results.phage_equilibrium_nonzero = phage_equilibrium_nonzero;
results.Lysogen_equilibrium_nonzero = Lysogen_equilibrium_nonzero;
results.L_traj = L_traj;
results.V_traj = V_traj;
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

