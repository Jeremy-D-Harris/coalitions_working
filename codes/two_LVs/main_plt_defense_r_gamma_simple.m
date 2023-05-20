% void = main_get_ss_LVmodel_022823(void)
% steady-state analysis of single LV-model

%%
clear all; close all; clc;

%% save figure?
save_ans_Fig = 1;
% 0: don't save
% 1: save

figure_name= 'LV_defense_r_gamma_R0equalsone';

my_rgb_colors = [78 132 193; 209 109 106; 236 180 118; 127 127 127]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];

%% system parameters (units of hours, micrograms and mL).
r_fixed = 0.9995;
gam_fixed = 1e-6;

r_pt = 1;
gam_pt = 0;

% r_pt_line = 1.0005;
% gam_pt_line = 1e-4

%% initial conditions
VB0 = 10^0;
LB0 = 10^0;

%% system parameters (units of hours, micrograms and mL).
% Assumes the system is a 500 mL flask running for ~ 24hr;
% conversion_efficiency = 5e-7; %ug/cell
% d_R = 0; % per hour
% mu_max = 1.2; % growth rate (per hour)
r_A = 1; % growth rate (per hour)
gam_A = 1e-4; % lysis rate (per hour)

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



%% range over r and gamma
updown_amt = 0.1;
r_vals_range = linspace(r_pt-updown_amt,r_pt+updown_amt,100);
gam_vals_range = 10.^linspace(-6,-2,200);

% r_fixed_vector = r_fixed*ones(size(gam_vals_range));
r_corresponding_tradeoff = zeros(size(gam_vals_range));
% r_corresponding_tradeoff_pt = zeros(size(gam_vals_range));

% fixed equilibrium: from which to invade
lysogen_equilibrium_fixed = (K*(r_fixed-gam_fixed-d)/r_fixed);
% lysogen_equilibrium_fixed_pt = (K*(r_pt-gam_pt-d)/r_pt);
% lysogen_equilibrium_fixed_pt_drawline = (K*(r_pt_line-gam_pt_line-d)/r_pt_line);

% lysogen_equilibrium_fixed_nophage = (K*(r_fixed-d)/r_fixed);
phage_equilibrium_fixed = bet*gam_fixed*lysogen_equilibrium_fixed/(phi*lysogen_equilibrium_fixed+m);

% for count_r = 1:length(r_vals_range)
%
%     this_r = r_vals_range(count_r);

for count_gam = 1:length(gam_vals_range)
    
    this_gam = gam_vals_range(count_gam);
    
    this_phage_equilibrium = bet*this_gam*lysogen_equilibrium_fixed/(phi*lysogen_equilibrium_fixed+m);
    
    r_corresponding_lysogenicR0equalsone(1,count_gam) =  (d+phi*this_phage_equilibrium+this_gam)/(1-lysogen_equilibrium_fixed/K);
    
    
end




%% plot R0 wrt r vs. gamma
f1 = figure(1); set(f1, 'Position', [100 800 600 450]);
semilogx(gam_vals_range,r_corresponding_lysogenicR0equalsone,'k--','linewidth',4); hold on;
semilogx(10^-3,1.01,'.','Color',my_rgb_colors(4,:),'MarkerSize',80); hold on;
semilogx(10^-5,1.01,'.','Color',my_rgb_colors(1,:),'MarkerSize',80); hold on;


xlim([1e-6 1e-2]);
ylim([r_pt-0.01 r_pt+.03]);
yticks([0.99 1 1.01 1.02 1.03]);
% % axis([x_values(1) x_values(end) 0 2*10^5]);
% xlabel('(Resident) Lysogen Induction Rate, $\gamma_A$ (hr$^{-1}$)','interpreter','latex');
% ylabel('(Invading) Lysogen Growth Rate, $r_B$ (hr$^{-1}$)','interpreter','latex');
% title('Lysogenic Reproduction Number');
f2=gca;
% f2.CLim = [1 max(max(lysogen_reproduction_number_lysogenic))];
% f2.XScale = 'log';
% f1.YScale = 'log';
f2.LineWidth = 1.5;
f2.FontSize = 28;
f2.FontWeight = 'bold';
f2.FontName = 'Times New Roman';
% f2.ColorScale = 'linear';
% f2.ColorScale = 'linear';


if save_ans_Fig
    
    folder_location = './figures/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n'); 
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n'); 
    fprintf(strcat(folder_location,'\n\n'));
    
    else
    
    fprintf('Figure not saved.\n');
    
end

