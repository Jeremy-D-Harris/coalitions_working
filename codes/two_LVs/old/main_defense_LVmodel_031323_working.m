% void = main_get_ss_LVmodel_022823(void)
% steady-state analysis of single LV-model

%%
clear all; close all; clc;

%% want to save?
save_ans = 0;
% 0: don't save
% 1: save

filename = 'LV_defense_r_gamma.mat';

my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];

%% system parameters (units of hours, micrograms and mL).
r_fixed = 0.9995;
gam_fixed = 0;

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

%% simulate gamma = 10^-4, r = 1.0005

if 1
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
gam_vals_range = [gam_pt,10^-4];

% init_fraction_LB = LB0/(LA0+LB0);
% init_fraction_VB = VB0/(VA0+VB0);
% init_conds_range = [2e4, 0.9e6;200,200];

for count = 1:length(gam_vals_range)
    
%     this_init_conds = init_conds(:,count);
    
    % this_L_traj = zeros(1,length(t_span));
    % this_V_traj = zeros(1,length(t_span));
    this_gam_A = gam_vals_range(count);
    params.gam_A = this_gam_A;
    
    [t_traj,y_traj] = ode45(@(t,y)simulate_2LVmodel(t,y,params), params.t_span, init_conds, options);
    
    this_LA_traj = y_traj(:,1)';
    this_VA_traj = y_traj(:,2)';
    this_LB_traj = y_traj(:,3)';
    this_VB_traj = y_traj(:,4)';
    
    LA_traj(count,:) = this_LA_traj;
    VA_traj(count,:) = this_VA_traj;
    LB_traj(count,:) = this_LB_traj;
    VB_traj(count,:) = this_VB_traj;
    
    this_fraction_LB = this_LB_traj./(this_LA_traj+this_LB_traj);
    this_fraction_VB = this_VB_traj./(this_VA_traj+this_VB_traj);
    
    fraction_LB(count,:) = this_fraction_LB;
    fraction_VB(count,:) = this_fraction_VB;
    
end

end

%% range over r and gamma
updown_amt = 0.002;
r_vals_range = linspace(r_pt-updown_amt,r_pt+updown_amt,100);
gam_vals_range = 10.^linspace(-6,-3,200);

% r_fixed_vector = r_fixed*ones(size(gam_vals_range));
r_corresponding_tradeoff = zeros(size(gam_vals_range));
% r_corresponding_tradeoff_pt = zeros(size(gam_vals_range));

% fixed equilibrium: from which to invade
lysogen_equilibrium_fixed = (K*(r_fixed-gam_fixed-d)/r_fixed);
% lysogen_equilibrium_fixed_pt = (K*(r_pt-gam_pt-d)/r_pt);
% lysogen_equilibrium_fixed_pt_drawline = (K*(r_pt_line-gam_pt_line-d)/r_pt_line);

% lysogen_equilibrium_fixed_nophage = (K*(r_fixed-d)/r_fixed);
phage_equilibrium_fixed = bet*gam_fixed*lysogen_equilibrium_fixed/(phi*lysogen_equilibrium_fixed+m);

for count_r = 1:length(r_vals_range)
    
    this_r = r_vals_range(count_r);
    
    for count_gam = 1:length(gam_vals_range)
        
        this_gam = gam_vals_range(count_gam);
        
        this_phage_equilibrium = bet*this_gam*lysogen_equilibrium_fixed/(phi*lysogen_equilibrium_fixed+m);
        phage_equilibrium_vary(count_r,count_gam) = this_phage_equilibrium;
        this_lysogen_reproduction_number_lysogenic = this_r*(1-lysogen_equilibrium_fixed/K)/(d+this_gam+phi*this_phage_equilibrium);
        lysogen_reproduction_number_lysogenic(count_r,count_gam) = this_lysogen_reproduction_number_lysogenic;
%         this_lysogen_reproduction_number_lysogenic_nophage = this_r*(1-lysogen_equilibrium_fixed_nophage/K)/(d+this_gam);
%         lysogen_reproduction_number_lysogenic_nophage(count_r,count_gam) = this_lysogen_reproduction_number_lysogenic_nophage;
        
%         this_lysogen_reproduction_number = max(this_lysogen_reproduction_number_lysogenic,this_lysogen_reproduction_number_lytic);
%         lysogen_reproduction_number(count_r,count_gam) = this_lysogen_reproduction_number;
        
        
        if count_r == 1
            
            % condition s.t. L_A^* is fixed
            r_corresponding_tradeoff(count_r,count_gam) =  (this_gam + d)/(1-lysogen_equilibrium_fixed/K);
%             r_corresponding_tradeoff_pt(count_r,count_gam) =  (this_gam + d)/(1-lysogen_equilibrium_fixed_pt/K);
%             r_corresponding_tradeoff_pt_drawline(count_r,count_gam) =  (this_gam + d)/(1-lysogen_equilibrium_fixed_pt_drawline/K);
            
            r_corresponding_lysogenicR0equalsone(count_r,count_gam) =  (d+phi*this_phage_equilibrium+this_gam)/(1-lysogen_equilibrium_fixed/K);
            r_corresponding_lysogenicR0equalsone_nophage(count_r,count_gam) =  r_fixed;
            
        end
        
    end
end



%% plot R0 wrt r vs. gamma
f1 = figure(1); set(f1, 'Position', [100 800 1200 450]);
subplot(1,2,1);
imagesc(gam_vals_range,r_vals_range,lysogen_reproduction_number_lysogenic); hold on;
% plot(gam_vals_range,r_fixed_vector,'Color',my_rgb_colors(1,:),'linewidth',2); hold on;
plot(gam_vals_range,r_corresponding_lysogenicR0equalsone,'k--','linewidth',4); hold on;
plot(gam_vals_range,r_corresponding_tradeoff,'Color',[0.5 0.5 0.5],'linewidth',4); hold on;
% plot(gam_vals_range,r_corresponding_tradeoff_pt_drawline,'--','Color',[0.5 0.5 0.5],'linewidth',2); hold on;
% plot(gam_vals_range,r_corresponding_tradeoff_pt,'Color',[0 0 0],'linewidth',2); hold on;
plot(10^-6,r_fixed,'.','Color',[0 0 0],'MarkerSize',45);
plot(10^-6,r_pt,'.','Color',my_rgb_colors(2,:),'MarkerSize',45);
plot(10^-4,r_pt,'.','Color',[0.5 0.5 0.5],'MarkerSize',50);
plot(10^-4,r_pt,'.','Color',my_rgb_colors(3,:),'MarkerSize',35);

% plot(gam_vals_range,r_corresponding_lysogenicR0equalsone_nophage,'Color',[0 0 0],'linewidth',2); hold on;

% if isempty(ind)~=1
% plot(gam_vals_range,r_corresponding_equalR0s,'--','Color',[0 0 0 ],'linewidth',2); hold on;
% end

set(gca,'YDir','normal');
cb1 = colorbar;
% cb1.Label.String = 'CFU/mL';
xlim([1e-6 1e-3]);
ylim([r_pt-updown_amt+0.001 r_pt+updown_amt]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('(Resident) Lysogen Induction Rate, $\gamma_A$ (hr$^{-1}$)','interpreter','latex');
ylabel('(Invading) Lysogen Growth Rate, $r_B$ (hr$^{-1}$)','interpreter','latex');
title('Lysogenic Reproduction Number');
f2=gca;
% f2.CLim = [1 max(max(lysogen_reproduction_number_lysogenic))];
f2.CLim = [1 1.002];
f2.XScale = 'log';
% f1.YScale = 'log';
f2.LineWidth = 1;
f2.FontSize = 20;
f2.FontWeight = 'normal';
f2.FontName = 'Times New Roman';
f2.ColorScale = 'linear';


subplot(1,2,2);
lg(1) = semilogy(t_span,LA_traj(2,:),'Color',[0.5 0.5 0.5],'linewidth',4); hold on;
lg(2) = semilogy(t_span,LB_traj(2,:),'Color',my_rgb_colors(3,:),'linewidth',4); hold on;
lg(3) = semilogy(t_span,LA_traj(1,:),'k','linewidth',4); hold on;
lg(4) = semilogy(t_span,LB_traj(1,:),'Color',my_rgb_colors(2,:),'linewidth',4); hold on;



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



legend_char1 = ['Resident Lysogen, $L_A$ ($\gamma_A > 0$)'];
legend_char2 = ['Invading Lysogen, $L_B$ ($\gamma_A > 0$)'];
legend_char3 = ['Resident Lysogen, $L_A$ ($\gamma_A = 0$)'];
legend_char4 = ['Invading Lysogen, $L_B$ ($\gamma_A = 0$)'];

legend(lg,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Location','SouthWest','FontSize',13);

% %% collect results
% results.r_vals_range = r_vals_range;
% results.gam_vals_range = gam_vals_range;
% results.r_corresponding_tradeoff =r_corresponding_tradeoff;
% 
% results.lysogen_equilibrium_fixed = lysogen_equilibrium_fixed;
% results.phage_equilibrium_fixed = phage_equilibrium_fixed;
% 
% results.lysogen_equilibrium = lysogen_equilibrium;
% results.phage_equilibrium = phage_equilibrium;
% 
% 
% results.lysogen_reproduction_number_lysogenic = lysogen_reproduction_number_lysogenic;
% results.lysogen_reproduction_number_lytic = lysogen_reproduction_number_lytic;
% results.lysogen_reproduction_number = lysogen_reproduction_number;
% 
% results.r_corresponding_lysogenicR0equalsone = r_corresponding_lysogenicR0equalsone;
% 
% results.ind = ind;
% results.updown_amt = updown_amt;
% r_corresponding_equalR0s

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

