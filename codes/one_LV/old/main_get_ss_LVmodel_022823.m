% void = main_get_ss_LVmodel_022823(void)
% steady-state analysis of single LV-model

%%
clear all; close all; clc;

%% want to save?
save_ans = 0;
% 0: don't save
% 1: save

filename = 'LV_steadystate.mat';
% figure_name = 'LV_phaseplane_notraj.eps';

my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];

%% system parameters (units of hours, micrograms and mL).
r_fixed = 0.5;
gam_fixed = 1e-4;
r_minus_gam_const = r_fixed - gam_fixed;
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
gam = 1e-6; % lysis rate (per hour)
bet = 10; % burst size
phi = 3.4e-10; %3.4e-10; % adsorption rate (mL/hr)
m = 1/24; % virus washout (per hour)
% alpha_s = 0; % selection coefficient: alpha_s>1 corresponds to advantage of lysogen over susceptible
% J = 0; %ug/mL-h

% rng(1);

%simulation parameters:
% dt = 0.05; % hours
% t_end = 40; % hours
% t_span = transpose(0:dt:t_end); % time
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
% params.dt = dt;
% params.t_span = t_span;
% params.t_end = t_end;
% params.flask_volume = flask_volume; %%flask volume in mL
% params.L0 = L0;
% params.V0 = V0;


%% range over r and eta

r_vals_range = linspace(0.1,1.5,20);
gam_vals_range = 10.^linspace(-6,-2,20);
r_fixed_corresponding = zeros(size(gam_vals_range));

% r_minus_eta_const = r_const - eta_const;
lysogen_equilibrium_fixed = (K*(r_fixed-gam_fixed-d)/r_fixed);
phage_equilibrium_fixed = bet*gam_fixed*lysogen_equilibrium_fixed/(phi*lysogen_equilibrium_fixed+m);

lysogen_equilibrium = zeros(length(r_vals_range),length(gam_vals_range));
phage_equilibrium = zeros(length(r_vals_range),length(gam_vals_range));

for count_r = 1:length(r_vals_range)
    this_r = r_vals_range(count_r);
    for count_gam = 1:length(gam_vals_range)
        
        this_gam = gam_vals_range(count_gam);

        this_lysogen_equilibrium = max((K*(this_r-this_gam-d)/this_r),0);
        lysogen_equilibrium(count_r,count_gam) = this_lysogen_equilibrium;
        phage_equilibrium(count_r,count_gam) = bet*this_gam*this_lysogen_equilibrium/(phi*this_lysogen_equilibrium+m);
        
        if count_r == 1
        r_fixed_corresponding(count_gam) = r_minus_gam_const + (1/gam_fixed)*this_gam;
%         r_fixed_corresponding(count_gam) = (this_gam+d)*K/(K-lysogen_equilibrium_fixed);
        end
        
    end
end


%% plot phase L-V phase plane
f1 = figure(1); set(f1, 'Position', [200 500 1200 450]);
subplot(1,2,1);
imagesc(gam_vals_range,r_vals_range,max(lysogen_equilibrium,10^0)); hold on;
plot(gam_vals_range,r_fixed_corresponding,'k','linewidth',2); hold on;
plot(gam_fixed,r_fixed,'k.','MarkerSize',20);

set(gca,'YDir','normal');
cb1 = colorbar;
% cb1.Label.String = 'log10(C)';
xlim([1e-6 1e-2]);
% ylim([0.1 2]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Lysogen Induction Rate, $\gamma$ (hr$^{-1}$)','interpreter','latex'); 
ylabel('Lysogen Growth Rate, $r$ (hr$^{-1}$)','interpreter','latex');
title('Lysogen Steady-state, $L^*$','interpreter','latex');
f1=gca;
% f1.CLim = [1e7 1.72e8];
f1.XScale = 'log';
% f1.YScale = 'log';
f1.LineWidth = 1;
f1.FontSize = 18;
f1.FontWeight = 'normal';
f1.FontName = 'Times New Roman';
f1.ColorScale = 'log';




subplot(1,2,2);
imagesc(gam_vals_range,r_vals_range,max(phage_equilibrium,10^4)); hold on;
plot(gam_vals_range,r_fixed_corresponding,'k','linewidth',2);
plot(gam_fixed,r_fixed,'k.','MarkerSize',20);
set(gca,'YDir','normal');
cb1 = colorbar;
% cb1.Label.String = 'log10(C)';
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Lysogen Induction Rate, $\gamma$ (hr$^{-1}$)','interpreter','latex'); 
ylabel('Lysogen Growth Rate, $r$ (hr$^{-1}$)','interpreter','latex');
title('Phage Steady-state, $V^*$','interpreter','latex');
f1=gca;
% f1.CLim = [1e4 1e7];
f1.XScale = 'log';
% f1.YScale = 'log';
f1.LineWidth = 1;
f1.FontSize = 18;
f1.FontWeight = 'normal';
f1.FontName = 'Times New Roman';
f1.ColorScale = 'log';



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

