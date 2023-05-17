% void = main_LVmodel_012423(void)
% simulate Lysogen-phage model 

%%
clear all; close all; clc;

%% want to save?
save_ans_fig = 0;
% 0: don't save
% 1: save

figure_name = 'LV_steadystate_eta_r.eps';
% figure_name = 'LV_phaseplane_notraj.eps';

%  = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];
my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];

%% system parameters (units of hours, micrograms and mL).
% Assumes the system is a 500 mL flask running for ~ 24hr;
% conversion_efficiency = 5e-7; %ug/cell
% d_R = 0; % per hour
% mu_max = 1.2; % growth rate (per hour)
r = 1; % growth rate (per hour)
K = 3e5;
% R_in = 5; %ug/mL
% d_S = .2; % death rate susceptibles (per hour)
% d_E = .2; % death rate exposed (per hour)
d = 0.05; % death rate lysogens (per hour)
% d_I = .2; % death rate infected (per hour)
% lam = 2; % commitment rate (per hour)
eta = 1e-5; % lysis rate (per hour)
bet = 10; % burst size
phi = 1e-10; %3.4e-10; % adsorption rate (mL/hr)
m = 1/24; % virus washout (per hour)
% alpha_s = 0; % selection coefficient: alpha_s>1 corresponds to advantage of lysogen over susceptible
% J = 0; %ug/mL-h

rng(1);

%simulation parameters:
dt = 0.1; % hours
t_end = 48; % hours
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
params.eta = eta;
params.bet = bet;
params.phi = phi;
params.m = m;
% params.alpha_s = alpha_s;
% params.J = J;
params.dt = dt;
params.t_span = t_span;
% params.flask_volume = flask_volume; %%flask volume in mL
% params.L0 = L0;
% params.V0 = V0;


%% range over r and eta

r_vals_range = linspace(0.1,2,200);
inv_eta_vals_range = 10.^linspace(0,3,200);
r_const_corresponding = zeros(size(inv_eta_vals_range));

r_const = 1;
eta_const = 0.05;
% r_minus_eta_const = r_const - eta_const;
lysogen_equilibrium_const = (K*(r_const-eta_const-d)/r_const);
phage_equilibrium_const = bet*eta_const*lysogen_equilibrium_const/(phi*lysogen_equilibrium_const+m);

lysogen_equilibrium = zeros(length(r_vals_range),length(inv_eta_vals_range));
phage_equilibrium = zeros(length(r_vals_range),length(inv_eta_vals_range));

for count_r = 1:length(r_vals_range)
    this_r = r_vals_range(count_r);
    for count_eta = 1:length(inv_eta_vals_range)
        
        this_eta = 1/inv_eta_vals_range(count_eta);

        this_lysogen_equilibrium = (K*(this_r-this_eta-d)/this_r);
        lysogen_equilibrium(count_r,count_eta) = this_lysogen_equilibrium;
        phage_equilibrium(count_r,count_eta) = bet*this_eta*this_lysogen_equilibrium/(phi*this_lysogen_equilibrium+m);
        
        if count_r == 1
%         r_const_corresponding(count_eta) = r_minus_eta_const + this_eta;
        r_const_corresponding(count_eta) = (this_eta+d)*K/(K-lysogen_equilibrium_const);
        end
        
    end
end


%% plot phase L-V phase plane
f1 = figure(1); set(f1, 'Position', [200 500 1200 450]);
subplot(1,2,1);
imagesc(inv_eta_vals_range,r_vals_range,lysogen_equilibrium); hold on;
plot(inv_eta_vals_range,r_const_corresponding,'k','linewidth',2); hold on;
plot(1/eta_const,r_const,'k.','MarkerSize',20);

set(gca,'YDir','normal');
cb1 = colorbar;
% cb1.Label.String = 'log10(C)';
% xlim([1e5 1e8]);
ylim([0.1 2]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Lysogen Stability, $1/\eta$','interpreter','latex'); 
ylabel('Lysogen Growth, $r$','interpreter','latex');
title('Lysogen Steady-state');
f1=gca;
f1.CLim = [1e5 3.e5];
f1.XScale = 'log';
% f1.YScale = 'log';
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
% f1.ColorScale = 'normal';




subplot(1,2,2);
imagesc(inv_eta_vals_range,r_vals_range,phage_equilibrium); hold on;
plot(inv_eta_vals_range,r_const_corresponding,'k','linewidth',2);
plot(1/eta_const,r_const,'k.','MarkerSize',20);
set(gca,'YDir','normal');
cb1 = colorbar;
% cb1.Label.String = 'log10(C)';
% xlim([1e5 1e8]);
ylim([0.1 2]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Lysogen Stability, $1/\eta$','interpreter','latex'); 
ylabel('Lysogen Growth, $r$','interpreter','latex');
title('Phage Steady-state');
f1=gca;
f1.CLim = [1e2 1e7];
f1.XScale = 'log';
% f1.YScale = 'log';
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.ColorScale = 'log';



%% save figure
if save_ans_fig
    
    folder_location = './figures/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('Figure not saved.\n');
    
end

