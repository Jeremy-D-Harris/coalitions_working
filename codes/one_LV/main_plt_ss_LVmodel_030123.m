% void = main_plt_single_LVmodel_022823(void)
% plots for single LV-model

%%
clear all; close all; clc;

%% want to save?
save_fig_ans = 0;
% 0: don't save
% 1: save

filename = 'LV_steadystate_r_gamma.mat';
folder_location = './sim_data/';
load(strcat(folder_location,filename));
    
figure_name = 'LV_steadystate_r_gamma.eps';

my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];


% read in parameter structure
% parameters
% params.conversion_efficiency = conversion_efficiency;
% params.d_R = d_R;
r_fixed = params.r_fixed;
gam_fixed = params.gam_fixed;
alph = params.alph;
r = params.r;
K = params.K;
% params.R_in = R_in;
% params.d_S = d_S;
% params.d_E = d_E;
 d = params.d;
% params.d_I = d_I;
% params.lambda = lam;
gam = params.gam;
bet = params.bet;
phi = params.phi;
m = params.m;
% params.alpha_s = alpha_s;
% params.J = J;
% dt = params.dt;
% t_span = params.t_span;
% t_end = params.t_end;
% params.flask_volume = flask_volume; %%flask volume in mL
% params.L0 = L0;
% params.V0 = V0;


%% read in results
r_vals_range = results.r_vals_range;
gam_vals_range = results.gam_vals_range;
r_fixed_corresponding =results.r_fixed_corresponding;

lysogen_equilibrium_fixed = results.lysogen_equilibrium_fixed;
phage_equilibrium_fixed = results.phage_equilibrium_fixed;

lysogen_equilibrium = results.lysogen_equilibrium;
phage_equilibrium = results.phage_equilibrium;



%% plot phase L-V phase plane
f1 = figure(1); set(f1, 'Position', [200 500 1200 450]);
subplot(1,2,1);
imagesc(gam_vals_range,r_vals_range,max(lysogen_equilibrium,10^0)); hold on;
plot(gam_vals_range,r_fixed_corresponding,'k','linewidth',5); hold on;
plot(gam_fixed,r_fixed,'k.','MarkerSize',50);

set(gca,'YDir','normal');
cb1 = colorbar;
% cb1.Label.String = 'CFU/mL';
xlim([1e-6 1e-2]);
ylim([0.98 1.02]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Lysogen Induction Rate, $\gamma_A$ (hr$^{-1}$)','interpreter','latex'); 
ylabel('Lysogen Growth Rate, $r_A$ (hr$^{-1}$)','interpreter','latex');
title('Lysogen Steady State, $L_A^*$','interpreter','latex');
f1=gca;
% f1.CLim = [1e4 1.75e8];
f1.XScale = 'log';
% f1.YScale = 'log';
f1.LineWidth = 1;
f1.FontSize = 24;
f1.FontWeight = 'bold';
f1.FontName = 'Times New Roman';
f1.ColorScale = 'log';

% xticks([10^-6 10^-5 10^-4 10^-3 10^-2]);
% yticks([0.4 0.8 1.2 1.6 2]);


subplot(1,2,2);
imagesc(gam_vals_range,r_vals_range,max(phage_equilibrium,10^4)); hold on;
plot(gam_vals_range,r_fixed_corresponding,'k','linewidth',5);
plot(gam_fixed,r_fixed,'k.','MarkerSize',50);
set(gca,'YDir','normal');
cb1 = colorbar;
xlim([1e-6 1e-2]);
ylim([0.98 1.02]);
% cb1.Label.String = 'PFU/mL';
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Lysogen Induction Rate, $\gamma_A$ (hr$^{-1}$)','interpreter','latex'); 
ylabel('Lysogen Growth Rate, $r_A$ (hr$^{-1}$)','interpreter','latex');
title('Phage Steady State, $V_A^*$','interpreter','latex');
f1=gca;
f1.CLim = [1e4 1.75e8];
f1.XScale = 'log';
% f1.YScale = 'log';
f1.LineWidth = 1;
f1.FontSize = 24;
f1.FontWeight = 'bold';
f1.FontName = 'Times New Roman';
f1.ColorScale = 'log';

% xticks([10^-6 10^-5 10^-4 10^-3 10^-2]);
% yticks([0.4 0.8 1.2 1.6 2]);
% xticklabels({'x = 0','x = 5','x = 10'});


%% save figure
if save_fig_ans
    
    folder_location = './figures/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('Figure not saved.\n');
    
end

