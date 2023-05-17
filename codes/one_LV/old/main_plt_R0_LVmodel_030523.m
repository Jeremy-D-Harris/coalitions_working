% void = main_plt_single_LVmodel_022823(void)
% plots for single LV-model

%%
clear all; close all; clc;

%% want to save?
save_fig_ans = 0;
% 0: don't save
% 1: save

filename = 'LV_R0_r_gamma.mat';
folder_location = './sim_data/';
load(strcat(folder_location,filename));
    
figure_name = 'LV_R0_r_gamma.eps';

my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];


% read in parameter structure
% parameters
% params.conversion_efficiency = conversion_efficiency;
% params.d_R = d_R;
r_fixed = params.r_fixed;
gam_fixed = params.gam_fixed;
% alph = params.alph;
% r = params.r;
K = params.K;
% params.R_in = R_in;
% params.d_S = d_S;
% params.d_E = d_E;
 d = params.d;
% params.d_I = d_I;
% params.lambda = lam;
% gam = params.gam;
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
r_corresponding_tradeoff = results.r_corresponding_tradeoff;

lysogen_equilibrium_fixed = results.lysogen_equilibrium_fixed;
phage_equilibrium_fixed = results.phage_equilibrium_fixed;

lysogen_equilibrium = results.lysogen_equilibrium;
phage_equilibrium = results.phage_equilibrium;


lysogen_reproduction_number_lysogenic = results.lysogen_reproduction_number_lysogenic;
lysogen_reproduction_number_lytic = results.lysogen_reproduction_number_lytic;
lysogen_reproduction_number = results.lysogen_reproduction_number;

r_corresponding_lysogenicR0equalsone = results.r_corresponding_lysogenicR0equalsone;
ind = results.ind;
updown_amt = results.updown_amt;


if 0
%% plot phase L-V phase plane
f1 = figure(1); set(f1, 'Position', [200 500 1200 450]);
subplot(1,2,1);
imagesc(gam_vals_range,r_vals_range,max(lysogen_equilibrium,10^0)); hold on;
plot(gam_vals_range,r_fixed_corresponding,'k','linewidth',3); hold on;
plot(gam_fixed,r_fixed,'k.','MarkerSize',40);

set(gca,'YDir','normal');
cb1 = colorbar;
% cb1.Label.String = 'CFU/mL';
xlim([1e-6 1e-2]);
ylim([0.95 1.05]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Lysogen Induction Rate, $\gamma$ (hr$^{-1}$)','interpreter','latex'); 
ylabel('Lysogen Growth Rate, $r$ (hr$^{-1}$)','interpreter','latex');
title('Lysogen Steady State, $L^*$','interpreter','latex');
f1=gca;
f1.CLim = [1e4 1.75e8];
f1.XScale = 'log';
% f1.YScale = 'log';
f1.LineWidth = 1;
f1.FontSize = 24;
f1.FontWeight = 'normal';
f1.FontName = 'Times New Roman';
f1.ColorScale = 'log';

% xticks([10^-6 10^-5 10^-4 10^-3 10^-2]);
% yticks([0.4 0.8 1.2 1.6 2]);


subplot(1,2,2);
imagesc(gam_vals_range,r_vals_range,max(phage_equilibrium,10^4)); hold on;
plot(gam_vals_range,r_fixed_corresponding,'k','linewidth',3);
plot(gam_fixed,r_fixed,'k.','MarkerSize',40);
set(gca,'YDir','normal');
cb1 = colorbar;
xlim([1e-6 1e-2]);
ylim([0.95 1.05]);
% cb1.Label.String = 'PFU/mL';
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Lysogen Induction Rate, $\gamma$ (hr$^{-1}$)','interpreter','latex'); 
ylabel('Lysogen Growth Rate, $r$ (hr$^{-1}$)','interpreter','latex');
title('Phage Steady State, $V^*$','interpreter','latex');
f1=gca;
f1.CLim = [1e4 1.75e8];
f1.XScale = 'log';
% f1.YScale = 'log';
f1.LineWidth = 1;
f1.FontSize = 24;
f1.FontWeight = 'normal';
f1.FontName = 'Times New Roman';
f1.ColorScale = 'log';

end
% xticks([10^-6 10^-5 10^-4 10^-3 10^-2]);
% yticks([0.4 0.8 1.2 1.6 2]);
% xticklabels({'x = 0','x = 5','x = 10'});


%% plot R0 wrt r vs. gamma
f2 = figure(2); set(f2, 'Position', [200 800 1200 450]);
subplot(1,2,1);
imagesc(gam_vals_range,r_vals_range,lysogen_reproduction_number_lysogenic); hold on;
plot(gam_vals_range,r_corresponding_lysogenicR0equalsone,'k--','linewidth',1); hold on;
plot(gam_vals_range,r_corresponding_tradeoff,'Color',[0.5 0.5 0.5],'linewidth',2); hold on;
plot(gam_fixed,r_fixed,'.','Color',[0.5 0.5 0.5],'MarkerSize',40);
if isempty(ind)~=1
plot(gam_vals_range,r_corresponding_equalR0s,'--','Color',[0 0 0 ],'linewidth',2); hold on;
end

set(gca,'YDir','normal');
cb1 = colorbar;
% cb1.Label.String = 'CFU/mL';
xlim([1e-6 1e-2]);
ylim([r_fixed-updown_amt r_fixed+updown_amt]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Lysogen Induction Rate, $\gamma$ (hr$^{-1}$)','interpreter','latex');
ylabel('Lysogen Growth Rate, $r$ (hr$^{-1}$)','interpreter','latex');
% title('Basic Reproduction Number - lysogenic path, $\mathcal R_0$','interpreter','latex');
f2=gca;
% f2.CLim = [1 max(max(lysogen_reproduction_number_lysogenic))];
f2.CLim = [1 1.02];
f2.XScale = 'log';
% f1.YScale = 'log';
f2.LineWidth = 1;
f2.FontSize = 22;
f2.FontWeight = 'normal';
f2.FontName = 'Times New Roman';
f2.ColorScale = 'linear';


subplot(1,2,2);
imagesc(gam_vals_range,r_vals_range,lysogen_reproduction_number_lytic); hold on;
% plot(gam_vals_range,r_corresponding_tradeoff,'k','linewidth',2); hold on;
plot(gam_vals_range,r_corresponding_lysogenicR0equalsone,'k--','linewidth',1); hold on;
plot(gam_vals_range,r_corresponding_tradeoff,'Color',[0.5 0.5 0.5],'linewidth',2); hold on;
plot(gam_fixed,r_fixed,'.','Color',[0.5 0.5 0.5],'MarkerSize',40);

set(gca,'YDir','normal');
cb1 = colorbar;
% cb1.Label.String = 'CFU/mL';
xlim([1e-6 1e-2]);
ylim([r_fixed-updown_amt r_fixed+updown_amt]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Lysogen Induction Rate, $\gamma$ (hr$^{-1}$)','interpreter','latex');
ylabel('Lysogen Growth Rate, $r$ (hr$^{-1}$)','interpreter','latex');
% title('Basic Reproduction Number - Lytic path, $\mathcal R_0$','interpreter','latex');
f2=gca;
f2.CLim = [1 4];
f2.XScale = 'log';
% f1.YScale = 'log';
f2.LineWidth = 1;
f2.FontSize = 22;
f2.FontWeight = 'normal';
f2.FontName = 'Times New Roman';
f2.ColorScale = 'linear';

if 0
%% plot R0 wrt r vs. gamma
f3 = figure(3); set(f3, 'Position', [400 200 600 450]);
% subplot(1,2,1);
imagesc(gam_vals_range,r_vals_range,lysogen_reproduction_number); hold on;
plot(gam_vals_range,r_corresponding_lysogenicR0equalsone,'k--','linewidth',1); hold on;

plot(gam_vals_range,r_corresponding_tradeoff,'Color',[0.5 0.5 0.5],'linewidth',2); hold on;
plot(gam_fixed,r_fixed,'.','Color',[0.5 0.5 0.5],'MarkerSize',40);
if isempty(ind)~=1
plot(gam_vals_range,r_corresponding_equalR0s,'--','Color',[0 0 0],'linewidth',2); hold on;
end
% plot(gam_vals_range,r_corresponding_equalR0s,'k','linewidth',2.5); hold on;
% if max(r_lysogenic_lytic_R0_intersection)>0
% plot(gam_vals_range,r_lysogenic_lytic_R0_intersection,'k','linewidth',2.5); hold on;
% end

set(gca,'YDir','normal');
cb1 = colorbar;
% cb1.Label.String = 'CFU/mL';
% xlim([1e-6 1e-2]);
ylim([r_fixed-updown_amt r_fixed+updown_amt]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Lysogen Induction Rate, $\gamma$ (hr$^{-1}$)','interpreter','latex');
ylabel('Lysogen Growth Rate, $r$ (hr$^{-1}$)','interpreter','latex');
title('Basic Reproduction Number, $\mathcal R_0$','interpreter','latex');
f3=gca;
f3.CLim = [1 4];
f3.XScale = 'log';
% f1.YScale = 'log';
f3.LineWidth = 1;
f3.FontSize = 22;
f3.FontWeight = 'normal';
f3.FontName = 'Times New Roman';
f3.ColorScale = 'linear';

end

%% save figure
if save_fig_ans
    
    folder_location = './figures/';
    saveas(f2,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('Figure not saved.\n');
    
end

