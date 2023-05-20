% void = main_get_ss_LVmodel_022823(void)
% steady-state analysis of single LV-model

%%
clear all; close all; clc;


%% save figure?
save_ans_Fig = 1;
% 0: don't save
% 1: save

% figure_name = 'LV_defense_r_gamma_LAinductionHigh';
figure_name = 'LV_defense_r_gamma_LAinductionLow';

%% load file
% infile = 'LV_defense_r_gamma_LAinductionHigh.mat';
infile = 'LV_defense_r_gamma_LAinductionLow.mat';
load(strcat('./sim_data/',infile));

my_rgb_colors = [78 132 193; 209 109 106; 236 180 118; 127 127 127]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];

%% variables
y_traj = results.y_traj;
t_span = params.t_span;
t_end = params.t_end;

LA_traj = y_traj(:,1)';
VA_traj = y_traj(:,2)';
LB_traj = y_traj(:,3)';
VB_traj = y_traj(:,4)';


%% plot R0 wrt r vs. gamma
f1 = figure(1); set(f1, 'Position', [200 500 600 450]);
g(1) = semilogy(t_span,LA_traj,'k','linewidth',4); hold on;
g(2) = semilogy(t_span,VA_traj,'k--','linewidth',4); hold on;
g(3) = semilogy(t_span,VB_traj,'--','Color',my_rgb_colors(1,:),'linewidth',4); hold on;
g(4) = semilogy(t_span,LB_traj,'Color',my_rgb_colors(1,:),'linewidth',4); hold on;
xlim([0 t_end]);
ylim([10^0 10^10]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
% xlabel('Time (hr)');
% ylabel('Population Density');
% title('Basic Reproduction Number - lysogenic path, $\mathcal R_0$','interpreter','latex');
f2=gca;
% f2.CLim = [1 max(max(lysogen_reproduction_number_lysogenic))];
% f2.XScale = 'log';
% f1.YScale = 'log';
f2.LineWidth = 1.5;
f2.FontSize = 28;
f2.FontWeight = 'bold';
f2.FontName = 'Times New Roman';
f2.ColorScale = 'linear';

% legend_char1 = ['Resident Lysogen, $L_A$ ($\gamma_A > 0$)'];
% legend_char2 = ['Invading Lysogen, $L_B$ ($\gamma_A > 0$)'];
% legend_char3 = ['Resident Lysogen, $L_A$ ($\gamma_A = 0$)'];
% legend_char4 = ['Invading Lysogen, $L_B$ ($\gamma_A = 0$)'];

% legend(lg,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Location','SouthWest','FontSize',13);

%% collect results
results.y_traj = y_traj; 

%% save figure
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

