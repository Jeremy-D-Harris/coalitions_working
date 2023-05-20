% void = main_get_ss_LVmodel_022823(void)
% steady-state analysis of single LV-model

%%
clear all; close all; clc;


%% save figure?
save_ans_Fig = 1;
% 0: don't save
% 1: save

% figure_name= 'LV_offense_lagtime50hrs';
figure_name= 'LV_offense_lagtime100hrs';

my_rgb_colors = [78 132 193; 209 109 106; 236 180 118; 127 127 127]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];

% filename = 'LV_offense_lagtime50hrs.mat';
filename = 'LV_offense_lagtime100hrs.mat';

load(strcat('./sim_data/',filename));

y_traj1 = results.y_traj1;
y_traj2 = results.y_traj2;

y_traj(1,:) = [y_traj1(:,1)',y_traj2(:,1)'];
y_traj(2,:) = [y_traj1(:,2)',y_traj2(:,2)'];
y_traj(3,:) = [zeros(size(y_traj1(:,3)))',y_traj2(:,3)'];
y_traj(4,:) = [y_traj1(:,4)',y_traj2(:,4)'];

LA_traj = y_traj(1,:);
VA_traj = y_traj(2,:);
LB_traj = y_traj(3,:);
VB_traj = y_traj(4,:);



%% collect results
t_span = params.t_span;
t_end = params.t_end;

results.y_traj1 = y_traj1;
results.y_traj2 = y_traj2;


%% plot R0 wrt r vs. gamma
f1 = figure(1); set(f1, 'Position', [200 800 600 450]);

% subplot(2,1,1);
lg(1) = semilogy(t_span,LA_traj,'Color',default_rgb_colors(1,:),'linewidth',4); hold on;
lg(2) = semilogy(t_span,VA_traj,'--','Color',default_rgb_colors(1,:),'linewidth',4); hold on;
lg(3) = semilogy(t_span,LB_traj,'Color',my_rgb_colors(1,:),'linewidth',4); hold on;
lg(4) = semilogy(t_span,VB_traj,'--','Color',my_rgb_colors(1,:),'linewidth',4); hold on;
semilogy(t_span(202),LB_traj(202),'.','Color',my_rgb_colors(1,:),'MarkerSize',60); hold on;
% semilogy(time_pt*ones(1,10),10.^linspace(0,10,10),'k--','linewidth',2); hold on;
% time_pt


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
f2.LineWidth = 1;
f2.FontSize = 25;
f2.FontWeight = 'bold';
f2.FontName = 'Times New Roman';
f2.ColorScale = 'linear';




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

