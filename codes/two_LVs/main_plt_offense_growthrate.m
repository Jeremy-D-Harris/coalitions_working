% void = main_get_ss_LVmodel_022823(void)
% steady-state analysis of single LV-model

clear all; close all; clc;

%% save figure?
save_ans_Fig = 1;
% 0: don't save
% 1: save

% figure_name= 'LV_offense_lagtime50hrs';
figure_name= 'LV_offense_growthrate_lagtime';

my_rgb_colors = [78 132 193; 209 109 106; 236 180 118; 127 127 127]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];



% filename = 'LV_offense_lagtime50hrs.mat';
filename = 'LV_offense_growthrate_lagtime.mat';
load(strcat('./sim_data/',filename));


t_span = params.t_span;
t_end = params.t_end;
r_corresponding_lysogenicR0equalsone_withVB = results.r_corresponding_lysogenicR0equalsone_withVB;



%% plot R0 wrt r vs. gamma
f1 = figure(1); set(f1, 'Position', [200 800 600 450]);
plot(t_span,r_corresponding_lysogenicR0equalsone_withVB,'k','linewidth',4); hold on;
plot(50,0.99,'.','Color',my_rgb_colors(4,:),'MarkerSize',60); hold on;
plot(100,0.99,'.','Color',my_rgb_colors(1,:),'MarkerSize',60); 
% plot(t_span(ind(1))*ones(size(t_span)),linspace(0.99,1.01,length(t_span)),'Color',[0.5 0.5 0.5],'linewidth',2);
plot(t_span,1.01*ones(size(t_span)),'--','Color',[0.5 0.5 0.5],'linewidth',4); hold on;


xlim([0 t_end]);
ylim([0.99 1.03]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
% xlabel('Time (hr)');
% ylabel('Population Density');
% title('Basic Reproduction Number - lysogenic path, $\mathcal R_0$','interpreter','latex');
f2=gca;
% f2.CLim = [1 max(max(lysogen_reproduction_number_lysogenic))];
% f2.XScale = 'log';
% f1.YScale = 'log';
f2.LineWidth = 1;
f2.FontSize = 28;
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

