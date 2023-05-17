% void = main_plt_single_LVmodel_022823(void)
% plots for single LV-model

%%
clear all; close all; clc;

%% want to save?
save_fig_ans = 0;
% 0: don't save
% 1: save

filename = 'LV_phaseplane_overtime.mat';
folder_location = './sim_data/';
load(strcat(folder_location,filename));
    
figure_name = 'LV_phaseplane_overtime.eps';

my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;
default_rgb_colors = [0, 0, 0; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330; 0.4660, 0.6740, 0.1880];


% read in parameter structure
% parameters
% params.conversion_efficiency = conversion_efficiency;
% params.d_R = d_R;
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
dt = params.dt;
t_span = params.t_span;
t_end = params.t_end;
% params.flask_volume = flask_volume; %%flask volume in mL
% params.L0 = L0;
% params.V0 = V0;


%% read in results
Lysogen_equilibrium_nonzero_vector = results.Lysogen_equilibrium_nonzero_vector;
phage_nullcline = results.phage_nullcline;
phage_equilibrium_nonzero = results.phage_equilibrium_nonzero;
Lysogen_equilibrium_nonzero = results.Lysogen_equilibrium_nonzero;
L_traj = results.L_traj;
V_traj = results.V_traj;
x_values = results.x_values;
y_values = results.y_values;
t_span = results.t_span;
init_conds_range = results.init_conds_range;



%% plot L-V phase plane
f1 = figure(1); set(f1, 'Position', [200 500 1200 450]);
subplot(1,2,1);
% semilogy(x_values,Lysogen_equilibrium_zero,'Color',my_rgb_colors(1,:)); hold on;
h(1)=loglog(x_values,Lysogen_equilibrium_nonzero_vector,'--','Color',my_rgb_colors(2,:),'linewidth',2); hold on;
h(2)=loglog(x_values,phage_nullcline,'--','Color',my_rgb_colors(1,:),'linewidth',2); hold on;
% loglog(phage_nullcline_asymptote*ones(size(y_values)),y_values,'k--','linewidth',1); hold on;
loglog(V_traj(1,:),L_traj(1,:),'k','linewidth',4); hold on;
loglog(V_traj(2,:),L_traj(2,:),'Color',default_rgb_colors(4,:),'linewidth',4); hold on;
loglog(V_traj(3,:),L_traj(3,:),'Color',default_rgb_colors(3,:),'linewidth',4); hold on;
loglog(V_traj(4,:),L_traj(4,:),'Color',default_rgb_colors(2,:),'linewidth',4); hold on;
h(3)=loglog(init_conds_range(2,:),init_conds_range(1,:),'ko','MarkerSize',12,'linewidth',2.5); hold on;
h(4)=loglog(phage_equilibrium_nonzero,Lysogen_equilibrium_nonzero,'.','Color',[0.5,0.5,0.5],'MarkerSize',50); hold on;

xlim([10^3 1e7]);
ylim([10^4 10^10]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Virus Population, $V_A$ (PFU/mL)','interpreter','latex');
ylabel('Lysogen Population, $L_A$ (CFU/mL)','interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 24;
f1.FontName = 'Times New Roman';
f1.FontWeight = 'bold';


% include legend
legend_char1 = ['$L$-Nullcline'];
legend_char2 = ['$V$-Nullcline'];
legend_char3 = ['Initial Condition'];
legend_char4 = ['$L$-$V$-Equilibirum'];
legend(h,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Location','SouthEast','FontSize',20);

%% plot L,V over time
subplot(1,2,2);
% g(3) = semilogy(t_span,Lysogen_equilibrium_nonzero*ones(size(t_span)),'Color',[0.5,0.5,0.5],'linewidth',2); hold on;
% semilogy(t_span,phage_equilibrium_nonzero*ones(size(t_span)),'Color',[0.5,0.5,0.5],'linewidth',2); hold on;
g(1)=semilogy(t_span,L_traj(1,:),'Color',default_rgb_colors(1,:),'linewidth',4); hold on;
g(2)=semilogy(t_span,V_traj(1,:),'--','Color',default_rgb_colors(1,:),'linewidth',4); hold on;
semilogy(t_span,L_traj(2,:),'Color',default_rgb_colors(4,:),'linewidth',4); hold on;
semilogy(t_span,V_traj(2,:),'--','Color',default_rgb_colors(4,:),'linewidth',4); hold on;
semilogy(t_span,L_traj(3,:),'Color',default_rgb_colors(3,:),'linewidth',4); hold on;
semilogy(t_span,V_traj(3,:),'--','Color',default_rgb_colors(3,:),'linewidth',4); hold on;
semilogy(t_span,L_traj(4,:),'Color',default_rgb_colors(2,:),'linewidth',4); hold on;
semilogy(t_span,V_traj(4,:),'--','Color',default_rgb_colors(2,:),'linewidth',4); hold on;

xlim([0 t_end]);
ylim([10^4 10^10]);

xlabel('Time (hr)');
ylabel('Population Density');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 24;
f1.FontName = 'Times New Roman';
f1.FontWeight = 'bold';



% include legend
% legend_char3_g = ['Equilibria'];
legend_char1_g = ['Lysogen population'];
legend_char2_g = ['Virus population'];
% legend(g,{legend_char1_g,legend_char2_g,legend_char3_g}, 'Interpreter','Latex','Location','NorthEast','FontSize',16);
% legend(g,{legend_char1_g,legend_char2_g}, 'Interpreter','Latex','Location','NorthEast','FontSize',24);


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

