% void = main_LVmodel_012423(void)
% simulate Lysogen-phage model 

%%
clear all; close all; clc;

%% want to save?
save_ans_fig = 0;
% 0: don't save
% 1: save

figure_name = 'LV_phaseplane_overtime.eps';
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


%% Simulate model
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% L0 = 1e4; 
% V0 = 200; 
% L0,V0 - 4 initial conditions
init_conds_range = [2e4,6e4,0.8e6,0.5e6;200,800,900,400];
% init_conds_range = [2e4, 0.9e6;200,200];

for count = 1:length(init_conds_range)

this_init_conds = init_conds_range(:,count);
    
% this_L_traj = zeros(1,length(t_span));
% this_V_traj = zeros(1,length(t_span));

[t,y_traj] = ode45(@(t,y)simulate_LVmodel(t,y,params), params.t_span, this_init_conds, options);
% 
% % reshape from vector to matrix form
% S_traj_discrete = reshape(y_traj(:,1:(m*n)),length(t_span),m,n);
% E_traj_discrete = reshape(y_traj(:,(m*n+1):(2*m*n)),length(t_span),m,n);
% I_traj_discrete = reshape(y_traj(:,(2*m*n+1):(3*m*n)),length(t_span),m,n);
this_L_traj = y_traj(:,1)';
this_V_traj = y_traj(:,2)';

L_traj(count,:) = this_L_traj;
V_traj(count,:) = this_V_traj;

end

%% nullclines
x_values = linspace(0,1e3,100);%10.^linspace(1,5,100);

Lysogen_equilibrium_zero = zeros(size(x_values));
Lysogen_equilibrium_nonzero = (K*(r-eta-d)/r);
Lysogen_equilibrium_nonzero_vector = Lysogen_equilibrium_nonzero*ones(size(x_values));
phage_nullcline = m*x_values./(bet*eta-phi*x_values);%bet*eta*x_values./(phi*x_values+m);
phage_equilibrium_nonzero = bet*eta*Lysogen_equilibrium_nonzero/(phi*Lysogen_equilibrium_nonzero+m);

%% plot phase L-V phase plane
f1 = figure(1); set(f1, 'Position', [200 500 1200 450]);
subplot(1,2,1);
% semilogy(x_values,Lysogen_equilibrium_zero,'Color',my_rgb_colors(1,:)); hold on;
h(1)=semilogy(x_values,Lysogen_equilibrium_nonzero_vector,'Color',my_rgb_colors(2,:),'linewidth',2); hold on;
h(2)=semilogy(x_values,phage_nullcline,'Color',my_rgb_colors(1,:),'linewidth',2); hold on;
h(3)=semilogy(phage_equilibrium_nonzero,Lysogen_equilibrium_nonzero,'k.','MarkerSize',30); hold on;
h(4)=semilogy(V_traj(1,:),L_traj(1,:),'k','linewidth',2); hold on;
semilogy(V_traj(2:end,:)',L_traj(2:end,:)','linewidth',2); hold on;

xlim([0 1000]);
ylim([10^4 10^6]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Virus Population, $V$','interpreter','latex'); 
ylabel('Lysogen Population, $L$','interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';


% include legend
legend_char1 = ['$L$-nullcline'];
legend_char2 = ['$V$-nullcline'];
legend_char3 = ['$L$-$V$-equilibirum'];
legend_char4 = ['Trajectories'];
%         legend_char4 = ['$T_s = ', num2str(1/gamma_s4),'$, $T_a = ', num2str(1/gamma_a4),'$'];
legend(h,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Location','SouthEast','FontSize',16);


subplot(1,2,2);
% semilogy(x_values,Lysogen_equilibrium_zero,'Color',my_rgb_colors(1,:)); hold on;
g(1)=semilogy(t_span,L_traj(1,:),'Color',default_rgb_colors(1,:),'linewidth',2); hold on;
g(2)=semilogy(t_span,V_traj(1,:),'--','Color',default_rgb_colors(1,:),'linewidth',2); hold on;
% semilogy(t_span,L_traj,'Color',default_rgb_colors,'linewidth',2); hold on;
% semilogy(t_span,V_traj,'--','Color',default_rgb_colors,'linewidth',2); hold on;
semilogy(t_span,L_traj(2,:),'Color',default_rgb_colors(4,:),'linewidth',2); hold on;
semilogy(t_span,V_traj(2,:),'--','Color',default_rgb_colors(4,:),'linewidth',2); hold on;
semilogy(t_span,L_traj(3,:),'Color',default_rgb_colors(3,:),'linewidth',2); hold on;
semilogy(t_span,V_traj(3,:),'--','Color',default_rgb_colors(3,:),'linewidth',2); hold on;
semilogy(t_span,L_traj(4,:),'Color',default_rgb_colors(2,:),'linewidth',2); hold on;
semilogy(t_span,V_traj(4,:),'--','Color',default_rgb_colors(2,:),'linewidth',2); hold on;

xlim([0 t_end]);
ylim([10^2 10^6]);
% axis([x_values(1) x_values(end) 0 2*10^5]);
xlabel('Time (hours)'); 
ylabel('Population Density');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';


% include legend
legend_char1_g = ['Lysogen'];
legend_char2_g = ['Virus'];
%         legend_char4 = ['$T_s = ', num2str(1/gamma_s4),'$, $T_a = ', num2str(1/gamma_a4),'$'];
legend(g,{legend_char1_g,legend_char2_g}, 'Interpreter','Latex','Location','SouthEast','FontSize',16);

% 
% % S,E,I,L,V Trajectories
% S_traj = sum(S_traj_discrete,[2,3]); % population sizes
% E_traj = sum(E_traj_discrete,[2,3]);
% I_traj = sum(I_traj_discrete,[2,3]);
% L_traj = sum(L_traj_discrete,[2,3]);
% 
% 

% 
% %% Plotting
% f2=figure(2); set(f2, 'Position', [900   50   400   930]);
% 
% this_p(3)=semilogy(params.t_span, L_traj,'Color',my_rgb_colors(3,:),'LineWidth',2); hold on;
% this_p(4)=semilogy(params.t_span, V_traj,'k','LineWidth',2); hold on;
% 
% 
% axis([0 t_end 10^0 10^10]);
% xlabel('Time (hours)'); ylabel({'Population Size'});
% title('Dynamics')
% legend(this_p,{'S','I','L', 'V'},'Location','SouthWest');
% set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');


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
