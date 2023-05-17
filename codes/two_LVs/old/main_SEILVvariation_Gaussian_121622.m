% simulate Lysogen-phage model 

%%

clear all; close all; clc;

%% want to save?
save_ans = 0;
% 0: don't save
% 1: save

write_video_ans = 1;
% 0: don't save
% 1: save
video_filename = 'variation_susc_trans.gif';

% filename = 'negbinomial_Poissonlike.mat';
filename = 'negbinomial_overdispersed_uncorr.mat';
% filename = 'negbinomial_overdispersed_poscorr.mat';

%  = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];
my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;


%% system parameters (units of hours, micrograms and mL).
% Assumes the system is a 500 mL flask running for ~ 24hr;
% conversion_efficiency = 5e-7; %ug/cell
% d_R = 0; % per hour
mu_max = 1.2; % growth rate (per hour)
K = 10^8;
% R_in = 5; %ug/mL
d_S = .2; % death rate susceptibles (per hour)
d_E = .2; % death rate exposed (per hour)
d_L = .2; % death rate lysogens (per hour)
d_I = .2; % death rate infected (per hour)
lam = 2; % commitment rate (per hour)
eta = 1; % lysis rate (per hour)
bet = 50; % burst size
phi = 3.4e-10; % adsorption rate (mL/hr)
m_V = 1/24; % virus washout (per hour)
% alpha_s = 0; % selection coefficient: alpha_s>1 corresponds to advantage of lysogen over susceptible
% J = 0; %ug/mL-h

rng(1);

%simulation parameters:
dt = 1/60; % hours
t_end = 48; % hours
t_span = transpose(0:dt:t_end); % time
% NRuns = 100;
% p = linspace(0,1,6); %dilution factor

% R0 = [10 logspace(2,8,6)]; %initial resource amount in ug/mL ( 500 mL flask)
S0 = 1e8; %Initial concentration of susceptibles in flask (per mL)
V0 = 1e6; %initial concentration of virus in flask (per mL)
% flask_volume = 500; %volume in mL


%% set up parameter structure
% parameters
% params.conversion_efficiency = conversion_efficiency;
% params.d_R = d_R;
params.mu_max = mu_max;
params.K = K;
% params.R_in = R_in;
params.d_S = d_S;
params.d_E = d_E;
params.d_L = d_L;
params.d_I = d_I;
params.lambda = lam;
params.eta = eta;
params.bet = bet;
params.phi = phi;
params.m_V = m_V;
% params.alpha_s = alpha_s;
% params.J = J;
params.dt = dt;
params.t_span = t_span;
% params.flask_volume = flask_volume; %%flask volume in mL
params.S0 = S0;
params.V0 = V0;


options = odeset('AbsTol',1e-8,'RelTol',1e-8,'NonNegative',1:6);


%% Initialize Joint Distributions & calculate Marginal
% can use other bivariation distributions, including nonparametric !
% x1 = gamma, x2 = q

% % discrete mesh
n = 100; % number induction classes
m = n; % number lysogeny probability classes

% gamma_vector = 10.^linspace(-3,0,26); %per hour
% n = 100;
gamma_vector = linspace(10^-3,1,n); %per hour
q_vector = linspace(0,1,m); % probability
dgamma = (1-10^-3)/n ; dq = 1/n;
q_vector_plt=q_vector;
q_vector_plt(1)=10^-3;

[gamma_mesh,q_mesh] = meshgrid(gamma_vector,q_vector);
X = [gamma_mesh(:) q_mesh(:)];

% lysogen parameters
% specify means of marginals
mean_gamma = 10^(-2); % gamma
mean_q = 0; % q
mean_S = [mean_gamma mean_q]; %means for S

variance_gamma = 0.0005;
variance_q = 0.001;
rho_S = 0.2; % 0, 0.6
Sigma_matrix = [variance_gamma rho_S*sqrt(variance_gamma*variance_q); rho_S*sqrt(variance_gamma*variance_q) variance_q]; %change off-diagonals for correlation

init_joint_S = mvnpdf(X,mean_S,Sigma_matrix);
init_joint_S = reshape(init_joint_S,length(q_vector),length(gamma_vector));
%normalize here
init_joint_S = init_joint_S/sum(sum(init_joint_S));

init_marginal_gamma_S = sum(init_joint_S);
init_marginal_q_S = sum(init_joint_S,2)';

% calculated mean values
init_mean_gamma_S = sum(gamma_vector.*init_marginal_gamma_S);
init_mean_q_S = sum(q_vector.*init_marginal_q_S);


%calculated values
init_covariance = sum(sum(((q_vector' - init_mean_gamma_S*ones(size(q_vector')))*(gamma_vector- init_mean_q_S*ones(size(gamma_vector)))).*init_joint_S));

init_variance_gamma = sum((gamma_vector - init_mean_gamma_S*ones(size(gamma_vector))).^2.*init_marginal_gamma_S);
init_sd_gamma = sqrt(init_variance_gamma);

init_variance_q = sum((q_vector- init_mean_q_S*ones(size(q_vector))).^2.*init_marginal_q_S);
init_sd_q = sqrt(init_variance_q);

corrcoef = init_covariance/(init_sd_gamma*init_sd_q);

% distribution parameters
params.gamma_mesh=gamma_mesh;
params.q_mesh=q_mesh;
params.mean_gamma = mean_gamma;
params.mean_q = mean_q;
params.Sigma_matrix = Sigma_matrix;
params.init_joint_S = init_joint_S;


% plot initial distributions - if you want!
if 1
    
    
    f1=figure(1); set(f1, 'Position', [200   400   1220   350]);
    subplot(1,3,1);
    
    imagesc(gamma_vector,q_vector,init_joint_S);
    % pcolor(eps,del,joint_S);
    %axis xy;
    set(gca,'YDir','normal');
    colorbar;
    xlim([0 0.2]); ylim([0 0.5]);
    
    xlabel('Induction rate, $\gamma$','interpreter','latex');
    ylabel({'Lysogeny Probability, $q$'},'interpreter','latex');
    %     zlabel('Population Fraction');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    title('Initial trait distribution');
    
    % cbh = colorbar();
    % set color range
    % caxis([0,.2])
    % set ticks
    % set(cbh, 'YTick', [0.001, 0.01, 0.05, .1], ...
    %     'YTickLabel', {'p=0.001', 'p=0.01', 'p=0.05', 'p=.1'})
    % cbh.Label.String = 'Population Density';
    
    
    % f1=figure(1);
    subplot(1,3,2);
    this_q(1)=plot(gamma_vector,init_marginal_gamma_S,'.-','Color',my_rgb_colors(3,:),'LineWidth',2.5,'MarkerSize',20); hold on;
    %     this_q(2)=plot(q_vector_plt,init_marginal_q_S,'.-','Color',my_rgb_colors(2,:),'LineWidth',2.5,'MarkerSize',20); hold on;
    % cat_margs = [marg_eps_S;marg_delta_S];
    % q=bar(eps,cat_margs);
    %     axis([0 0.05 0 0.2]);
    xlim([0 0.2]);
    
    %     axis([0 4 0 1]);
    title('Marginal Distribution');
    
    xlabel('Induction rate $\gamma$','interpreter','latex');
    ylabel({'Population Fraction'});
    
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    subplot(1,3,3);
    %     this_q(1)=plot(gamma_vector,init_marginal_gamma_S,'.-','Color',my_rgb_colors(3,:),'LineWidth',2.5,'MarkerSize',20); hold on;
    this_q(2)=plot(q_vector_plt,init_marginal_q_S,'.-','Color',my_rgb_colors(2,:),'LineWidth',2.5,'MarkerSize',20); hold on;
    % cat_margs = [marg_eps_S;marg_delta_S];
    % q=bar(eps,cat_margs);
    %     axis([0 0.5 0 0.2]);
    xlim([0 0.5]);
    %     axis([0 4 0 1]);
    title('Marginal Distribution');
    
    xlabel('Lysogeny Probability, $q$','interpreter','latex');
    ylabel({'Population Fraction'});
    
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    %     legend(this_q,{'Induction rate','Probability lysogeny'},'Location','NorthEast');
    
end

%%
% pause;
% small perturbation in direction of eigenvector
% eps_perturb = 5e-8;

% params.eps_perturb = eps_perturb;

% params.mu_eps_S = init_mean_gamma_S;
% params.mu_delta_I = init_mu_delta_I;

% eigen_direction_SIR = get_eigendirection_SIRdelta_eps(params);
% 
% S_init = 1 - eps_perturb*abs(eigen_direction_SIR(1));
% I_init = eps_perturb*abs(eigen_direction_SIR(2));
% R_init = eps_perturb*abs(eigen_direction_SIR(3));

% pause;

init_dist_S = init_joint_S.*S0; 
params.init_dist_S = init_dist_S;

init_dist_E = zeros(size(init_dist_S));
params.init_dist_E = init_dist_E;

init_dist_I = zeros(size(init_dist_S));
params.init_dist_I = init_dist_I;

init_dist_L = zeros(size(init_dist_S));
params.init_dist_L = init_dist_L;

%check should equal 1
%sum(sum((init_dist_S + init_dist_E + init_dist_I + init_dist_L),2));

% need to reshape from matrix to vector from in order to pass to ODE function
init_conds = [reshape(init_dist_S,m*n,1); reshape(init_dist_E,m*n,1); ...
    reshape(init_dist_I,m*n,1); reshape(init_dist_L,m*n,1); V0];

%% Simulate model
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

S_traj_discrete =zeros(length(t_span),m,n);
E_traj_discrete = zeros(length(t_span),m,n);
E_traj_discrete = zeros(length(t_span),m,n);
L_traj_discrete = zeros(length(t_span),m,n);
V_traj = zeros(length(t_span),1);

[t,y_traj] = ode45(@(t,y)simulate_SEILV_gamma_q(t,y,params), params.t_span, init_conds, options);

% reshape from vector to matrix form
S_traj_discrete = reshape(y_traj(:,1:(m*n)),length(t_span),m,n);
E_traj_discrete = reshape(y_traj(:,(m*n+1):(2*m*n)),length(t_span),m,n);
I_traj_discrete = reshape(y_traj(:,(2*m*n+1):(3*m*n)),length(t_span),m,n);
L_traj_discrete = reshape(y_traj(:,(3*m*n+1):(4*m*n)),length(t_span),m,n);
V_traj = y_traj(:,(4*m*n+1));

% S,E,I,L,V Trajectories
S_traj = sum(S_traj_discrete,[2,3]); % population sizes
E_traj = sum(E_traj_discrete,[2,3]);
I_traj = sum(I_traj_discrete,[2,3]);
L_traj = sum(L_traj_discrete,[2,3]);


marginal_gamma_L_traj = reshape(sum(L_traj_discrete,2),length(t_span),m);
marginal_q_L_traj = reshape(sum(L_traj_discrete,3),length(t_span),n);

for kk=1:length(t_span)
    % mean epsilon_S(t), delta_I(t)
    mu_gamma_L_traj(kk) = sum(gamma_vector.*marginal_gamma_L_traj(kk,:))./L_traj(kk);
    mu_q_L_traj(kk) = sum(q_vector.*marginal_q_L_traj(kk,:))./L_traj(kk);
    
end


SEILV_gamma_q_traj = zeros(length(t_span),5);
SEILV_gamma_q_traj(:,1) = S_traj;
SEILV_gamma_q_traj(:,2) = E_traj;
SEILV_gamma_q_traj(:,3) = I_traj;
SEILV_gamma_q_traj(:,4) = L_traj;
SEILV_gamma_q_traj(:,5) = V_traj;

% Rt_SIR_traj = get_Rt_SIR_delta_eps(params,SEILV_gamma_q_traj);


for kk = 1:length(t_span)
    
    %     marginal_eps_S_traj(kk,:) = marginal_eps_S_traj(kk,:)/sum(marginal_eps_S_traj(kk,:));
    %     marginal_delta_S_traj(kk,:) = marginal_delta_S_traj(kk,:)/sum(marginal_delta_S_traj(kk,:));
    %
    %     marginal_eps_I_traj(kk,:) = marginal_eps_I_traj(kk,:)/sum(marginal_eps_I_traj(kk,:));
    %     marginal_delta_I_traj(kk,:) = marginal_delta_I_traj(kk,:)/sum(marginal_delta_I_traj(kk,:));
    
    variance_gamma_L_traj(kk) = sum((gamma_vector- mu_gamma_L_traj(kk)*ones(size(gamma_vector))).^2.*marginal_gamma_L_traj(kk,:));
    variance_q_L_traj(kk) = sum(((q_vector- mu_q_L_traj(kk)*ones(size(q_vector))).^2).*marginal_q_L_traj(kk,:));
    
    sd_gamma_L_traj(kk) = sqrt(variance_gamma_L_traj(kk));
    sd_q_L_traj(kk) = sqrt(variance_q_L_traj(kk));
    
end

%dispersion parameters
% dispSe_traj = std_Se_traj./mu_epsilon_S_traj;
% dispId_traj = std_Id_traj./mu_delta_I_traj;

CV_gamma_traj = variance_gamma_L_traj./(mu_gamma_L_traj.^2);
CV_q_traj = variance_q_L_traj./(mu_q_L_traj.^2);


%% Plotting
f2=figure(2); set(f2, 'Position', [900   50   400   930]);

subplot(3,1,1);
this_p(1)=semilogy(params.t_span, S_traj,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;
this_p(2)=semilogy(params.t_span, I_traj,'Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
this_p(3)=semilogy(params.t_span, L_traj,'Color',my_rgb_colors(3,:),'LineWidth',2); hold on;
this_p(4)=semilogy(params.t_span, V_traj,'k','LineWidth',2); hold on;


axis([0 t_end 10^0 10^10]);
xlabel('Time (hours)'); ylabel({'Population Size'});
title('Dynamics')
legend(this_p,{'S','I','L', 'V'},'Location','SouthWest');
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

subplot(3,1,2);
semilogy(params.t_span, mu_gamma_L_traj,'-','Color',my_rgb_colors(1,:),'LineWidth',2); hold on;
% r(1)=plot(params.t_span, mu_gamma_L_traj,'-','Color',my_rgb_colors(1,:),'LineWidth',2); hold on;
% r(2)=plot(params.t_span, mu_q_L_traj,'-','Color',my_rgb_colors(2,:),'LineWidth',2); hold on;

% axis([0 t_end 10^(-10) 10^-2]);
% ylim([0 2])
% title('Mean Trajectories');
xlabel('Time (hours)'); ylabel({'Mean Induction Rate, $\langle \gamma \rangle$'},'interpreter','latex');
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

% legend(r,{'$\bar{\varepsilon}$','$\bar{\delta}_I$'},'Interpreter','Latex','Location','NorthEast');

subplot(3,1,3);
plot(params.t_span, mu_q_L_traj,'-','Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
% plot(params.t_span, Rt_SIR_traj,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;
% q(2) = plot(params.t_span, Rt_SIR_iv,'--','Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
% q(3) = plot(params.t_span, Rt_SIR_v,':','Color',my_rgb_colors(3,:),'LineWidth',2); hold on;
% axis([0 t_end 0 4]);
% ylim([0 4]);

xlabel('Time (hours)'); ylabel({'Mean Probability Lysogeny, $\langle q \rangle$'},'interpreter','latex');
% legend(q,{'Correlated','independent','classic'},'Location','SouthEast');
%'correlated',
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');


%%
if 1
    
    %% Movies
    % create this video
%     f3= figure(3);
%     set(f3, 'Position', [100   450   840   700]);
    
    
    
    %     name_videofile = strcat('joint_S_and_I_dist');
    %     v1 = VideoWriter(name_videofile,'MPEG-4');
    %     v1.FrameRate = 1;
    
    %     open(v1);
    
    %     this_ind = 2;
    for this_ind = 1:100:length(t_span)
        
        fig = figure('Visible','on');
        set(fig, 'Position', [100   450   840   700]);
        
        subplot(2,2,1);
        this_p(1)=semilogy(params.t_span, S_traj,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;
        this_p(2)=semilogy(params.t_span, I_traj,'Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
        this_p(3)=semilogy(params.t_span, L_traj,'Color',my_rgb_colors(3,:),'LineWidth',2); hold on;
        this_p(4)=semilogy(params.t_span, V_traj,'k','LineWidth',2); hold on;
        plot(t_span(this_ind)*ones(size(t_span)),linspace(10^0,10^10,length(t_span)),'k--','LineWidth',1); hold on;

        axis([0 t_end 10^0 10^10]);
        xlabel('Time (hours)'); ylabel({'Population Size'});
        title('Dynamics')
        legend(this_p,{'S','I','L', 'V'},'Location','SouthEast');
        set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

        
%         plot(t_span,mu_gamma_L_traj,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;
%         
%         
%         xlabel('Time (hours)');
%         ylabel({'Mean Susceptibility, $\bar{\varepsilon}$'},'interpreter','latex');
%         axis([0 t_end 0 1.2]);
%         set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');
        
        %         f=gca;
        %         f.LineWidth = 1;
        %         f.FontSize = 14;
        %         f.FontWeight = 'normal';
        %         f.FontName = 'Times';
        
        
        subplot(2,2,2);
        this_q(1)=plot(t_span,mu_gamma_L_traj,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;
        this_q(2)=plot(t_span,mu_q_L_traj,'Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
        plot(t_span(this_ind)*ones(size(t_span)),linspace(0,2,length(t_span)),'k--','LineWidth',1); hold on;
        
        xlabel('Time (hours)');
        ylabel({'Means of Lysogen traits'},'interpreter','latex');
        axis([0 t_end 0 0.05]);
        set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');
        %         f=gca;
        %         f.LineWidth = 1;
        %         f.FontSize = 14;
        %         f.FontWeight = 'normal';
        %         f.FontName = 'Times';
        legend(this_q,{'Induction rate, $\langle \gamma \rangle$','Lysogeny Probability, $\langle q \rangle$'},'interpreter','latex','Location','SouthEast');
        
        subplot(2,2,3);
        this_joint_S = reshape(S_traj_discrete(this_ind,:,:)/S_traj(this_ind),m,n);
        imagesc(gamma_vector,q_vector_plt,this_joint_S);
        % pcolor(eps,del,joint_S);
        
        set(gca,'YDir','normal');
        colorbar;
        xlim([0 0.1]); ylim([0 0.1]);
        caxis([0 0.06]);
        
        xlabel('Induction rate, $\gamma$','interpreter','latex');
        ylabel({'Lysogeny Probability, $q$'},'interpreter','latex');
        set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');
        %         f=gca;
        %         f.LineWidth = 1;
        %         f.FontSize = 14;
        %         f.FontWeight = 'normal';
        %         f.FontName = 'Times';
        
        title('Joint Distribution in $S$','interpreter','latex');
        
        subplot(2,2,4);
        this_joint_L = reshape(L_traj_discrete(this_ind,:,:)/L_traj(this_ind),m,n);
        imagesc(gamma_vector,q_vector_plt,this_joint_L);
        % pcolor(eps,del,joint_S);
        %axis xy;
        set(gca,'YDir','normal');
        colorbar;
        xlim([0 0.1]); ylim([0 0.1]);
        % clim([]);
        caxis([0 0.07]);
        
        xlabel('Induction rate, $\gamma$','interpreter','latex');
        ylabel({'Lysogeny Probability, $q$'},'interpreter','latex');
        set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');
        %         f=gca;
        %         f.LineWidth = 1;
        %         f.FontSize = 14;
        %         f.FontWeight = 'normal';
        %         f.FontName = 'Times';
        
        title('Joint Distribution in $L$','interpreter','latex');
        
        
        %         frame = getframe(fig);
        %         writeVideo(v1,frame);
        %         close(fig);
        
        frame = getframe(fig);
        im{this_ind} = frame2im(frame);
        [imind,cm] = rgb2ind(im{this_ind},256);
        % Write to the GIF File
        
        if write_video_ans ==1
            if this_ind == 1
                imwrite(imind,cm,video_filename,'gif','DelayTime',.1, 'Loopcount',inf);
            else
                imwrite(imind,cm,video_filename,'gif','DelayTime',.1, 'WriteMode','append');
            end
        end
        
    end
    
    %     fig = figure('Visible','off');
    
    
end


%% Save results

results.S_traj = S_traj;
results.S_traj_discrete = S_traj_discrete;
% results.S_traj_v = S_traj_v;

results.I_traj = I_traj;
results.I_traj_discrete = E_traj_discrete;
% results.I_traj_v = I_traj_v;

% results.R_traj = R_traj;
% results.R_traj_discrete = I_traj_discrete;
% results.R_traj_v = R_traj_v(1:s:end);

% results.Rt_SIR_traj = Rt_SIR_traj;
% results.Rt_SIR_v = Rt_SIR_v(1:s:end);

% results.FOS_traj = FOS_traj;
% results.FOS_SIR_traj_v = FOS_SIR_traj_v;

% marginals over time
% results.marg_eps_I_traj = marginal_eps_I_traj;
% results.marg_delta_I_traj = marginal_delta_I_traj;
% results.marg_delta_S_traj = marginal_delta_S_traj;
% results.marg_eps_S_traj = marginal_gamma_L_traj;

% mean value trajectories
results.mu_epsilon_S_traj = mu_gamma_L_traj;
results.mu_delta_I_traj = mu_q_L_traj;

% results.mu_delta_S = mu_delta_S;
% results.mu_eps_I = mu_eps_I;
% results.mu_delta_I = mu_delta_I;

% variance
results.var_eps_S_traj = variance_gamma_L_traj;
results.var_delta_I_traj = variance_q_L_traj;
% results.dispSe_traj = dispSe_traj;
% results.dispId_traj = dispId_traj;

results.joint_s = init_joint_S;
% results.joint_i = init_joint_I;

% results.marg_eps_s = init_marginal_gamma_S;
% results.marg_delta_S=init_marginal_q_S;
% results.marg_eps_I=init_marg_eps_I;
% results.marg_delta_I=init_marg_delta_I;


%% homogeneous
% results.FOS_traj_v = FOS_SIR_traj_v;

%%
% save simulated results
if save_ans==1
    
    folder_location = './sim_results/';
    save(strcat(folder_location,filename),'params','results');
    
    fprintf('Saved to file: \n');
    fprintf(strcat(filename,'\n'));
    
else
    
    fprintf('Not Saved. \n');
    
end

