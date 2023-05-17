% simulate SIR model with transmissibility & susceptibility variation
% discretized version with correlation
%using negative binomial distributions

%%

clear all; close all; clc;

%% want to save?
save_ans = 0;
% 0: don't save
% 1: save

% filename = 'negbinomial_Poissonlike.mat';
filename = 'negbinomial_overdispersed_uncorr.mat';
% filename = 'negbinomial_overdispersed_poscorr.mat';

default_colors = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];

%% parameters
gam=1/10;

% dispersion parameters - kappa->infty = Poisson; kappa=1 = Geometric;
% kappa<1 overdispersed
kappa1 = 50; % dispersion (susceptibility)
kappa2 = 0.5; % dispersion (transmissibility)

% For the Singapore outbreak, the maximum-likelihood estimate k̂ is 0.16
% (90% confidence interval 0.11–0.64), indicating highly overdispersed

% correlation coefficient for initial joint in S population
rho_S = 0; % 0, 0.6 

intended_R0 = 3;
% match basic reproduction numbers: R0=3
bet=intended_R0*gam; %Poisson-like, rho=0
% bet = intended_R0*gam*.5; % overdispersed, rho=0.6

mean_joint = 1;

N = 1; % starting population;

% discrete mesh
n = 100; % number transmissibility classes
m = n; % susceptibilty classes

eps = 0:1:m;
del = 0:1:n;

eps_plt = eps;
eps(1)=0.05;
% del(1)=0.1;

[Eps,Del] = meshgrid(eps,del);
% [Del,Eps] = meshgrid(del,eps);
% Eps=Eps';
% Del=Del';

params.beta = bet;
params.gamma = gam;
params.mean = mean_joint;
params.N = N;
% params
params.eps=eps;
params.del=del;
params.E=Eps;
params.D=Del;


t_start = 0; t_end = 300; % 200 days is about 6-7 months; 250 is 8-9 months
t_span = t_start:1:t_end; %0.2
params.t_span = t_span;


% Initialize Joint Distributions & calculate Marginal
% x1 = eps, x2 = delta

m1 = 1/kappa1; % reciprocal of dispersion (susceptibility)
m2 = 1/kappa2; % reciprocal of dispersion (transmissibility)


%theta values - theta_t = mt*mean_joint/(mt*mean_joint +1)
% want means of 1
mu1 = 1;
mu2 = mu1;
theta1 = m1*mu1/(m1*mu2+1);
theta2 = m2*mu2/(m2*mu2+1);
% rearranging: means of marginals
% mu1 = theta1/(1-theta1)/m1;
% mu2 = theta2/(1-theta2)/m2;

var1 = theta1/(1-theta1)^2/m1; % variances of marginals
var2 = theta2/(1-theta2)^2/m2;

c1 = ((1-theta1)/(1-theta1*exp(-1)))^(1/m1);
c2 = ((1-theta2)/(1-theta2*exp(-1)))^(1/m2);

A1 = m1^(-1)*theta1*exp(-1)/(1-theta1*exp(-1))-m1^(-1)*theta1/(1-theta1);
A2 = m2^(-1)*theta2*exp(-1)/(1-theta2*exp(-1))-m2^(-1)*theta2/(1-theta2);

lambda_S = rho_S*(sqrt(var1)*sqrt(var2))/(c1*c2*A1*A2); % calculate lambda - multiplicative factor

%bivariate probability function for S compartment
joint_S = biNegBinomialPdf(Eps,Del,theta1,theta2,m1,m2,lambda_S);
% joint_S = biNegBinomialPdf(Del,Eps,theta2,theta1,m2,m1,lambda_S);
joint_I = joint_S;

%normalize here
joint_S = joint_S./sum(sum(joint_S));
joint_I = joint_I./sum(sum(joint_I));

marg_eps_S = sum(joint_S);
marg_delta_S = sum(joint_S,2)';

marg_eps_I = sum(joint_I);
marg_delta_I = sum(joint_I,2)';

% calculated mean values
mu_eps_S = sum(eps.*marg_eps_S);
mu_delta_S = sum(del.*marg_delta_S);
mu_eps_I = sum(eps.*marg_eps_I);
mu_delta_I = sum(del.*marg_delta_I);



%calculated values
cov = sum(sum(((eps' - mu_eps_S*ones(size(eps')))*(del- mu_delta_S*ones(size(del)))).*joint_S));

var_eps = sum((eps - mu_eps_S*ones(size(eps))).^2.*marg_eps_S);
sd_eps = sqrt(var_eps);

var_delta = sum((del- mu_delta_S*ones(size(del))).^2.*marg_delta_S);
sd_delta = sqrt(var_delta);

corrcoef = cov/(sd_eps*sd_delta);

% plot initial distributions
f1=figure(1); set(f1, 'Position', [10   50   840   350]);
subplot(1,2,1);

imagesc(eps_plt,del,joint_S);
% pcolor(eps,del,joint_S);
%axis xy;
set(gca,'YDir','normal');
colorbar;
xlim([-0.5 4.5]); ylim([-0.5 4.5]);

xlabel('susceptibility $\varepsilon$','interpreter','latex');
ylabel({'transmissibility $\delta$'},'interpreter','latex');
%     zlabel('Population Fraction')
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

title('Joint Distribution in S');

f1=figure(1);
subplot(1,2,2);
q(1)=plot(eps_plt,marg_eps_S,'.-','Color',default_colors(3,:),'LineWidth',2.5,'MarkerSize',20); hold on;
q(2)=plot(del, marg_delta_S,'.-','Color',default_colors(2,:),'LineWidth',2.5,'MarkerSize',20); hold on;
% cat_margs = [marg_eps_S;marg_delta_S];
% q=bar(eps,cat_margs);
axis([-0.5 4.5 0 1]);

axis([0 4 0 1]);
title('Marginal Distributions');
xlabel('transmissibility d');

xlabel('susceptibility $\varepsilon$ or transmissibility $\delta$','interpreter','latex');
ylabel({'Population Fraction'});


f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

legend(q,{'Susceptibility','Potential Transmissibility'},'Location','NorthEast');

%

% %change epsilon as needed
eps_perturb = 5e-8;
%
% %reset m & n
% % m = length(marg_delta_S);
% % n = length(marg_eps_S);
%
params.eps_perturb = eps_perturb;

params.mu_eps_S = mu_eps_S;
params.mu_delta_I = mu_delta_I;

eigen_direction_SIR = get_eigendirection_SIRdelta_eps(params);

S_init = 1 - eps_perturb*abs(eigen_direction_SIR(1));
I_init = eps_perturb*abs(eigen_direction_SIR(2));
R_init = eps_perturb*abs(eigen_direction_SIR(3));

% S_init = 1;
% I_init = 0;
% R_init = 0;


%
init_dist_S = joint_S.*(S_init); %joint probability matrix S_{i,j}
params.init_dist_S = init_dist_S;


init_dist_I = joint_I.*I_init;
params.init_dist_I = init_dist_I;

init_dist_R = joint_S.*R_init;

%check should equal 1
%sum(sum((init_dist_S + init_dist_I + init_dist_R),2));

init_conds = [reshape(init_dist_S, (m+1)*(n+1),1); reshape(init_dist_I,(m+1)*(n+1),1); reshape(init_dist_R, (m+1)*(n+1),1)];

% Simulate model
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

S_traj_discrete =zeros(length(t_span),(m+1),(n+1));
I_traj_discrete = zeros(length(t_span),(m+1),(n+1));
R_traj_discrete = zeros(length(t_span),(m+1),(n+1));

[t,y_traj] = ode45(@(t,y)simulate_SIRdelta_eps(t,y,params), params.t_span, init_conds, options);

S_traj_discrete = reshape(y_traj(:,1:(m+1)*(n+1)), length(t_span),(m+1),(n+1));
I_traj_discrete = reshape(y_traj(:,((m+1)*(n+1)+1):(2*(m+1)*(n+1))),length(t_span),(m+1),(n+1));
R_traj_discrete = reshape(y_traj(:,(2*(m+1)*(n+1)+1):(3*(m+1)*(n+1))),length(t_span),(m+1),(n+1));



% S,I,R & Marginal eps & delta Distribution Trajectories
S_traj = sum(S_traj_discrete,[2,3]); % population size
I_traj = sum(I_traj_discrete,[2,3]);
R_traj = sum(R_traj_discrete,[2,3]);

marg_eps_S_traj = reshape(sum(S_traj_discrete,2),length(t_span),m+1);
marg_delta_S_traj = reshape(sum(S_traj_discrete,3),length(t_span),n+1);
marg_eps_I_traj = reshape(sum(I_traj_discrete,2),length(t_span),m+1);
marg_delta_I_traj = reshape(sum(I_traj_discrete,3),length(t_span),n+1);

for kk=1:length(t_span)
    
    mu_epsilon_S_traj(kk) = sum(eps.*marg_eps_S_traj(kk,:))./S_traj(kk);
    mu_delta_I_traj(kk) = sum(del.*marg_delta_I_traj(kk,:))./I_traj(kk);
    
end


SIRed_traj = zeros(length(t_span),5);
SIRed_traj(:,1) = S_traj;
SIRed_traj(:,2) = I_traj;
SIRed_traj(:,3) = R_traj;
SIRed_traj(:,4) = mu_epsilon_S_traj;
SIRed_traj(:,5) = mu_delta_I_traj;

Rt_SIR_traj = get_Rt_SIR_delta_eps(params,SIRed_traj);

for kk = 1:length(t_span)
    
    marg_eps_S_traj(kk,:) = marg_eps_S_traj(kk,:)/sum(marg_eps_S_traj(kk,:));
    marg_delta_S_traj(kk,:) = marg_delta_S_traj(kk,:)/sum(marg_delta_S_traj(kk,:));
    
    marg_eps_I_traj(kk,:) = marg_eps_I_traj(kk,:)/sum(marg_eps_I_traj(kk,:));
    marg_delta_I_traj(kk,:) = marg_delta_I_traj(kk,:)/sum(marg_delta_I_traj(kk,:));
    
    
    var_eps_S_traj(kk) = sum((eps- mu_epsilon_S_traj(kk)*ones(size(eps))).^2.*marg_eps_S_traj(kk,:));
    var_delta_I_traj(kk) = sum(((del- mu_delta_I_traj(kk)*ones(size(del))).^2).*marg_delta_I_traj(kk,:));
    
    std_Se_traj(kk) = sqrt(var_eps_S_traj(kk));
    std_Id_traj(kk) = sqrt(var_delta_I_traj(kk));
    
end

%dispersion parameters
% dispSe_traj = std_Se_traj./mu_epsilon_S_traj;
% dispId_traj = std_Id_traj./mu_delta_I_traj;

CV_Se_traj = var_eps_S_traj./(mu_epsilon_S_traj.^2);
CV_Id_traj = var_delta_I_traj./(mu_delta_I_traj.^2);


% Calculate final size of outbreak
% FOSed_traj = getFinalOutbreakSize(I_traj, R_traj);
FOS_traj = I_traj + R_traj;

% Simulate classic SIR model
% [S_traj_v, I_traj_v, R_traj_v, Rt_SIR_v, FOS_SIR_traj_v, params_v] = simulateClassic();

% simulate Independent var model
% [S_traj_iv, I_traj_iv, R_traj_iv, Rt_SIR_iv, FOSiv_traj,params_iv] = simulateINDvarnegbinomial();


% Plotting
f2=figure(2); set(f2, 'Position', [900   50   400   930]);

subplot(3,1,1);
q(1)=plot(params.t_span, S_traj,'Color',default_colors(1,:),'LineWidth',2); hold on;
q(2)=plot(params.t_span, I_traj,'Color',default_colors(2,:),'LineWidth',2); hold on;
q(3)=plot(params.t_span, R_traj,'Color',default_colors(3,:),'LineWidth',2); hold on;
% q(4)=semilogy(params.t_span, S_traj_iv,'--','Color',default_colors(1,:),'LineWidth',2); hold on;
% q(5)=semilogy(params.t_span, I_traj_iv,'--','Color',default_colors(2,:),'LineWidth',2); hold on;
% q(6)=semilogy(params.t_span, R_traj_iv,'--','Color',default_colors(3,:),'LineWidth',2); hold on;
% q(7)=semilogy(params.t_span, S_traj_v,':','Color',default_colors(1,:),'LineWidth',2); hold on;
% q(8)=semilogy(params.t_span, I_traj_v,':','Color',default_colors(2,:),'LineWidth',2); hold on;
% q(9)=semilogy(params.t_span, R_traj_v,':','Color',default_colors(3,:),'LineWidth',2); hold on;

axis([0 t_end 10^-5 0.1*10]);
xlabel('Time (days)'); ylabel({'Population'; 'Fraction'});
title('Dynamics')
legend(q,{'S','I','R'},'Location','SouthEast');
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

subplot(3,1,2);
r(1)=plot(params.t_span, mu_epsilon_S_traj,'-','Color',default_colors(2,:),'LineWidth',2); hold on;
r(2)=plot(params.t_span, params.beta*mu_delta_I_traj,'-','Color',default_colors(3,:),'LineWidth',2); hold on;

axis([0 t_end 0 1.1]);
% ylim([0 2])
title('Mean Trajectories');
xlabel('Time (days)'); ylabel({'Population'; 'Fraction'});
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

legend(r,{'$\bar{\varepsilon}$','$\bar{\delta}_I$'},'Interpreter','Latex','Location','SouthEast');

subplot(3,1,3);
plot(params.t_span, Rt_SIR_traj,'Color',default_colors(1,:),'LineWidth',2); hold on;
% q(2) = plot(params.t_span, Rt_SIR_iv,'--','Color',default_colors(2,:),'LineWidth',2); hold on;
% q(3) = plot(params.t_span, Rt_SIR_v,':','Color',default_colors(3,:),'LineWidth',2); hold on;
axis([0 t_end 0 4]);
% ylim([0 4]);

xlabel('Time (days)'); ylabel([{'Effective'; 'Reproduction'; 'Number'}]);
% legend(q,{'Correlated','independent','classic'},'Location','SouthEast');
%'correlated',
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');


% plot joint in S and I at certain time point
f3=figure(3); set(f3, 'Position', [100   450   840   350]);

this_ind = 50;
subplot(1,2,1);
this_joint_S = reshape(S_traj_discrete(this_ind,:,:),(m+1),(n+1))/S_traj(this_ind);
imagesc(eps_plt,del,this_joint_S);
% pcolor(eps,del,joint_S);

set(gca,'YDir','normal');
colorbar;
xlim([-0.5 6.5]); ylim([-0.5 6.5]);
caxis([0 0.35]);

xlabel('susceptibility $\varepsilon$','interpreter','latex');
ylabel({'transmissibility $\delta$'},'interpreter','latex');

f2=gca;
f2.LineWidth = 1;
f2.FontSize = 14;
f2.FontWeight = 'normal';
f2.FontName = 'Times';

subplot(1,2,2);
this_joint_I = reshape(I_traj_discrete(this_ind,:,:),(m+1),(n+1))/I_traj(this_ind);
imagesc(eps_plt,del,this_joint_I);
% pcolor(eps,del,joint_S);
%axis xy;
set(gca,'YDir','normal');
colorbar;
xlim([-0.5 6.5]); ylim([-0.5 6.5]);
% clim([]);
caxis([0 0.35]);

xlabel('susceptibility $\varepsilon$','interpreter','latex');
ylabel({'transmissibility $\delta$'},'interpreter','latex');

f2=gca;
f2.LineWidth = 1;
f2.FontSize = 14;
f2.FontWeight = 'normal';
f2.FontName = 'Times';



%% Save results 

results.S_traj = S_traj;
results.S_traj_discrete = S_traj_discrete;
% results.S_traj_v = S_traj_v;

results.I_traj = I_traj;
results.I_traj_discrete = I_traj_discrete;
% results.I_traj_v = I_traj_v;

results.R_traj = R_traj;
results.R_traj_discrete = R_traj_discrete;
% results.R_traj_v = R_traj_v(1:s:end);

results.Rt_SIR_traj = Rt_SIR_traj;
% results.Rt_SIR_v = Rt_SIR_v(1:s:end);

results.FOS_traj = FOS_traj;
% results.FOS_SIR_traj_v = FOS_SIR_traj_v;

% marginals over time
results.marg_eps_I_traj = marg_eps_I_traj;
results.marg_delta_I_traj = marg_delta_I_traj;
results.marg_delta_S_traj = marg_delta_S_traj;
results.marg_eps_S_traj = marg_eps_S_traj;

% mean value trajectories
results.mu_epsilon_S_traj = mu_epsilon_S_traj;
results.mu_delta_I_traj = mu_delta_I_traj;

% results.mu_delta_S = mu_delta_S;
% results.mu_eps_I = mu_eps_I;
% results.mu_delta_I = mu_delta_I;

% variance 
results.var_eps_S_traj = var_eps_S_traj;
results.var_delta_I_traj = var_delta_I_traj;
% results.dispSe_traj = dispSe_traj;
% results.dispId_traj = dispId_traj;

results.joint_s = joint_S;
results.joint_i = joint_I;

results.marg_eps_s = marg_eps_S;
results.marg_delta_S=marg_delta_S;
results.marg_eps_I=marg_eps_I;
results.marg_delta_I=marg_delta_I;


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

