function [S_traj_iv, I_traj_iv, R_traj_iv, Rt_SIR_iv, FOSiv_traj, params] = simulateINDvarnegbinomial()
    default_colors = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];  

    gam=1/10;

    %change c0 as needed to match growth rates 
    c0 = 0; 
    bet=2*gam;%2*gam; 


    N = 1; %starting population; 

    %discretely valued mesh

    n = 20; % n different transmissibility classes
    m = n; % m different susceptibilty classes

    d = 0:1:n;
    e = 0:1:m;

    [E,D] = meshgrid(e,d);

    % D = d.*ones(n+1,n+1); % same across columns
    % E = e'.*ones(m+1,m+1); % same across rows
    % 

    params.beta = bet;
    params.gamma = gam;
    params.c0 = c0;
    params.mean = 1; %change later ******************

    params.N = N;

    %params.end_val = end_val;

    t_start = 0; t_end = 400; % 200 days is about 6-7 months; 250 is 8-9 months
    t_span = t_start:0.2:t_end;
    params.t_span = t_span;




    %% Initialize Joint Distributions & calculate Marginal

    %manually input parameters

    %dispersion parameters
    m1 = 1;
    m2 = 1/20;

    %theta values - tt = mt*mean/(mt*mean +1)
    theta1 = 1/2;
    theta2 = 1/21; %0.099;

    mu1 =m1^(-1)*theta1/(1-theta1); %what means theoretically should be
    mu2 = m2^(-1)*theta2/(1-theta2);

    sig1 =m1^(-1)*theta1/(1-theta1)^2; %what variances theoretically should be
    sig2 = m2^(-1)*theta2/(1-theta2)^2;

    c1 = ((1-theta1)/(1-theta1*exp(-1)))^(m1^(-1));
    c2 = ((1-theta2)/(1-theta2*exp(-1)))^(m2^(-1));

    A1 = m1^(-1)*theta1*exp(-1)/(1-theta1*exp(-1))-m1^(-1)*theta1/(1-theta1);
    A2 = m2^(-1)*theta2*exp(-1)/(1-theta2*exp(-1))-m2^(-1)*theta2/(1-theta2);

    rho_S = 0.0; %input correlation coefficient
    lambda_S = rho_S*(sqrt(sig1)*sqrt(sig2))/(c1*c2*A1*A2); %calculate lambda

    %probability function for S compartment
    pdf_S = biNegBinomialPdf(E,D,theta1,theta2,m1,m2,lambda_S);


    theta1 = 1/2; %mean should theoretically be 10 in this case
    theta2 = 1/21; %0.1238;

    mu1 =m1^(-1)*theta1/(1-theta1); %what means theoretically should be
    mu2 = m2^(-1)*theta2/(1-theta2);

    sig1 =m1^(-1)*theta1/(1-theta1)^2; %what variances theoretically should be
    sig2 = m2^(-1)*theta2/(1-theta2)^2;

    c1 = ((1-theta1)/(1-theta1*exp(-1)))^(m1^(-1));
    c2 = ((1-theta2)/(1-theta2*exp(-1)))^(m2^(-1));

    A1 = m1^(-1)*theta1*exp(-1)/(1-theta1*exp(-1))-m1^(-1)*theta1/(1-theta1);
    A2 = m2^(-1)*theta2*exp(-1)/(1-theta2*exp(-1))-m2^(-1)*theta2/(1-theta2);

    rho_I = 0.0;
    lambda_I = rho_I*(sqrt(sig1)*sqrt(sig2))/(c1*c2*A1*A2);

    %pdf for I compartment 
    pdf_I =biNegBinomialPdf(E,D,theta1,theta2,m1,m2,lambda_I) ;

    %plot joint S
%     imagesc(e,d,pdf_S');
%     set(gca,'YDir','normal');

    %interpolate here for a finer mesh
    % 
    end_val = n;
    % d_interp = 0:end_val/200:end_val;
    % e_interp = 0:end_val/200:end_val;
    % 
    % [E_interp,D_interp] = meshgrid(e_interp,d_interp);
    % 
    % joint_s = interp2(E,D,pdf_S,E_interp,D_interp,'spline');
    % 
    % joint_i = interp2(E,D,pdf_I,E_interp,D_interp,'spline');
    % 


    % % 
    d_interp = d;
    e_interp = e;
    E_interp = E;
    D_interp = D;
    joint_s = pdf_S;
    joint_i = pdf_I;

    D= D_interp';
    E = E_interp';
    params.D = D;
    params.E = E;
    d = d_interp;
    e = e_interp;
    params.end_val = end_val;


    %normalzie here
    joint_s = joint_s./sum(sum(joint_s),2); 
    joint_i = joint_i./sum(sum(joint_i),2);

    %sum(sum(joint_s),2) = 1

    marg_d_s = sum(joint_s); 
    marg_e_s = sum(joint_s,2)'; 

    marg_d_i = sum(joint_i);
    marg_e_i = sum(joint_i,2)'; %check if marginals are correctly labeled!!


    %check the means
    mu_delta_S = sum(d.*marg_d_s)/sum(marg_d_s); % close to 1, but not 1
    mu_eps_S = sum(e.*marg_e_s)/sum(marg_e_s);
    mu_delta_I = sum(d.*marg_d_i)/sum(marg_d_i);
    mu_eps_I = sum(e.*marg_e_i)/sum(marg_e_i);

    params.mu_eps_S = mu_eps_S;
    params.mu_delta_I = mu_delta_I;
    params.mu_delta_S = mu_delta_S;
    params.mu_eps_I = mu_eps_I;
    params.joint_s = joint_s;
    params.joint_i = joint_i;
    params.marg_e_s = marg_e_s;
    

    %% plotting initial distributions
    

%     plotInitDists(marg_e_s, marg_d_s, marg_e_i, marg_d_i, e,d, joint_s, joint_i,E,D, default_colors);
%     close;
    %reset m & n
    m = length(marg_d_s);
    n = length(marg_e_s);

    %check
    % marg_d_s = sum(joint_s);
    % marg_e_s = sum(joint_s,2)';
    % marg_d_i = sum(joint_i);
    % marg_e_i = sum(joint_i,2)';


    eps = 1e-6;
    params.eps = eps;

    
    [eigdir_S, eigdir_I, eigdir_R,eigen_direction_SIR_ed] = get_eigendirection_SIRdelta_eps(params);


    S_init = 1 + eps*eigen_direction_SIR_ed(1);
    I_init = eps*eigen_direction_SIR_ed(2);
    R_init = eps*eigen_direction_SIR_ed(3);

    init_dist_S = joint_s*(S_init); %joint probability matrix S_{i,j}
    params.init_dist_S = init_dist_S;

    init_dist_I = joint_i*I_init; 
    params.init_dist_I = init_dist_I;

    init_dist_R = joint_s*R_init;



    init_conds = [reshape(init_dist_S, (m)*(n),1); reshape(init_dist_I,(m)*(n),1); reshape(init_dist_R, (m)*(n),1)];



    %% Simulate model

    options = odeset('RelTol',1e-12,'AbsTol',1e-12);

    S_traj_discrete =zeros(length(t_span),m,n);
    I_traj_discrete = zeros(length(t_span),m,n);
    R_traj_discrete = zeros(length(t_span),m,n);

    %find some way to matrix ode solve - has to be done simulataneously - maybe
    %vectorize to mn vector instead of a mxn matrix

    [t,y_traj] = ode45(@(t,y)simulate_SIRdelta_eps(t,y,params), params.t_span, init_conds, options);

    S_traj_discrete(:,:,:) = abs(reshape(y_traj(:,1:m*n), length(t_span),m,n));
    I_traj_discrete(:,:,:) = abs(reshape(y_traj(:,m*n+1:2*m*n),length(t_span),m,n));
    R_traj_discrete(:,:,:) = abs(reshape(y_traj(:,2*m*n+1:3*m*n),length(t_span),m,n));

    %taking abs so that everything is positive - doesn't help though

    %% S,I,R & Marginal eps & delta Distribution Trajectories

    S_traj = sum(S_traj_discrete,[2,3]);
    I_traj = sum(I_traj_discrete,[2,3]);
    R_traj = sum(R_traj_discrete,[2,3]);

    marg_S_e_traj = reshape(sum(S_traj_discrete,3),length(t_span),m); 
    marg_S_d_traj = reshape(sum(S_traj_discrete,2),length(t_span),n);
    marg_I_e_traj = reshape(sum(I_traj_discrete,3),length(t_span),m);
    marg_I_d_traj = reshape(sum(I_traj_discrete,2),length(t_span),n);

    delta_S_dist_traj = sum(d.*marg_S_d_traj,2)./S_traj;
    mu_delta_S_traj = sum(delta_S_dist_traj,2);

    epsilon_S_dist_traj = sum(e.*marg_S_e_traj,2)./S_traj;
    mu_epsilon_S_traj = sum(epsilon_S_dist_traj,2);

    delta_I_dist_traj = sum(d.*marg_I_d_traj,2)./I_traj; %fishy
    mu_delta_I_traj = sum(delta_I_dist_traj,2);

    epsilon_I_dist_traj = sum(e.*marg_I_e_traj,2)./I_traj; %fishy
    mu_epsilon_I_traj = sum(epsilon_I_dist_traj,2);

    SIRed_traj = zeros(length(t_span),5);
    SIRed_traj(:,1) = S_traj;
    SIRed_traj(:,2) = I_traj;
    SIRed_traj(:,3) = R_traj;
    SIRed_traj(:,4) = mu_epsilon_S_traj;
    SIRed_traj(:,5) = mu_delta_I_traj;

    Rt_SIR_iv = get_Rt_SIR_delta_eps(params,SIRed_traj);
    FOSiv_traj = getFinalOutbreakSize(I_traj, R_traj);

    for tt = 1:length(t_span)
        marg_S_e_traj(tt,:) = marg_S_e_traj(tt,:)/sum(marg_S_e_traj(tt,:));
        marg_I_d_traj(tt,:) = marg_I_d_traj(tt,:)/sum(marg_I_d_traj(tt,:));
        varSe_traj(tt) = sum((E(:,1)'- mu_epsilon_S_traj(tt)*ones(size(E(:,1)'))).^2.*marg_S_e_traj(tt,:));
        varId_traj(tt) = sum((D(1,:)- mu_delta_I_traj(tt)*ones(size(D(1,:)))).^2.*marg_I_d_traj(tt,:),2);
    end
    dispSe_traj = mu_epsilon_S_traj.^2./(varSe_traj'-mu_epsilon_S_traj);
    dispId_traj = mu_delta_I_traj.^2./(varId_traj'-mu_delta_I_traj);
    


    %% Joint Distribution Trajectories for Sus. & Infected

    joint_s_traj = zeros(length(t_span),m,n);
    joint_i_traj = zeros(length(t_span),m,n);
    for t = 1:length(t_span)
        joint_s_traj(t,:,:) = marg_S_e_traj(t,:,:)'*marg_S_d_traj(t,:,:); 
        joint_i_traj(t,:,:) = marg_I_e_traj(t,:,:)'*marg_I_d_traj(t,:,:);
    end

        %another plot
    joint_s_new = zeros(m,n);
    joint_i_new = joint_s_new;
    joint_s_new(:,:) = joint_s_traj(1,:,:);
    joint_i_new(:,:) = joint_i_traj(1,:,:);
%     plotInitDists(marg_S_e_traj(1,:), marg_S_d_traj(1,:), marg_I_e_traj(1,:), marg_I_d_traj(1,:), E(:,1),D(1,:), joint_s_new, joint_i_new,E,D, default_colors);
%     close;

    S_traj_iv = S_traj;
    I_traj_iv = I_traj;
    R_traj_iv = R_traj;
    
    %% simulate classic
    [S_traj_v, I_traj_v, R_traj_v, Rt_SIR_v, FOS_SIR_traj, params_v] = simulateClassic();


    %% Plotting
    
%     figure(31)
%     q(1)=semilogy(params.t_span, S_traj_iv,'--','Color',default_colors(1,:),'LineWidth',2); hold on;
%     q(2)=semilogy(params.t_span, I_traj_iv,'--','Color',default_colors(2,:),'LineWidth',2); hold on;
%     q(3)=semilogy(params.t_span, R_traj_iv,'--','Color',default_colors(3,:),'LineWidth',2); hold on;
%     q(4) = semilogy(params.t_span, S_traj_v,':','Color',default_colors(1,:),'LineWidth',2); hold on;
%     q(5) = semilogy(params.t_span, I_traj_v,':','Color',default_colors(2,:),'LineWidth',2); hold on;
%     q(6) = semilogy(params.t_span, R_traj_v,':','Color',default_colors(3,:),'LineWidth',2); hold on;
%     axis([0 t_end 10^(-6) 10]);
%     xlabel('Time (days)'); ylabel({'Population'; 'Fraction'});
% 
%     legend(q,{'S-iv','I-iv','R-iv'},'Location','SouthEast');
%     
%     figure(32)
%     plot(params.t_span, Rt_SIR_iv,'--','Color',default_colors(2,:),'LineWidth',2);
%     hold on;
%     plot(params.t_span, Rt_SIR_v,':','Color',default_colors(3,:),'LineWidth',2);
%     xlabel('Time (days)'); ylabel([{'Effective'; 'Reproduction'; 'Number'}]);
%     axis([0 t_end 10^(-6) 10]);
%     ylim([0 4]);
%     
%     figure(33);
%     set(gca,'FontSize', 14, 'LineWidth',1,'FontWeight','normal','FontName','Times');
%     q(1)=plot(params.t_span, 100.*FOSiv_traj,'--','Color',default_colors(2,:),'LineWidth',2); hold on;
%     q(2)=plot(params.t_span, 100.*FOS_SIR_traj,':','Color',default_colors(3,:),'LineWidth',2); hold on;
%     axis([0 t_end 10^(-6) 10]);
%     ylim([0,105]);
%     xlabel('Time (days)'); ylabel({'Cumulated Infected (%)'});
%     legend(q,{'Independent', 'Classic'},'Location', 'SouthEast');
%     title('Cumulative Infected (I + R)');

    %% saving variables
    
    s = 20;
    params.s = s;
    params.t_span = t_span(1:s:end);
data.FOS_SIR_traj = FOS_SIR_traj(1:s:end);
data.FOSiv_traj = FOSiv_traj(1:s:end);
data.I_traj_discrete = I_traj_discrete(1:s:end,:,:);
data.I_traj = I_traj_iv(1:s:end);
data.I_traj_v = I_traj_v(1:s:end);
data.marg_I_e_traj = marg_I_e_traj(1:s:end,:);
data.marg_I_d_traj = marg_I_d_traj(1:s:end,:);
data.marg_S_d_traj = marg_S_d_traj(1:s:end,:);
data.marg_S_e_traj = marg_S_e_traj(1:s:end,:);
data.mu_delta_I_traj = mu_delta_I_traj(1:s:end);
data.mu_delta_S_traj = mu_delta_S_traj(1:s:end);
data.mu_epsilon_I_traj = mu_epsilon_I_traj(1:s:end);
data.mu_epsilon_S_traj = mu_epsilon_S_traj(1:s:end);
data.R_traj_discrete = R_traj_discrete(1:s:end,:,:);
data.R_traj = R_traj_iv(1:s:end);
data.R_traj_v = R_traj_v(1:s:end);
data.Rt_SIR_iv = Rt_SIR_iv(1:s:end);
data.Rt_SIR_v = Rt_SIR_v(1:s:end);
data.S_traj_discrete = S_traj_discrete(1:s:end,:,:);
data.S_traj = S_traj_iv(1:s:end);
data.S_traj_v = S_traj_v(1:s:end);
data.varSe_traj = varSe_traj(1:s:end);
data.varId_traj = varId_traj(1:s:end);
data.dispSe_traj = dispSe_traj(1:s:end);
data.dispId_traj = dispId_traj(1:s:end);


save('dataIndnegBinomial.mat','params','data', 'params_v');

    
   
%     name1 = strcat('joint_S_and_I_distIV','.avi');
%     v1 = VideoWriter(name1);
%     v1.FrameRate = 2; 
% 
%     open(v1); 
%     for f = [1,76,101,151,251,351, 551,751,951, 1151,1201, 1251, 1301, 1351, 1401, 1451, 1501, 1651, 1851,2001, 2251, 2451, 2651, 2851, 3001, 3251, 3451, 3651]
%         f8 = figure(8);
%         set(gcf, 'Position', get(0, 'Screensize'));
%         dist = zeros(m,n);
%         dist(:,:) = S_traj_discrete(f,:,:);
%         sgtitle(strcat('Susceptible & Infectious Populations at time: ', num2str(t_span(f))));
% 
%         subplot(3,2,1)
%         imagesc(e,d,dist'); %maybe use pcolor function 
%         axis xy; 
%         xlabel('eps');
%         ylabel('delta'); 
%         zlabel('pop frac');
%         title('joint S');
% 
%         subplot(3,2,2)
%         dist = zeros(m,1);
%         dist(:,:) = marg_S_e_traj(f,:,:);
%         plot(e,dist); hold on;
%         xline(mu_epsilon_S_traj(f));
%         xlabel('eps');
%         ylabel('pop frac');
%         title('marg eps S')
% 
% 
%         dist = zeros(m,n);
%         dist(:,:) = I_traj_discrete(f,:,:);
% 
%         subplot(3,2,3)
%         imagesc(e,d,dist');
%         axis xy;
%         %surf(E,D,dist); 
%         xlabel('eps');
%         ylabel('delta');
%         zlabel('pop frac');
%         title('joint I');
% 
%         subplot(3,2,4)
%         dist = zeros(m,1);
%         dist(:,:) = marg_I_d_traj(f,:,:);
%         plot(d,dist); hold on;
%         xline(mu_delta_I_traj(f));
%         xlabel('delta');
%         ylabel('pop frac');
%         title('marg delta I')
% 
%         subplot(3,2,5)
%         q(1)=semilogy(params.t_span, S_traj,'--','Color',default_colors(1,:),'LineWidth',2); hold on;
%         q(2)=semilogy(params.t_span, I_traj,'--','Color',default_colors(2,:),'LineWidth',2); hold on;
%         q(3)=semilogy(params.t_span, R_traj,'--','Color',default_colors(3,:),'LineWidth',2); hold on;
%         xline(t_span(f));
%         axis([0 t_end 10^(-6) 10]);
%         xlabel('Time (days)'); ylabel({'Population'; 'Fraction'});
%         legend(q,{'S','I','R'},'Location', 'SouthEast');
%         
%         subplot(3,2,6)
%     
%         plot(params.t_span, 100.*FOSiv_traj,'--','Color',default_colors(2,:),'LineWidth',2);
%         axis([0 t_end 10^(-6) 10]);
%         ylim([0,105]);
%         xlabel('Time (days)'); ylabel({'Cumulated Infected (%)'});
%         xline(t_span(f));
%     
% 
%         frame = getframe(f8);
%         writeVideo(v1,frame);
%         close (f8);
%     end
% 
%     close(v1);
%     
    
end