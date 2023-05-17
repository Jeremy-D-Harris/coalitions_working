function dydt = simulate_SEILV_gamma_q(t,y,params)

% parameters to local variables
% parameters
% conversion_efficiency = params.conversion_efficiency;
% d_R = params.d_R;
mu_max = params.mu_max;
K = params.K;
% R_in = params.R_in;
d_S = params.d_S;
d_E = params.d_E;
d_L = params.d_L;
d_I = params.d_I;
lam = params.lambda;
eta = params.eta;
bet = params.bet;
phi = params.phi;
m_V = params.m_V;
% gam = params.gamma;
% q = params.q;

% alpha_s = params.alpha_s;
% J = params.J;
% flask_volume = params.flask_volume;

gamma_mesh = params.gamma_mesh;
q_mesh = params.q_mesh;

n = length(gamma_mesh); % gamma goes across columns
m = length(q_mesh); % q goes down rows

% S, E, I, L, V
dist_S = y(1:(m*n)); dist_E = y((m*n+1):(2*m*n)); 
dist_I = y((2*m*n+1):(3*m*n)); dist_L = y((3*m*n+1):(4*m*n));
V = y(4*m*n+1);

%reshaping
dist_S = reshape(dist_S,m,n);
dist_E = reshape(dist_E,m,n);
dist_I = reshape(dist_I,m,n);
dist_L = reshape(dist_L,m,n);
N = sum(sum(dist_S+dist_E+dist_I+dist_L));

%differential equations
dSdt = -phi*V*dist_S + mu_max*dist_S*(1-N/K) - d_S*dist_S;

dEdt = phi*V*dist_S - lam*dist_E - d_E*dist_E;

dIdt = lam*(1-q_mesh).*dist_E + gamma_mesh.*dist_L - eta*dist_I - - d_I*dist_I; 

dLdt = lam*(q_mesh.*dist_E) - gamma_mesh.*dist_L + mu_max*dist_L*(1-N/K) - d_L*dist_L; 

dVdt = bet*eta*sum(sum(dist_I))-phi*V*sum(sum(dist_S)) - m_V*V;


%reshaping
dSdt = reshape(dSdt,[m*n,1]);
dEdt = reshape(dEdt,[m*n,1]);
dIdt = reshape(dIdt,[m*n,1]);
dLdt = reshape(dLdt,[m*n,1]);

dydt = [dSdt; dEdt; dIdt; dLdt; dVdt];

        
end




