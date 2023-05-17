function dy = ODE_onelysogen_onevirus_vf(y,params)

L = y(1);
V = y(2);

% parameters
K= params.K;
r = params.r;
phi = params.phi;
bet = params.bet;
d = params.d;
m = params.m;
eta = params.eta;

% Variables
% ODE
dy = zeros(2,1);


% ODE system
dy(1) = r*L*(1-L/K) - eta*L - d*L;
dy(2) = bet*eta*L -phi*L*V - m*V;

