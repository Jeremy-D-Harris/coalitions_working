function dy = simulate_LVmodel(t,y,params)

% parameters to local variables
r = params.r;
K = params.K;
d = params.d;
gam = params.gam;
bet = params.bet;
phi = params.phi;
m = params.m;

% L, V
L = y(1);
V = y(2);

% Variables
% ODE
dy = zeros(2,1);

% ODE system
dy(1) = r*L*(1-L/K) - gam*L - d*L;
dy(2) = bet*gam*L -phi*L*V - m*V;

end
