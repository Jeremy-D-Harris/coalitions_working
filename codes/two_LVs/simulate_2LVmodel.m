function dy = simulate_2LVmodel(t,y,params)

% parameters to local variables
r_A = params.r_A;
gam_A = params.gam_A;
r_B = params.r_B;
gam_B = params.gam_B;

K = params.K;
d = params.d;
bet = params.bet;
phi = params.phi;
m = params.m;

% ODE Variables: LA, VA, LB, VB
LA = y(1);
VA = y(2);
LB = y(3);
VB = y(4);

N = LA + LB;

dy = zeros(4,1);

% ODE system: 
dy(1) = r_A*LA*(1-N/K) - phi*LA*VB  - gam_A*LA - d*LA;

dy(2) = bet*phi*LB*VA + bet*gam_A*LA - phi*N*VA - m*VA;

dy(3) = r_B*LB*(1-N/K) - phi*LB*VA - gam_B*LB  - d*LB;

dy(4) = bet*phi*LA*VB + bet*gam_B*LB - phi*N*VB - m*VB;

if y(3) < 1
    
    y(3) == 0;
    dy(3) = 0;
    
end


end
