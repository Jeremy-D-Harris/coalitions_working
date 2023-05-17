function Rt_calc = get_Rt_SIR(params,y)

% calculate effective reproduction number, R_t

% parameters to local variables
bet = params.beta; 
gam = params.gamma; 
mean = params.mean;

for count=1:length(params.t_span)

% t = params.t_span(count);
S = y(count,1);

% transmissions
T = mean*bet*S;

% transitions
Sigma = -gam;

NGM = -T*(1/Sigma); % next generation matrix: -T*Sigma^-1

eigen_values = eig(NGM);

Rt_calc(count) = max(eigen_values);

end
