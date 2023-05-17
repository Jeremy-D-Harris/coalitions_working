function eigen_vector = get_eigendirection_SIRdelta(params)

% parameters to local variables
bet=params.beta; 
gam = params.gamma;
delta_S = params.delta_S;

% disease free state
S = 1; I = 0; R = 0; delta_I=1;

dSdt = [-bet*delta_I*I,-bet*delta_I*S,0,-bet*I*S];

dIdt = [bet*delta_I*I,bet*delta_I*S - gam,0,bet*I*S];

dRdt = [0,gam,0,0];

ddelta_Idt = [bet*delta_I*(delta_S-delta_I),0,0,bet*S*delta_S-2*bet*S*delta_I];

A = [dSdt;dIdt;dRdt;ddelta_Idt];

[eigen_directions, eigen_values] = eig(A); % get eigenvalues/eigenvectors
eigen_values
[val, ind] = max(diag(eigen_values));
eigen_vector = eigen_directions(:,ind); % corresponds to max eigenvalue
