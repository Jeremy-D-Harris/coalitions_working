function eigen_vector = get_eigendirection_SIR(params)

%% parameters
% parameters to local variables
bet=params.beta; 
gam = params.gamma;

%% SIR Version
S_init = 1; I_init = 0; R_init = 0; 

dSdt = [-bet*epsilon_S*delta_I*I_init,-bet*epsilon_S*delta_I*S_init,0];

dIdt = [bet*epsilon_S*delta_I*I_init,bet*epsilon_S*delta_I*S_init - gam,0];

dRdt = [0,gam,0];


A = [dSdt;dIdt;dRdt];

[eigen_directions, eigen_values] = eig(A); % get eigenvalues/eigenvectors
% eigen_directions
% eigen_values
[val, ind] = max(diag(eigen_values));
eigen_vector = eigen_directions(:,ind);% corresponds to max eigenvalue
end
