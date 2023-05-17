function dydt = simulateClassic(params)

%% parameters
% parameters to local variables
bet=params.beta;
gam = params.gamma;

%% SIR Version
S = y(1); I = y(2); R = y(3);

dSdt = -bet*I*S;

dIdt = bet*I*S - gam*I;

dRdt = gam*I;


A = [dSdt;dIdt;dRdt];

[eigen_directions, eigen_values] = eig(A); % get eigenvalues/eigenvectors
% eigen_directions
% eigen_values
[val, ind] = max(diag(eigen_values));
eigen_vector = eigen_directions(:,ind);% corresponds to max eigenvalue
end
