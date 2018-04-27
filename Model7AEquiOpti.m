% This function returns the value to be minimized during the optimization
% procedure. 

function y = Model7AEquiOpti(p)

equiV = Model7AEqui(p);

% Target values
target = [6*10^8, 1*10^8, 4*10^4, 0.1, 5, 1.5*10^9, 9*10^9];
% Weighted difference between equilibrium and target
y = sum(((equiV-target)./target).^2);

end