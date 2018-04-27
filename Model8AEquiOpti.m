% This function returns the value to be minimized during the optimization
% procedure. 

function y = Model8AEquiOpti(p)

equiVfunction = Model8AEqui(p);

% Target values
target = [6*10^8, 1*10^8, 4*10^4, 7*10^-3, 0.1, 9*10^9, 1.5*10^9];
% Weighted difference between equilibrium and target
y = sum(((equiVfunction-target)./target).^2);

end