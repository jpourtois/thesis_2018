% This function returns the rate of change for each variable at
% equilibrium. All rates should thus be zero (if calculations are correct). 

% VARIABLES
% H = Heterotrophic bacteria (H)
% C = Cyanobacteria (C)
% Z = Zooplankton (Z)
% Nin = Inorganic nitrogen (NIn)
% Norg = Organic nitrogen (NOn)

function xprime = Model4B_Verif(p,equi)

% Assign value to parameters
rH = p(1);
rC = p(2);
kOrg = p(3);
kIn = p(4);
gammaH = p(5);
gammaC = p(6);
sigmaH = p(7);
sigmaC = p(8);
lambdaH = p(9);
lambdaC = p(10);
pG = p(11);
bH = p(12);
bZ = p(13); 
lambdaZ = p(14);
gammaZ = p(15);
pEx = p(16);
omega = p(17);
Nsub = p(18);
phiVH = p(19);
phiVC = p(20);
bV = p(21);
thetaH = p(22);
thetaC = p(23);
betaH = p(24);
betaC = p(25);
bC = p(26);
epsiH = p(27);

% Assign equilibrium values to variables
H = equi(1);
C = equi(2);
Z = equi(3);
Nin = equi(4);
Norg = equi(5);

% Evaluate rates for variables at equilibrium
xprime1 = rH*H*Norg/kOrg - gammaH*H*Z - lambdaH*H - sigmaH*H;
xprime2 = rC*C*Nin/kIn - gammaC*C*Z - lambdaC*C - sigmaC*C;
xprime3 = pG*Z*(bH*gammaH*H + bC*gammaC*C)/bZ - lambdaZ*Z - gammaZ*Z^2;
xprime4 = - omega*(Nin - Nsub) - bC*rC*C*Nin/kIn + bH*(1 - epsiH)*rH*H*Norg/(epsiH*kOrg) ...
    + bZ*lambdaZ*Z + bH*lambdaH*H + bC*lambdaC*C + Z*(1 - pG - pEx)*(bH*gammaH*H + bC*gammaC*C);
xprime5 = -bH*rH*H*Norg/(kOrg*epsiH) + bH*sigmaH*H + bC*sigmaC*C + pEx*Z*(bH*gammaH*H + bC*gammaC*C);

xprime = [xprime1, xprime2, xprime3, xprime4, xprime5];