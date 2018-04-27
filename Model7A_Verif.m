% This function returns the rate of change for each variable at
% equilibrium. All rates should thus be zero (if calculations are correct).

% VARIABLES
% H = Heterotrophic bacteria (H)
% C = Cyanobacteria (C)
% Z = Zooplankton (Z)
% Vh = Hetero phages (Vh)
% Vc = Cyanophages (Vc)
% Norg = Organic nitrogen (xOn)
% Nin = Inorganic nitrogen (xIn)


function xprime = Model7A_Verif(p,equi)

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
Vc = equi(6);
Vh = equi(7);

% Evaluate rates for variables at equilibrium
xprime1 = rH*H*Norg/(kOrg) - phiVH*H*Vh - gammaH*H*Z - sigmaH*H - lambdaH*H;
xprime2 = rC*C*Nin/(kIn) - phiVC*C*Vc - gammaC*C*Z - sigmaC*C - lambdaC*C;
xprime3 = pG*(bH*gammaH*H*Z/bZ + bC*gammaC*C*Z/bZ) - lambdaZ*Z - gammaZ*Z*Z;
xprime4 = betaH*phiVH*H*Vh - thetaH*Vh;
xprime5 = betaC*phiVC*C*Vc - thetaC*Vc;
xprime6 = -bH/epsiH*rH*H*Norg/(kOrg)+ bV*thetaH*Vh + bV*thetaC*Vc + (bH - bV*betaH)*phiVH*H*Vh...
    + (bC - bV*betaC)*phiVC*C*Vc + bH*sigmaH*H + bC*sigmaC*C + pEx*bH*gammaH*H*Z + pEx*bC*gammaC*C*Z;
xprime7 = -omega*(Nin - Nsub) + bH*(1-epsiH)*rH*H*Norg/(epsiH*(kOrg)) - bC*rC*C*Nin/(kIn) +bZ*lambdaZ*Z...
    + bH*lambdaH*H + bC*lambdaC*C + (1-pEx-pG)*bH*gammaH*H*Z + (1-pEx-pG)*bC*gammaC*C*Z;
xprime = [xprime1, xprime2, xprime3, xprime7, xprime6, xprime5,xprime4];
