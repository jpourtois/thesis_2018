% This function returns the rate of change for each variable at
% equilibrium. All rates should thus be zero (if calculations are correct).

% VARIABLES
% H = Heterotrophic bacteria (H)
% C = Cyanobacteria (C)
% Z = Zooplankton (Z)
% Vh = Hetero phages (Vh)
% Vc = Cyanophages (Vc)
% Porg = Organic nitrogen (Porg)
% Pin = Inorganic nitrogen (PIn)

function xprime = Model8A_Verif(p,equi)

    % Assign value to parameters
    rH = p(1);
    kInH = p(2);
    lambdaH = p(3);
    sigmaH = p(4);
    gammaH = p(5);
    phiH = p(6);
    thetaH = p(7);
    pG = p(8);
    qH = p(9);
    qZ = p(10);
    lambdaZ = p(11);
    gammaZ = p(12);
    omega = p(13);
    pSub = p(14);
    pEx = p(15);
    psiH = p(16);
    qV = p(17);
    betaH = p(18);
    rC = p(19);
    kInC = p(20);
    lambdaC = p(21);
    sigmaC = p(22);
    gammaC = p(23);
    phiC = p(24);
    betaC = p(25);
    thetaC = p(26);
    psiC = p(27);
    qC = p(28);
    kOrgH = p(29);
    kOrgC = p(30);
    epsiH = p(31);
    epsiC = p(32);
    
    % Assign equilibrium values to variables
    H = equi(1);
    C = equi(2);
    Z = equi(3);
    Pin = equi(4);
    Porg = equi(5);
    Vh = equi(6);
    Vc = equi(7);
    
    % Evaluate rates for variables at equilibrium
    xprime1 = rH*H*Pin/kInH + rH*H*Porg/kOrgH - phiH*H*Vh - gammaH*H*Z - lambdaH*H - sigmaH*H;
    xprime2 = rC*C*Pin/kInC + rC*C*Porg/kOrgC - phiC*C*Vc - gammaC*C*Z - lambdaC*C - sigmaC*C;
    xprime3 = pG*Z*(qH*gammaH*H + qC*gammaC*C)/qZ - lambdaZ*Z - gammaZ*Z^2;
    xprime4 = -omega*(Pin - pSub) - qH*rH*H*Pin/kInH - qC*rC*C*Pin/kInC + qH*(1-epsiH)*rH*H*Porg/(epsiH*kOrgH) + qC*(1-epsiC)*rC*C*Porg/(epsiC*kOrgC)...
        + qZ*lambdaZ*Z + qH*lambdaH*H + qC*lambdaC*C + (1 - pG - pEx)*(qH*gammaH*H*Z + qC*gammaC*C*Z);
    xprime5 = - qH*rH*H*Porg/(epsiH*kOrgH) - qC*rC*C*Porg/(epsiC*kOrgC) + qV*thetaH*Vh + qV*thetaC*Vc + (qH - qV*betaH)*phiH*H*Vh ...
        + (qC - qV*betaC)*phiC*C*Vc + qH*sigmaH*H + qC*sigmaC*C + pEx*Z*(qH*gammaH*H + qC*gammaC*C);
    xprime6 = betaH*phiH*H*Vh - thetaH*Vh;
    xprime7 = betaC*phiC*C*Vc - thetaC*Vc;
    
    xprime = [xprime1, xprime2, xprime3, xprime4, xprime5, xprime6, xprime7];
    
    