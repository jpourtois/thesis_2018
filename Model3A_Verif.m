% This function returns the rate of change for each variable at
% equilibrium. All rates should thus be zero (if calculations are correct). 

% VARIABLES
% H = Heterotrophic bacteria 
% Z = Zooplankton 
% Pin = Inorganic phosphorus 
% Porg = Organic phosphorus 

function xprime = Model3A_Verif(p,equi)

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
    Z = equi(2);
    Pin = equi(3);
    Porg = equi(4); 
    Vh = equi(5);
    
    % Evaluate rates for variables at equilibrium
    xprime1 = rH*H*Pin/(Pin + kInH) - phiH*H*Vh - gammaH*H*Z - lambdaH*H - sigmaH*H;
    xprime2 = pG*Z*qH*gammaH*H/qZ - lambdaZ*Z - gammaZ*Z^2;
    xprime3 = - omega*(Pin - pSub) - qH*rH*H*Pin/(Pin + kInH) + qZ*lambdaZ*Z + qH*lambdaH*H + (1- pG - pEx)*qH*gammaH*H*Z + psiH*H*Porg;
    xprime4 = qV*thetaH*Vh + (qH - qV*betaH)*phiH*H*Vh + qH*sigmaH*H + pEx*qH*gammaH*H*Z - psiH*H*Porg;
    xprime5 = betaH*phiH*H*Vh - thetaH*Vh;
    
    xprime = [xprime1, xprime2, xprime3, xprime4, xprime5];