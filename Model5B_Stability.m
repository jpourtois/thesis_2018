% This function returns the eigen values of the Jacobian matrix for an
% equilibrium

% VARIABLES
% H = Heterotrophic bacteria 
% C = Cyanobacteria 
% Z = Zooplankton 
% Pin = Inorganic phosphorus 
% Porg = Organic phosphorus 

function eigenValues = Model5B_Stability(p,equi)

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
    
    
    syms H C Z Pin Porg 
    
    % Differential equations
    xprime1 = rH*H*Pin/kInH + rH*H*Porg/kOrgH - gammaH*H*Z - lambdaH*H - sigmaH*H;
    xprime2 = rC*C*Pin/kInC + rC*C*Porg/kOrgC - gammaC*C*Z - lambdaC*C - sigmaC*C;
    xprime3 = pG*Z*(qH*gammaH*H + qC*gammaC*C)/qZ - lambdaZ*Z - gammaZ*Z^2;
    xprime4 = -omega*(Pin - pSub) - qH*rH*H*Pin/kInH - qC*rC*C*Pin/kInC + qH*(1-epsiH)*rH*H*Porg/(epsiH*kOrgH) + qC*(1-epsiC)*rC*C*Porg/(epsiC*kOrgC)...
        + qZ*lambdaZ*Z + qH*lambdaH*H + qC*lambdaC*C + (1 - pG - pEx)*(qH*gammaH*H*Z + qC*gammaC*C*Z);
    xprime5 = - qH*rH*H*Porg/(epsiH*kOrgH) - qC*rC*C*Porg/(epsiC*kOrgC) + qH*sigmaH*H + qC*sigmaC*C + pEx*Z*(qH*gammaH*H + qC*gammaC*C);
    
    xprime = [xprime1, xprime2, xprime3, xprime4, xprime5];
    
    variable = [H, C, Z, Pin, Porg];
    
    % Make jacobian matrix
    jacMatrix = jacobian(xprime, variable);

    % Evaluate at equilibrium
    jacMatrixEval = subs(jacMatrix, variable, equi);

    % Find eigen values
    eigenValues = eig(jacMatrixEval);