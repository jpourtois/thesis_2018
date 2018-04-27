%This function calculates the equilibrium values for Model 4A. 

function model4A = Model4AEqui(p)
    
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
    
    % Evaluate equilibrium expressions for parameters
    equiH = thetaH/(betaH*phiVH);
    
    %Solve quadratic expression to find Z
    c2 = -rC*bZ*gammaZ/(kIn*omega);
    c1 = -gammaC;
    c0 = -lambdaC - sigmaC + rC*Nsub/kIn;
    
    pRoots = [c2 c1 c0];
    equiZMult = roots(pRoots);
    equiZMult = equiZMult(imag(equiZMult)==0);
    
    model4A = zeros(length(equiZMult),6);
    
    for i = 1:length(equiZMult)
        
        equiZ = equiZMult(i);
        equiNin = Nsub - bZ*gammaZ*equiZ*equiZ/omega;
        equiC = (bZ*gammaZ*equiZ + bZ*lambdaZ - pG*bH*gammaH*equiH)/(pG*bC*gammaC);
        numVh = bH*sigmaH*equiH + bC*sigmaC*equiC + pEx*equiZ*(bH*gammaH*equiH + bC*gammaC*equiC) - bH*(gammaH*equiH*equiZ + lambdaH*equiH + sigmaH*equiH)/epsiH;
        denomVh = bH*phiVH*equiH/epsiH - bV*thetaH - (bH - bV*betaH)*phiVH*equiH;
        equiVh = numVh/denomVh;
        equiNorg = kOrg*(phiVH*equiVh + gammaH*equiZ + lambdaH + sigmaH)/rH;
        
        model4A(i,:) = [equiH, equiC, equiZ, equiNin, equiNorg, equiVh];
    end
    