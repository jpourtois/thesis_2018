%This function calculates the equilibrium values for Model 3B.

function model3B = Model3BEqui(p)

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
    
    % Evaluate equilibrium expressions for parameters
    eCons = qZ*gammaZ/omega;
    
    c3 = -gammaH*eCons;
    c2 = eCons*(rH - lambdaH - sigmaH);
    c1 = gammaH*(pSub + kInH);
    c0 = lambdaH*(pSub + kInH) + sigmaH*(pSub + kInH) - rH*pSub;
    
    pRoots = [c3 c2 c1 c0];
    equiZMult = roots(pRoots);
    equiZMult = equiZMult(imag(equiZMult)==0);
    
    for i = 1:length(equiZMult)
        
        equiZ = equiZMult(i);
        equiPin = pSub - eCons*equiZ^2;
        equiPon = qH*sigmaH/psiH + pEx*gammaH*qH*equiZ/psiH;
        equiH = qZ*(lambdaZ + gammaZ*equiZ)/(pG*qH*gammaH);
        
        model3B(i,:) = [equiH, equiZ, equiPin, equiPon];
    end 
        