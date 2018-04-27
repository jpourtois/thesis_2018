%This function calculates the equilibrium values for Model 3A.

function model3A = Model3AEqui(p)

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
    equiH = thetaH/(betaH*phiH);
    equiZ = pG*qH*gammaH*equiH/(qZ*gammaZ) - lambdaZ/gammaZ;
    equiPin = pSub - qZ*gammaZ*equiZ^2/omega;
    equiVh = rH*equiPin/((equiPin + kInH)*phiH) - lambdaH/phiH - sigmaH/phiH - gammaH*equiZ/phiH;
    equiPorg = (qV*thetaH*equiVh + (qH - qV*betaH)*phiH*equiH*equiVh + qH*sigmaH*equiH + pEx*qH*gammaH*equiH*equiZ)/(psiH*equiH);
    
    model3A = [equiH, equiZ, equiPin, equiPorg, equiVh];
    
    
  