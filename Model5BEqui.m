%This function calculates the equilibrium values for Model 5B. 

function model5B = Model5BEqui(p)

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
    qVh = p(17);
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
    model5B = zeros(1,5);
    
    aleph = rH*kOrgC/(rC*kOrgH);
    
    c2 = -qZ*gammaZ/omega*(rH/kInH - rH*kOrgC/(kInC*kOrgH));
    c1 = aleph*gammaC - gammaH;
    c0 = pSub*(rH/kInH - rH*kOrgC/(kInC*kOrgH)) + aleph*lambdaC + aleph*sigmaC - lambdaH - sigmaH;
    
    pRoots = [c2 c1 c0];
    equiZMult = roots(pRoots);
    equiZMult = equiZMult(imag(equiZMult)==0);
    
    for i = 1:length(equiZMult)
        
        NVequiZ = equiZMult(i);
        NVequiPin = pSub - qZ*gammaZ*NVequiZ*NVequiZ/omega;
        NVequiPon = kOrgC/rC*(-rC*NVequiPin/kInC + lambdaC + sigmaC + gammaC*NVequiZ);
        numC = (-rH*NVequiPon*qZ/(epsiH*kOrgH*pG*gammaH) + sigmaH*qZ/(pG*gammaH) + pEx*NVequiZ*qZ/pG)*(lambdaZ + gammaZ*NVequiZ);
        denomC = rH*NVequiPon*qC*gammaC/(epsiH*kOrgH*gammaH) - qC*rC*NVequiPon/(epsiC*kOrgC) - sigmaH*qC*gammaC/gammaH + qC*sigmaC;
        NVequiC = -numC/denomC;
        NVequiH = qZ/(pG*qH*gammaH)*(-pG*qC*gammaC*NVequiC/qZ + lambdaZ + gammaZ*NVequiZ);
        
        NVequi = [NVequiH,NVequiC,NVequiZ,NVequiPin, NVequiPon];
       
        model5B(i,:) = NVequi;
    end