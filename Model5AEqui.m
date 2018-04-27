%This function calculates the equilibrium values for Model 5A. 

function model5A = Model5AEqui(p)

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
    model5A = zeros(1,6);
    
    equiH = thetaH/(betaH*phiH);
    
    c3 = rC*qZ^2*gammaZ^2/(epsiC*kInC*gammaC*pG*omega);
    c2 = qZ*gammaZ + qH*rH*equiH*qZ*gammaZ/(kInH*omega) + rC*lambdaZ*qZ^2*gammaZ/(kInC*gammaC*pG*omega)...
        - rC*qH*gammaH*equiH*qZ*gammaZ/(kInC*gammaC*omega) + qH*(1-epsiH)*rH*kOrgC*qZ*gammaZ*equiH/(epsiH*kOrgH*kInC*omega) ...
        + (1-epsiC)*rC*qZ*gammaZ/(epsiC*kInC*omega*gammaC)*(lambdaZ*qZ/pG - qH*gammaH*equiH) + (1 - epsiC)*gammaZ*qZ/(epsiC*pG) + (1-pG-pEx)*gammaZ*qZ/pG;
    c1 = -rC*gammaZ*qZ*pSub/(kInC*gammaC*pG)+ qH*(1-epsiH)*rH*equiH*kOrgC*gammaC/(epsiH*kOrgH*rC) + (1 - epsiC)*gammaZ*qZ*(lambdaC+sigmaC - rC*pSub/kInC)/(epsiC*gammaC*pG) ...
        + (1 - epsiC)*gammaC/epsiC*(lambdaZ*qZ/(gammaC*pG) - qH*gammaH*equiH/gammaC) + qZ*lambdaZ + gammaZ*lambdaC*qZ/(gammaC*pG) + (1-pG-pEx)*(qH*gammaH*equiH + lambdaZ*qZ/pG - qH*gammaH*equiH);    
    c0 = - qH*rH*equiH*pSub/kInH - rC*lambdaZ*qZ*pSub/(kInC*gammaC*pG) + rC*qH*gammaH*equiH*pSub/(kInC*gammaC) ...
        + qH*(1 - epsiH)*rH*equiH*kOrgC/(epsiH*kOrgH*rC)*(lambdaC + sigmaC - rC*pSub/kInC) + (1-epsiC)/epsiC*(lambdaZ*qZ/gammaC*pG - qH*gammaH*equiH/gammaC)*(lambdaC+sigmaC-rC*pSub/kInC) ...
        + qH*lambdaH*equiH + lambdaC*lambdaZ*qZ/(gammaC*pG) - lambdaC*qH*gammaH*equiH/gammaC;
    
    pRoots = [c3 c2 c1 c0];
    equiZMult = roots(pRoots);
    equiZMult = equiZMult(imag(equiZMult)==0);
    
    for i = 1:length(equiZMult)
        
        equiZ = equiZMult(i);
        equiPin = pSub - qZ*gammaZ*equiZ*equiZ/omega;
        equiC = (lambdaZ*qZ + gammaZ*equiZ*qZ)/(qC*gammaC*pG) - qH*gammaH*equiH/(qC*gammaC);
        equiPorg = kOrgC*(gammaC*equiZ + lambdaC + sigmaC - rC*equiPin/kInC)/rC;
        equiVh = rH*equiPin/(kInH*phiH) + rH*equiPorg/(kOrgH*phiH) - gammaH*equiZ/phiH - lambdaH/phiH - sigmaH/phiH;
        
        equi = [equiH,equiC,equiZ,equiPin,equiPorg,equiVh];
       
        model5A(i,:) = equi;
    end