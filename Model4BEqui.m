%This function calculates the equilibrium values for Model 4B. 

function model4B = Model4BEqui(p)
    
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
    
    mC = lambdaC + sigmaC; 
 
    %Solve quadratic expression to find Z
    c2 = -rC*bZ*gammaZ/(kIn*omega);
    c1 = -gammaC;
    c0 = -sigmaC - lambdaC + rC*Nsub/kIn;
    pRoots = [c2 c1 c0];
    equiZMult = roots(pRoots);

    %Keep only the real roots
    equiZMult = equiZMult(imag(equiZMult)==0); 
    
    model4B = zeros(length(equiZMult),5);
    
    %For each value of Z, evaluate the other variables
    for n = 1:length(equiZMult)
                
        % Evaluate equilibria
        NVequiZ = equiZMult(n);
        NVequiXin = Nsub - bZ*gammaZ*NVequiZ*NVequiZ/omega;
        NVequiXon = kOrg*(gammaH*NVequiZ + sigmaH + lambdaH)/rH;
        numNVequiH = pEx*bZ*NVequiZ*(lambdaZ + gammaZ*NVequiZ)/pG + sigmaC*bZ*(lambdaZ + gammaZ*NVequiZ)/(pG*gammaC);
        denomNVequiH = bH*rH*NVequiXon/(epsiH*kOrg) - bH*sigmaH + sigmaC*bH*gammaH/gammaC;
        NVequiH = numNVequiH/denomNVequiH;
        NVequiC = (gammaZ*NVequiZ + lambdaZ - pG*bH*gammaH*NVequiH/bZ)*bZ/(pG*bC*gammaC);
        
        model4B(n,:) = [NVequiH, NVequiC, NVequiZ, NVequiXin, NVequiXon];
    end
    