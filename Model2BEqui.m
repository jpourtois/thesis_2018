%This function calculates the equilibrium values for Model 2B. 


function model2B = Model2BEqui(p)
    
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
    equiZ = (lambdaH + sigmaH - sigmaH*epsiH)/(epsiH*pEx*gammaH - gammaH);
    equiNin = Nsub - bZ*gammaZ*equiZ^2/omega;
    equiNorg = kOrg*(gammaH*equiZ + lambdaH + sigmaH)/rH;
    equiH = bZ*(gammaZ*equiZ + lambdaZ)/(pG*bH*gammaH);
      
    model2B = [equiH, equiZ, equiNin, equiNorg];
end
   
  