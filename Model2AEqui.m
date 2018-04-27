%This function returns the equilibrium values for Model 2A. 


function model2A = Model2AEqui(p)
    
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
    equiZ = pG*bH*gammaH*equiH/(bZ*gammaZ) - lambdaZ/gammaZ;
    numEquiVh = -(bH*equiH*(sigmaH + pEx*gammaH*equiZ - gammaH*equiZ/epsiH - lambdaH/epsiH - sigmaH/epsiH));
    denomEquiVh = bV*thetaH + (bH - bV*betaH)*phiVH*equiH - bH*phiVH*equiH/epsiH;
    equiVh = numEquiVh/denomEquiVh;
    equiNorg = kOrg*epsiH*(bV*thetaH*equiVh + (bH - bV*betaH)*phiVH*equiH*equiVh + bH*sigmaH*equiH + pEx*bH*gammaH*equiH*equiZ)/(bH*rH*equiH);
    equiNin = Nsub - bZ*gammaZ*equiZ^2/omega;
    
    model2A = [equiH, equiZ, equiNin, equiNorg, equiVh];
end
   
  