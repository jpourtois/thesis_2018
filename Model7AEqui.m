% This function calculates the equilibrium values for Model 7A. 

function model7A = Model7AEqui(p)

    % Assign value to parameters
    rH = p(1);
    rC = p(2);
    kOrg = p(3);
    kIn = p(4);
    gammaH = p(5);
    gammaC = p(6);
    sigmaH= p(7);
    sigmaC= p(8);
    lambdaH= p(9);
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

    % Evaluate equilibria
    equiH = thetaH/(betaH*phiVH);
    equiC = thetaC/(betaC*phiVC);
    equiZ = pG*(bH*gammaH*equiH + bC*gammaC*equiC)/(bZ*gammaZ) - lambdaZ/gammaZ;
    equiNin = Nsub - bZ*gammaZ*equiZ*equiZ/omega;
    equiVc = (rC*equiNin/kIn - gammaC*equiZ - lambdaC - sigmaC)/phiVC;
    numEquiVh = -bH*equiH*(gammaH*equiZ + lambdaH+ sigmaH)/epsiH + bV*thetaC*equiVc + (bC - bV*betaC)*phiVC*equiC*equiVc + pEx*bH*gammaH*equiH*equiZ + pEx*bC*gammaC*equiC*equiZ + bH*sigmaH*equiH + bC*sigmaC*equiC;
    equiVh = numEquiVh/(bH*phiVH*equiH/epsiH - bV*thetaH - (bH - bV*betaH)*phiVH*equiH);
    equiNorg = kOrg*(gammaH*equiZ + phiVH*equiVh + sigmaH + lambdaH)/rH;
    %vtob=(equiVh+equiVc)/(equiH+equiC);
    %fracvh=phiVH*equiVh/(phiVH*equiVh+gammaH*equiZ);
    %fracvc=phiVC*equiVc/(phiVC*equiVc+gammaC*equiZ);
    
    model7A = [equiH, equiC, equiZ, equiNin, equiNorg, equiVc, equiVh];


end 