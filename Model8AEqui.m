% This function calculates the equilibrium values for Model 7A.

function model8A = Model8AEqui(p)

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

    % Evaluate equilibria
    equiH = thetaH/(betaH*phiH);
    equiC = thetaC/(betaC*phiC);
    equiZ = pG*(qH*gammaH*equiH + qC*gammaC*equiC)/(qZ*gammaZ) - lambdaZ/gammaZ;
    equiPin = pSub - qZ*gammaZ*equiZ*equiZ/omega;
    sH = rH*equiPin/kInH - lambdaH - sigmaH - gammaH*equiZ;
    sC = rC*equiPin/kInC - lambdaC - sigmaC - gammaC*equiZ;
    sPon = qH*sigmaH*equiH + qC*sigmaC*equiC + pEx*equiZ*(gammaH*qH*equiH + gammaC*qC*equiC);
    num = sPon + sH/phiH*(thetaH*qV + (qH - betaH*qV)*phiH*equiH) + sC/phiC*(thetaC*qV + (qC - betaC*qV)*phiC*equiC);
    denom = qH*rH*equiH/(epsiH*kOrgH) + qC*rC*equiC/(epsiC*kOrgC) - rH*(thetaH*qV + (qH - betaH*qV)*phiH*equiH)/(kOrgH*phiH)...
        - rC*(thetaC*qV + (qC - betaC*qV)*phiC*equiC)/(kOrgC*phiC);
    equiPon = num/denom;
    equiVh = rH*equiPin/(kInH*phiH) + rH*equiPon/(kOrgH*phiH) - lambdaH/phiH - sigmaH/phiH - gammaH*equiZ/phiH;
    equiVc = rC*equiPin/(kInC*phiC) + rC*equiPon/(kOrgC*phiC) - lambdaC/phiC - sigmaC/phiC - gammaC*equiZ/phiC;


    model8A = [equiH, equiC, equiZ, equiPin, equiPon, equiVh, equiVc];

end
