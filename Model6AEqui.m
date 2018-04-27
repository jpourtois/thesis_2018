%This function calculates the equilibrium values for Model 6A. 

function model6A = Model6AEqui(p)

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
    equiC = thetaC/(betaC*phiC);
    equiZ = pG*(qH*gammaH*equiH + qC*gammaC*equiC)/(qZ*gammaZ) - lambdaZ/gammaZ;
    equiPin = pSub - qZ*gammaZ*equiZ*equiZ/omega;
    equiVh = rH*equiPin/((equiPin + kInH)*phiH) - lambdaH/phiH - sigmaH/phiH - gammaH*equiZ/phiH;
    equiVc = rC*equiPin/((equiPin + kInC)*phiC) - lambdaC/phiC - sigmaC/phiC - gammaC*equiZ/phiC;
    equiPorg = (qH*sigmaH*equiH + qC*sigmaC*equiC + pEx*equiZ*(gammaH*qH*equiH + gammaC*qC*equiC) + thetaH*equiVh*qV + thetaC*equiVc*qV + (qH - betaH*qV)*phiH*equiVh*equiH + (qC - betaC*qV)*phiC*equiVc*equiC)/(psiH*equiH + psiC*equiC);
    

    model6A = [equiH,equiC,equiZ,equiPin,equiPorg,equiVh,equiVc];

end
