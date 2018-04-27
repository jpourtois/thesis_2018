%This function returns the equilibrium values for Model 1A.

function equiValues = Model1AEqui(p)

% Assign value to parameters
    rHmax = p(1);
    kInH = p(2);
    kOrgH = p(3);
    lambdaH = p(4);
    sigmaH = p(5);
    qH = p(6);
    psiH = p(7);
    N = p(8);
    epsiH= p(9);
 
% Evaluate equilibrium expressions for parameters
    % Proportions
    propH = (kInH*lambdaH*psiH + lambdaH*N*psiH - N*psiH*rHmax + kInH*psiH*sigmaH + N*psiH*sigmaH - lambdaH*qH*sigmaH + qH*rHmax*sigmaH - qH*sigmaH*sigmaH)/(N*psiH*(lambdaH - rHmax + sigmaH));
    propPin = (-kInH*lambdaH - kInH*sigmaH)/(N*(lambdaH - rHmax + sigmaH));
    
    % Absolute units
    H = propH*N/qH;
    Pin = propPin*N;
    Porg = (1 - propH - propPin)*N;
    
    equiValues = [H, Pin,Porg];
    