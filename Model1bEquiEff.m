%This function returns the equilibrium values for Model 1B (with growth efficiency).

function equiValues = Model1bEquiEff(p)

% Assign value to parameters
    rHmax = p(1);
    kInH = p(2);
    kOrgH = p(3);
    lambdaH = p(4);
    sigmaH = p(5);
    qH = p(6);
    psiH = p(7);
    N = p(8);
    epsiH = p(9);
    
% Evaluate equilibrium expressions for parameters
    % Proportions    
    propH = (kInH*lambdaH*rHmax + lambdaH*N*rHmax - N*rHmax*rHmax - epsiH*kInH*lambdaH*sigmaH - epsiH*kOrgH*lambdaH*sigmaH - epsiH*lambdaH*N*sigmaH +kInH*rHmax*sigmaH - epsiH*kInH*rHmax*sigmaH + epsiH*kOrgH*rHmax*sigmaH ...
        + N*rHmax*sigmaH - epsiH*kInH*sigmaH^2 + epsiH^2*kInH*sigmaH^2-epsiH*kOrgH*sigmaH^2+epsiH^2*kOrgH*sigmaH^2 - epsiH*N*sigmaH^2 +epsiH^2*N*sigmaH^2)/(N*(lambdaH - rHmax+sigmaH-epsiH*sigmaH)*(rHmax - epsiH*sigmaH));
    propPin = (-kInH*lambdaH - kInH*sigmaH + epsiH*kInH*sigmaH)/(N*(lambdaH - rHmax + sigmaH - epsiH*sigmaH));
    
    % Absolute units
    H = propH*N/qH;
    Pin = propPin*N;
    Porg = (1 - propH - propPin)*N;
    
    equiValues = [H, Pin,Porg];