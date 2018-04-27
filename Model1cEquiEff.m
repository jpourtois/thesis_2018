%This function returns the equilibrium values for Model 1C (with growth efficiency).

function equiValues = Model1cEquiEff(p)

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
    propH = (kInH*lambdaH - N*rHmax + kInH*sigmaH - epsiH*kInH*sigmaH + epsiH*kOrgH*sigmaH)/(- N*rHmax);
    propPin = (kInH*lambdaH + kInH*sigmaH - epsiH*kInH*sigmaH)/(N*rHmax);
    
    % Absolute units
    H = propH*N/qH;
    Pin = propPin*N;
    Porg = (1 - propH - propPin)*N;
    
    equiValues = [H, Pin,Porg];