% Generate random sets of parameters to evaluate equilibrium of Model 2A and
% 2B

% Define parameters' minimum
pMin = (1:26);

pMin(1) = 0.5; %rH
pMin(2) = 0.1; %rC
pMin(3) = 0.25; %kOn
pMin(4) = 0.05; %kIn
pMin(5) = 10^(-6); %gammaH
pMin(6) = 10^(-6); %gammaC
pMin(7) = 0.0015; %sigmaH
pMin(8) = 0.0005; %sigmaC
pMin(9) = 0.005; %lambdaH
pMin(10) = 0.0005; %lambdaC
pMin(11) = 0.2; %pG
pMin(12) = 5*10^(-10); %bH
pMin(13) = 5*10^(-5); %bZ
pMin(14) = 0.0005; %lambdaZ
pMin(15) = 10^(-8); %gammaZ
pMin(16) = 0.2; %pEx
pMin(17) = 0.005; %omega
pMin(18) = 2.5; %Nsub
pMin(19) = 10^(-13); %phiVH
pMin(20) = 10^(-13); %phiVC
pMin(21) = 0.5*10^(-12) ; %bV
pMin(22) = 0.05; %thetaH
pMin(23) = 0.05; %thetaC
pMin(24) = 15; %betaH
pMin(25) = 20; %betaC
pMin(26) = 5*10^(-10); %bC
pMin(27) = 0.05; %epsiH


% Define parameters' maximum
pMax = (1:26);

pMax(1) = 5; %rH
pMax(2) = 1.5; %rC
pMax(3) = 1; %kOn
pMax(4) = 1; %kIn
pMax(5) = 10^(-4); %gammaH
pMax(6) = 10^(-4); %gammaC
pMax(7) = 0.003; %sigmaH
pMax(8) = 0.0015; %sigmaC
pMax(9) = 0.05; %lambdaH
pMax(10) = 0.005; %lambdaC
pMax(11) = 0.4; %pG
pMax(12) = 4*10^(-9); %bH
pMax(13) = 4*10^(-4); %bZ
pMax(14) = 0.01; %lambdaZ
pMax(15) = 10^(-4); %gammaZ
pMax(16) = 0.4; %pEx
pMax(17) = 0.02; %omega
pMax(18) = 10; %Nsub
pMax(19) = 10^(-10); %phiVH
pMax(20) = 10^(-10); %phiVC
pMax(21) = 20*10^(-12) ; %bV
pMax(22) = 1; %thetaH
pMax(23) = 1; %thetaC
pMax(24) = 60; %betaH
pMax(25) = 250; %betaC
pMax(26) = 4*10^(-9); %bC
pMax(27) = 0.2; %epsiH

%Log transform min and max
logMin = log(pMin);
logMax = log(pMax);

% Find range spread in log space
diffLog = logMax - logMin;

% Number of data sets generated
reps = 1000;

% Specify seed for reproducibility
rng(2018);

% Generate 'reps'random sets of numbers, ranging from 0 to 1
randomUnscaled = lhsdesign(reps,27);

% Scale the random number with the range spread in log space
randomdiffLog = zeros(reps,27);
 
for n = 1:27
    randomdiffLog(:,n) = randomUnscaled(:,n)*diffLog(n);    
end

randomLogScaled = logMin + randomdiffLog;

%Get out of log space
randomScaled = exp(randomLogScaled);

%Store all results 
results2A = zeros(1,5);
results2B = zeros(1,4);

% Store positive results when coexistence
results2APos = zeros(1,5);
results2BPos = zeros(1,4);

coexist2A = 0; % Number of all positive equilibria for 7A
coexist2B = 0; % Number of all positive equilibria for 4B

comp = 0; % Keep track of # coexistence in both A and B 

for j = 1:reps
    
    p = randomScaled(j,:);
    
    % Evaluate equilibria for 2A
    equiV = Model2AEqui(p);
    
    % Store results
    results2A(j,:) = equiV;
    
    %Check for coexistence
    coexistVLogic = (sum(equiV > 0) == length(equiV));  
    coexist2A = coexist2A + coexistVLogic;
    
    % Evaluate equilibria for 2B
    equiNV = Model2BEqui(p);
    
    % Store results
    results2B(j,:) = equiNV;
    
    %Check for coexistence
    coexistNVLogic = (sum(equiNV > 0) == length(equiNV));  
    coexist2B = coexist2B + coexistNVLogic;
    
    % If both show coexistence, store results
    if coexistVLogic && coexistNVLogic
        
        comp = comp + 1;
        
        results2APos(comp,:) = equiV;
        results2BPos(comp,:) = equiNV;
        
    end
   
end

% Create graphs

sz = 7;

% Heterotrophic bacteria
figure;
subplot(2,3,1)
z = 3*10^8;
plot(z,z,'g^')
hold on
x = [1, 10^11];
y = [1, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results2BPos(:,1),results2APos(:,1),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^6 10^11])
ylim([10^6 10^11])
xticks([10^6 10^8 10^10])
yticks([10^6 10^8 10^10])
title("Heterotrophic bacteria (L^{-1})")    
    
% Zooplankton

subplot(2,3,3)
z = 5*10^4;
plot(z,z,'g^')
hold on
x = [1, 10^11];
y = [1, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results2BPos(:,2),results2APos(:,2),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^1 10^6])
ylim([10^1 10^6])
xticks([10^1 10^3 10^5])
yticks([10^1 10^3 10^5])
title("Zooplankton (L^{-1})")

% Inorganic nitrogen
subplot(2,3,4)
z = 0.1;
plot(z,z,'g^')
hold on
x = [10^-6, 10^2];
y = [10^-6, 10^2];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results2BPos(:,3),results2APos(:,3),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-4 10^2])
ylim([10^-4 10^2])
xticks([10^-4 10^-2 1 10^2])
yticks([10^-4 10^-2 1 10^2])
title("Inorganic nitrogen (\muM)")    
    
% Organic nitrogen
subplot(2,3,5)
z = 5;
plot(z,z,'g^')
hold on
x = [10^-6, 10^4];
y = [10^-6, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results2BPos(:,4),results2APos(:,4),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-3 10^4])
ylim([10^-3 10^4])
xticks([10^-3 10^-1 10^1 10^3])
yticks([10^-3 10^-1 10^1 10^3])
title("Organic nitrogen (\muM)")    
    