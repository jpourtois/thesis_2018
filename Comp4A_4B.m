% Generate random sets of parameters to evaluate equilibrium of Model 4A and
% 4B

% Define parameters' minimum

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
pMin(28) = 0.05; %pF
pMin(29) = 4; %zCN
pMin(30) = 60; %zCP
pMin(31) = 4; %bCN
pMin(32) = 64; %bCP
pMin(33) = 0.1; %pRef
pMin(34) = 3; % vCN
pMin(35) = 18; % vCP


% Define parameters' maximum

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
pMax(21) = 3*10^(-12) ; %bV
pMax(22) = 1; %thetaH
pMax(23) = 1; %thetaC
pMax(24) = 50; %betaH
pMax(25) = 150; %betaC
pMax(26) = 4*10^(-9); %bC
pMax(27) = 0.2; %epsiH
pMax(28) = 0.3; %pF
pMax(29) = 8; %zCN
pMax(30) = 200; %zCP
pMax(31) = 7; %bCN
pMax(32) = 112; %bCP
pMax(33) = 0.3; %pRef
pMax(34) = 3.5; % vCN
pMax(35) = 45; % vCP


%Log transform min and max
logMin = log(pMin);
logMax = log(pMax);

% Find range spread in log space
diffLog = logMax - logMin;

% Number of data sets generated
reps = 5000;

% Specify seed for reproducibility
rng(2018);

% Generate 'reps'random sets of numbers, ranging from 0 to 1
randomUnscaled = lhsdesign(reps,35);

% Scale the random number with the range spread in log space
randomdiffLog = zeros(reps,35);
 
for n = 1:35
    randomdiffLog(:,n) = randomUnscaled(:,n)*diffLog(n);    
end

randomLogScaled = logMin + randomdiffLog;

%Get out of log space
randomScaled = exp(randomLogScaled);

%Store all results 
results4A = zeros(1,6);
results4B = zeros(1,5);

% Store positive results when coexistence
results4APos = zeros(1,6);
results4BPos = zeros(1,5);

coexist4A = 0; % Number of all positive equilibria for 7A
coexist4B = 0; % Number of all positive equilibria for 4B

comp4A = 0; % Number equilibria 4A
comp4B = 0; % Number equilibria 4B (!= reps because quadratic formula)
comp = 0; % Keep track of # coexistence 

for j = 1:reps
    
    p = randomScaled(j,:);
  
    %Get equilibrium values for 4A
    equiVMult = 0;
    equiVMult = Model4AEqui(p);
    
    for k = 1:length(equiVMult(:,1))
    
        comp4A = comp4A + 1;
        
        equiV = equiVMult(k,:);
        
        %Store results for 4A
        results4A(comp4A,:) = equiV;
    
        %Check for coexistence
        coexistVLogic = (sum(equiV > 0) == length(equiV));
    
        coexist4A = coexist4A + coexistVLogic;
    
        %Get equilibria for 4B
        equiNVMult = zeros(1,5);
        equiNVMult = Model4BEqui(p);
    
        for i = 1:length(equiNVMult(:,1))
         
            % Store equilibrium values for 4B
            comp4B = comp4B + 1;
        
            equiNV = equiNVMult(i,:);
            results4B(comp4B,:) = equiNV;
            coexistNVLogic = (sum(equiNV > 0) == length(equiNV));   
            coexist4B = coexist4B + coexistNVLogic;  
            
            if coexistVLogic & coexistNVLogic
                
                comp = comp + 1;
                 
                % Store positive results
                results4APos(comp,:) = equiV;
                results4BPos(comp,:) = equiNV;
                 
                % Put back equilibrium values into equations
                verif4A(comp,:) = Model4A_Verif(p,equiV);
                % See Comp7A_4B for verif4B 
                
                % Calculate fluxes
                
                % DOM recycling
                eco4A(comp,1) = p(16)*equiV(3)*(p(12)*p(5)*equiV(1) + p(26)*p(6)*equiV(2)) + (p(12)-p(21)*p(24))*p(19)*equiV(1)*equiV(6); 
                eco4B(comp,1) = p(16)*equiNV(3)*(p(12)*p(5)*equiNV(1) + p(26)*p(6)*equiNV(2));
                
                % Export
                eco4A(comp,2) = p(13)*p(15)*equiV(3)^2;
                eco4B(comp,2) = p(13)*p(15)*equiNV(3)^2;
                
                % Carbon sink
                eco4A(comp,3) = p(28)*eco4A(comp,2)*p(29) + (p(12)*p(31) - p(21)*p(34)*p(24))*p(19)*equiV(1)*equiV(6)*p(33);
                eco4B(comp,3) = p(28)*eco4B(comp,2)*p(29);
                
                % Primary Productivity
                eco4A(comp,4) = p(2)*equiV(2)*equiV(4)*p(26)/p(4);
                eco4B(comp,4) = p(2)*equiNV(2)*equiNV(4)*p(26)/p(4);
            end
        end
    end   
end 

sz = 7;

%%% Analysis %%%

logResults4A = log10(results4APos(:,1:5));
logResults4B = log10(results4BPos);

logEco4A = log10(eco4A);
logEco4B = log10(eco4B);


for i = 1:length(logResults4A(1,:))
    
    %Check for normality
    [nonNormal4A(i),p4A(i)] = kstest(logResults4A(:,i));
    [nonNormal4B(i),p4B(i)] = kstest(logResults4B(:,i));
    
end

% Calculate median
median4A = median(logResults4A);
median4B = median(logResults4B);

medianEco4A = median(logEco4A);
medianEco4B = median(logEco4B);

% Calculate Interquatile range
iqr4A = iqr(logResults4A);
iqr4B = iqr(logResults4B);

iqrEco4A = iqr(logEco4A);
iqrEco4B = iqr(logEco4B);

% Calculate difference between 8A and 5B

diffLogResults4A_4B = logResults4A - logResults4B;

diffLogEco4A_4B = logEco4A - logEco4B;

for i = 1:length(logResults4A(1,:))
    
    %Check for normality
    [diff4A_4BNormal(i),pNormallDiff4A_4B(i)] = kstest(diffLogResults4A_4B(:,i));
    
    %Median different from 0?
    pDiff4A_4B(i) = signrank(diffLogResults4A_4B(:,i));
    
    
end

for j = 1:length(logEco4A(1,:))
    
    % Check for normality
    [diffEco4A_4BNormal(j),pNormalDiffEco4A_4B(j)] = kstest(diffLogEco4A_4B(:,j));

    %Median different from 0?
    pDiffEco4A_4B(j) = signrank(diffLogEco4A_4B(:,j));
    
end

% Find median and IQR
medianDiff4A_4B = median(diffLogResults4A_4B);
iqrDiff4A_4B = iqr(diffLogResults4A_4B);

medianDiffEco4A_4B = median(diffLogEco4A_4B);
iqrDiffeEco4A_4B = iqr(diffLogEco4A_4B);

%%% Graph

% Heterotrophic bacteria
figure;
subplot(3,3,1)
z = 3*10^8;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [1, 10^11];
y = [1, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results4BPos(:,1),results4APos(:,1),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^6 10^11])
ylim([10^6 10^11])
xticks([10^6 10^8 10^10])
yticks([10^6 10^8 10^10])
title("Heterotrophic bacteria (L^{-1})")

% Cyanobacteria

subplot(3,3,2)
z = 3*10^8;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [1, 10^13];
y = [1, 10^13];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results4BPos(:,2),results4APos(:,2),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^6 10^12])
ylim([10^6 10^12])
xticks([10^6 10^8 10^10 10^12])
yticks([10^6 10^8 10^10 10^12])
title("Cyanobacteria (L^{-1})")

% Zooplankton

subplot(3,3,3)
z = 5*10^4;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [1, 10^11];
y = [1, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results4BPos(:,3),results4APos(:,3),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
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
subplot(3,3,4)
z = 0.1;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-6, 10^2];
y = [10^-6, 10^2];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results4BPos(:,4),results4APos(:,4),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
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
subplot(3,3,5)
z = 5;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-6, 10^4];
y = [10^-6, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results4BPos(:,5),results4APos(:,5),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-4 10^2])
ylim([10^-4 10^2])
xticks([10^-4 10^-2 10^0 10^2])
yticks([10^-4 10^-2 10^0 10^2])
title("Organic nitrogen (\muM)")

% Primary Productivity

subplot(3,3,6)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-6, 10^4];
y = [10^-6, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco4B(:,4),eco4A(:,4),sz,'black')
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-3 10^2])
ylim([10^-3 10^2])
xticks([10^-5 10^-3 10^-1 10^1 10^3])
yticks([10^-5 10^-3 10^-1 10^1 10^3])
title("Primary Productivity (\muM N d^{-1})")

% DOP recycling

subplot(3,3,7)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-6, 10^4];
y = [10^-6, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco4B(:,1),eco4A(:,1),sz,'black')
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-3 10^1])
ylim([10^-3 10^1])
xticks([10^-5 10^-3 10^-1 10^1 10^3])
yticks([10^-5 10^-3 10^-1 10^1 10^3])
title("DON recycling (\muM d^{-1})")


% Nutrient export

subplot(3,3,8)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-7, 10^4];
y = [10^-7, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco4B(:,2),eco4A(:,2),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-3 10])
ylim([10^-3 10])
xticks([10^-3 10^-1 10])
yticks([10^-3 10^-1 10])
title("Nutrient Export (\muM N d^{-1}) ")
set(gcf, 'Position', [800, 800, 450, 450])

% Carbon sink

subplot(3,3,9)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-13, 10^4];
y = [10^-13, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco4B(:,3),eco4A(:,3),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-3 10^1])
ylim([10^-3 10^1])
xticks([10^-7 10^-5 10^-3 10^-1 10^1 10^3])
yticks([10^-7 10^-5 10^-3 10^-1 10^1 10^3])
title("Carbon sink (\muM C d^{-1})")
set(gcf, 'Position', [800, 800, 950, 450])


set(gcf, 'Position', [800, 800, 800, 800])

saveas(gcf,'model4A_4B.png')
