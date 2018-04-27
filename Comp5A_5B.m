% Generate random sets of parameters to evaluate equilibrium of Model 5A and
% 5B

% Define parameters' minimum

pMin(1) = 0.5; %rH 
pMin(2) = 0.005; %kInH 
pMin(3) = 0.005; %lambdaH
pMin(4) = 0.0015; %sigmaH
pMin(5) = 10^(-6); %gammaH
pMin(6) = 10^(-13); %phiH 
pMin(7) = 0.05; %thetaH
pMin(8) = 0.2; %pG
pMin(9) = 10^(-11); %qH 
pMin(10) = 10^(-6); %qZ
pMin(11) = 0.0005; %lambdaZ
pMin(12) = 10^(-6); %gammaZ
pMin(13) = 0.005; %omega
pMin(14) = 0.05; %pSub
pMin(15) = 0.2; %pEx
pMin(16) = 2*10^(-13); %psiH
pMin(17) = 0.5*10^(-13); %qV
pMin(18) = 15; %betaH
pMin(19) = 0.1; % rC
pMin(20) = 0.005; % kInC
pMin(21) = 0.0005; % lambdaC
pMin(22) = 0.0005; % sigmaC
pMin(23) = 10^(-6); % gammaC
pMin(24) = 10^(-13); % phiC
pMin(25) = 20; % betaC
pMin(26) = 0.05; % thetaC
pMin(27) = 2*10^(-13); % psiC
pMin(28) = 3*10^(-11); %qC 
pMin(29) = 0.005; %kOrgH
pMin(30) = 0.005; %kOrgC
pMin(31) = 0.05; % epsiH
pMin(32) = 0.05; % epsiC
pMin(33) = 0.05; %pF
pMin(34) = 4; %zCN
pMin(35) = 60; %zCP
pMin(36) = 4; %bCN
pMin(37) = 64; %bCP
pMin(38) = 0.1; %pRef
pMin(39) = 3; % vCN
pMin(40) = 18; % vCP

% Define parameters' maximum

pMax(1) = 5;%rH
pMax(2) = 1.5;%kInH
pMax(3) = 0.05;%lambdaH
pMax(4) = 0.003;%sigmaH
pMax(5) = 10^(-4);%gammaH
pMax(6) = 10^(-10);%phiH
pMax(7) = 1;%thetaH
pMax(8) = 0.4;%pG
pMax(9) = 5*10^(-10);%qH
pMax(10) = 2*10^(-5);%qZ
pMax(11) = 0.01;%lambdaZ
pMax(12) = 10^(-3);%gammaZ
pMax(13) = 0.02;%omega
pMax(14) = 0.2;%pSub
pMax(15) = 0.4;%pEx
pMax(16) = 2*10^(-10); %psiH
pMax(17) = 2*10^(-13); %qV
pMax(18) = 50; %betaH
pMax(19) = 1.5; % rC
pMax(20) = 1.5; % kInC
pMax(21) = 0.005; % lambdaC
pMax(22) = 0.0015; % sigmaC
pMax(23) = 10^(-4); % gammaC
pMax(24) = 10^(-10); % phiC
pMax(25) = 150; % betaC
pMax(26) = 1; % thetaC
pMax(27) = 2*10^(-10); % psiC
pMax(28) = 2*10^(-10); %qC
pMax(29) = 2.5; %kOrgH
pMax(30) = 2.5; %kOrgC
pMax(31) = 0.2; % epsiH
pMax(32) = 0.2; % epsiC
pMax(33) = 0.3; %pF
pMax(34) = 8; %zCN
pMax(35) = 200; %zCP
pMax(36) = 7; %bCN
pMax(37) = 112; %bCP
pMax(38) = 0.3; %pRef
pMax(39) = 3.5; % vCN
pMax(40) = 45; % vCP

%Log transform min and max
logMin = log(pMin);
logMax = log(pMax);

% Find range spread in log space
diffLog = logMax - logMin;

% Number of data sets generated
reps = 10000;

% Specify seed for reproducibility
rng(2016);

% Generate 'reps'random sets of numbers, ranging from 0 to 1
randomUnscaled = lhsdesign(reps,40);

% Scale the random number with the range spread in log space
randomdiffLog = zeros(reps,40);
 
for n = 1:40
    randomdiffLog(:,n) = randomUnscaled(:,n)*diffLog(n);    
end

randomLogScaled = logMin + randomdiffLog;

%Get out of log space
randomScaled = exp(randomLogScaled);

%Store all results 
results5A = zeros(1,6);
results5B = zeros(1,5);

% Store positive results when coexistence
results5APos = zeros(1,6);
results5BPos = zeros(1,5);

coexist5A = 0; % Number of all positive equilibria for 5A
coexist5B = 0; % Number of all positive equilibria for 5B

comp5A = 0; % Number equilibria 5A (!= reps because quadratic formula)
comp5B = 0; % Number equilibria 5B
comp = 0; % Keep track of # coexistence 

for j = 1:reps
    
    p = randomScaled(j,:);
    
    % Evaluate 5A equilibria
    
    equiVMult = zeros(1,6);
    equiVMult = Model5AEqui(p);
    
    for i = 1:length(equiVMult(:,1))
        
        comp5A = comp5A + 1;
        equiV = equiVMult(i,:);
        
        % Store values
        results5A(comp5A,:) = equiV;
    
        %Check for coexistence
        if (length(equiV) == sum(equiV > 0))
            coexist5A = coexist5A + 1;
        end
    
        % Evaluate 5B equilibria
        equiNVMult = zeros(1,5);
        equiNVMult = Model5BEqui(p);
    
        % Store results and check for coexistence
        for k = 1:length(equiNVMult(:,1))
        
            comp5B = comp5B + 1;
        
            equiNV = equiNVMult(k,:);
            results5B(comp5B,:) = equiNV;
        
            %Check for coexistence
            if (length(equiNV) == sum(equiNV > 0))
 
                coexist5B = coexist5B + 1;
            
                if (length(equiV) == sum(equiV > 0))
            
                    comp = comp + 1;
                    
                    % Store positive results
                    results5APos(comp,:) = equiV;
                    results5BPos(comp,:) = equiNV;
                    parameters(comp,:) = p;
                    
                    %Put back equilibrium values into equations
                    verif5A(comp,:) = Model5A_Verif(p,equiV);                  
                    % See Comp8A_5B for verif5B 
                    
                    %%% Evaluate fluxes
                    
                    % DOM recycling
                    eco5A(comp,1) = equiV(3)*p(15)*(p(9)*p(5)*equiV(1)+p(28)*p(23)*equiV(2)) + (p(9)-p(17)*p(18))*p(6)*equiV(1)*equiV(6);
                    eco5B(comp,1) = equiNV(3)*p(15)*(p(9)*p(5)*equiNV(1)+p(28)*p(23)*equiNV(2));
                
                    % Export
                    eco5A(comp,2) = p(10)*p(12)*equiV(3)^2;
                    eco5B(comp,2) = p(10)*p(12)*equiNV(3)^2;
                
                    % Carbon sink
                    eco5A(comp,3) = p(33)*eco5A(comp,2)*p(35) + (p(9)*p(37)-p(17)*p(40)*p(18))*p(6)*equiV(1)*equiV(6)*p(38);
                    eco5B(comp,3) = p(33)*eco5B(comp,2)*p(35);
                    
                    %Primary productivity
                    eco5A(comp,4) = p(19)*equiV(2)*equiV(4)*p(28)/p(20);
                    eco5B(comp,4) = p(19)*equiNV(2)*equiNV(4)*p(28)/p(20);
                end
            end
  
        end
    end
                   
end

%%% Analysis %%%

logResults5A = log10(results5APos(:,1:5));
logResults5B = log10(results5BPos);

logEco5A = log10(eco5A);
logEco5B = log10(eco5B);


for i = 1:length(logResults5A(1,:))
    
    %Check for normality
    [nonNormal5A(i),p5A(i)] = kstest(logResults5A(:,i));
    [nonNormal5B(i),p5B(i)] = kstest(logResults5B(:,i));
    
end

% Calculate median
median5A = median(logResults5A);
median5B = median(logResults5B);

medianEco5A = median(logEco5A);
medianEco5B = median(logEco5B);

% Calculate Interquatile range
iqr5A = iqr(logResults5A);
iqr5B = iqr(logResults5B);

iqrEco5A = iqr(logEco5A);
iqrEco5B = iqr(logEco5B);

% Calculate difference between 5A and 5B

diffLogResults5A_5B = logResults5A - logResults5B;

diffLogEco5A_5B = logEco5A - logEco5B;

for i = 1:length(logResults5A(1,:))
    
    %Check for normality
    [diff5A_5BNormal(i),pNormallDiff5A_5B(i)] = kstest(diffLogResults5A_5B(:,i));
    
    %Median different from 0?
    pDiff5A_5B(i) = signrank(diffLogResults5A_5B(:,i));
    
    
end

for j = 1:length(logEco5A(1,:))
    
    % Check for normality
    [diffEco5A_5BNormal(j),pNormalDiffEco5A_5B(j)] = kstest(diffLogEco5A_5B(:,j));

    %Median different from 0?
    pDiffEco5A_5B(j) = signrank(diffLogEco5A_5B(:,j));
    
end

% Find median and IQR
medianDiff5A_5B = median(diffLogResults5A_5B);
iqrDiff5A_5B = iqr(diffLogResults5A_5B);

medianDiffEco5A_5B = median(diffLogEco5A_5B);
iqrDiffeEco5A_5B = iqr(diffLogEco5A_5B);

%%% Graph

sz = 7;

% Heterotrophic bacteria
figure;
subplot(3,3,1)
z = 3*10^8;
plot(z,z,'g^')
hold on
x = [1, 10^12];
y = [1, 10^12];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,1),results5APos(:,1),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^7 10^12])
ylim([10^7 10^12])
xticks([10^7 10^9 10^11])
yticks([10^7 10^9 10^11])
title("Heterotrophic bacteria (L^{-1})")

% Cyanobacteria

subplot(3,3,2)
z = 3*10^8;
plot(z,z,'g^')
hold on
x = [1, 10^13];
y = [1, 10^13];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,2),results5APos(:,2),sz,'black' )
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
plot(z,z,'g^')
hold on
x = [1, 10^11];
y = [1, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,3),results5APos(:,3),sz,'black' )
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

% Inorganic phosphorus
subplot(3,3,4)
z = 0.007;
plot(z,z,'g^')
hold on
x = [10^-6, 10^2];
y = [10^-6, 10^2];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,4),results5APos(:,4),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-5 10^2])
ylim([10^-5 10^2])
xticks([10^-4 10^-2 1 10^2])
yticks([10^-4 10^-2 1 10^2])
title("Inorganic phosphorus (\muM)")

% Organic phosphorus
subplot(3,3,5)
z = 0.1;
plot(z,z,'g^')
hold on
x = [10^-7, 10^4];
y = [10^-7, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,5),results5APos(:,5),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-7 10^1])
ylim([10^-7 10^1])
xticks([10^-7 10^-5 10^-3 10^-1 10^1])
yticks([10^-7 10^-5 10^-3 10^-1 10^1])
title("Organic phosphorus (\muM)")

% Primary Productivity

subplot(3,3,6)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-7, 10^4];
y = [10^-7, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco5B(:,4),eco5A(:,4),sz,'black')
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-7 10^-1])
ylim([10^-7 10^-1])
xticks([10^-7 10^-5 10^-3 10^-1])
yticks([10^-7 10^-5 10^-3 10^-1])
title("Primary Productivity (\muM P d^{-1})")

% DOP recycling

subplot(3,3,7)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-6, 10^4];
y = [10^-6, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco5B(:,1),eco5A(:,1),sz,'black')
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-5 10^2])
ylim([10^-5 10^2])
xticks([10^-5 10^-3 10^-1 10^1 10^3])
yticks([10^-5 10^-3 10^-1 10^1 10^3])
title("DOP recycling (\muM d^{-1})")


% Nutrient export

subplot(3,3,8)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-7, 10^4];
y = [10^-7, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco5B(:,2),eco5A(:,2),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-5 10^-1])
ylim([10^-5 10^-1])
xticks([10^-7 10^-5 10^-3 10^-1])
yticks([10^-7 10^-5 10^-3 10^-1])
title("Nutrient Export (\muM P d^{-1}) ")
set(gcf, 'Position', [800, 800, 450, 450])

% Carbon sink

subplot(3,3,9)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-13, 10^4];
y = [10^-13, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco5B(:,3),eco5A(:,3),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-5 10^3])
ylim([10^-5 10^3])
xticks([10^-7 10^-5 10^-3 10^-1 10^1 10^3])
yticks([10^-7 10^-5 10^-3 10^-1 10^1 10^3])
title("Carbon sink (\muM C d^{-1})")

set(gcf, 'Position', [800, 800, 800, 800])

saveas(gcf,'model5A_5B.png')