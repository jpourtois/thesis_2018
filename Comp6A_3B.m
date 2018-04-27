% Generate random sets of parameters to evaluate equilibrium of Model 6A and
% 3B

% Parameters' minimum
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

%Parameters' maximum
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
reps = 1000;

% Specify seed for reproducibility
rng(2018);

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
results6A = zeros(1,7);
results3B = zeros(1,4);

% Store positive results when coexistence
results6APos = zeros(1,7);
results3BPos = zeros(1,4);

coexist6A = 0; % Number of all positive equilibria for 6A
coexist3B = 0; % Number of all positive equilibria for 3B

comp3B = 0; % Number equilibria 3B (!= reps because quadratic formula)
comp = 0; % Keep track of # coexistence 

coexistV = 0;
coexistNV = 0; 

for j = 1:reps
    
    p = randomScaled(j,:);
    
    % Evaluate 6A equilibria
    equiV = Model6AEqui(p);
    
    % Store values
    results6A(j,:) = equiV;
    
    %Check for coexistence
    if (length(equiV) == sum(equiV > 0))
        coexist6A = coexist6A + 1;
    end
    
    % Evaluate 3B equilibria
    equiNVMult = zeros(1,5);
    equiNVMult = Model3BEqui(p);
    
    % Store results and check for coexistence
    
    for i = 1:length(equiNVMult(:,1))
        
        comp3B = comp3B + 1;
        equiNV = equiNVMult(i,:);
        results3B(comp3B,:) = equiNV;
        
        if (length(equiNV) == sum(equiNV > 0))
 
            coexist3B = coexist3B + 1;
            
            if (length(equiV) == sum(equiV > 0))
            
                comp = comp + 1;
                
                % Store positive results
                results6APos(comp,:) = equiV;
                results3BPos(comp,:) = equiNV;
                
                % Put back equilibrium values into equations
                 verif6A(comp,:) = Model6A_Verif(p,equiV);
                
                 % Calculate fluxes
                 
                    % DOM recycling
                    eco6A(comp,1) = equiV(3)*p(15)*(p(9)*p(5)*equiV(1)+p(28)*p(23)*equiV(2)) + (p(9)-p(17)*p(18))*p(6)*equiV(1)*equiV(6) + (p(28) - p(17)*p(25))*p(24)*equiV(2)*equiV(7);
                    eco3B(comp,1) = equiNV(2)*p(15)*p(9)*p(5)*equiNV(1);
                
                    % Export
                    eco6A(comp,2) = p(10)*p(12)*equiV(3)^2;
                    eco3B(comp,2) = p(10)*p(12)*equiNV(2)^2;
                
                    % Carbon sink
                    eco6A(comp,3) = p(33)*eco6A(comp,2)*p(35) + (p(9)*p(37)-p(17)*p(40)*p(18))*p(6)*equiV(1)*equiV(6)*p(38) + (p(28)*p(37) - p(17)*p(40)*p(25))*p(24)*equiV(2)*equiV(7)*p(38);
                    eco3B(comp,3) = p(33)*eco3B(comp,2)*p(35);
                    
                    parameters(comp,:) = p;
             
                
            end
  
        end
    end
                   
end

% Calculate total bacteria phosphorus content

bacPhosV = results6APos(:,1).*parameters(:,9) + results6APos(:,2).*parameters(:,28);
bacPhosNV = results3BPos(:,1).*parameters(:,9);

%%% Analysis %%%

logResults6A = log10([bacPhosV, results6APos(:,3:5)]);
logResults3B = log10([bacPhosNV, results3BPos(:,2:4)]);

logEco6A = log10(eco6A);
logEco3B = log10(eco3B);


for i = 1:length(logResults6A(1,:))
    
    %Check for normality
    [nonNormal6A(i),p6A(i)] = kstest(logResults6A(:,i));
    [nonNormal3B(i),p3B(i)] = kstest(logResults3B(:,i));
    
end

% Calculate median
median6A = median(logResults6A);
median3B = median(logResults3B);

medianEco6A = median(logEco6A);
medianEco3B = median(logEco3B);

% Calculate Interquatile range
iqr6A = iqr(logResults6A);
iqr3B = iqr(logResults3B);

iqrEco6A = iqr(logEco6A);
iqrEco3B = iqr(logEco3B);

% Calculate difference between 6A and 3B

diffLogResults6A_3B = logResults6A - logResults3B;

diffLogEco6A_3B = logEco6A - logEco3B;

for i = 1:length(logResults6A(1,:))
    
    %Check for normality
    [diff6A_3BNormal(i),pNormallDiff6A_3B(i)] = kstest(diffLogResults6A_3B(:,i));
    
    %Median different from 0?
    pDiff6A_3B(i) = signrank(diffLogResults6A_3B(:,i));
    
    
end

for j = 1:length(logEco6A(1,:))
    
    % Check for normality
    [diffEco6A_3BNormal(j),pNormalDiffEco6A_3B(j)] = kstest(diffLogEco6A_3B(:,j));

    %Median different from 0?
    pDiffEco6A_3B(j) = signrank(diffLogEco6A_3B(:,j));
    
end

% Find median and IQR
medianDiff6A_3B = median(diffLogResults6A_3B);
iqrDiff6A_3B = iqr(diffLogResults6A_3B);

medianDiffEco6A_3B = median(diffLogEco6A_3B);
iqrDiffeEco6A_3B = iqr(diffLogEco6A_3B);

%%% Graph

sz = 7;

% Phosphorus stored in bacteria
figure;
subplot(2,4,1)
z = 35*10^-3;
plot(z,z,'g^')
hold on
x = [10^-4, 10^11];
y = [10^-4, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(bacPhosNV,bacPhosV,sz,'black' )
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
title("Phosphorus in bacteria (\muM)")

% Zooplankton

subplot(2,4,2)
z = 5*10^4;
plot(z,z,'g^')
hold on
x = [10^-1, 10^11];
y = [10^-1, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results3BPos(:,2),results6APos(:,3),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-1 10^6])
ylim([10^-1 10^6])
xticks([10^-1 10^1 10^3 10^5])
yticks([10^-1 10^1 10^3 10^5])
title("Zooplankton (L^{-1})")

% Inorganic phosphorus
subplot(2,4,3)
z = 0.007;
plot(z,z,'g^')
hold on
x = [10^-6, 10^3];
y = [10^-6, 10^3];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results3BPos(:,3),results6APos(:,4),sz,'black' )
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
subplot(2,4,4)
z = 0.1;
plot(z,z,'g^')
hold on
x = [10^-6, 10^4];
y = [10^-6, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results3BPos(:,4),results6APos(:,5),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-5 10^4])
ylim([10^-5 10^4])
xticks([10^-5 10^-3 10^-1 10^1 10^3])
yticks([10^-5 10^-3 10^-1 10^1 10^3])
title("Organic phosphorus (\muM)")

% DOP recycling

subplot(2,4,5)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-6, 10^4];
y = [10^-6, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco3B(:,1),eco6A(:,1),sz,'black')
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
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

subplot(2,4,6)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-7, 10^4];
y = [10^-7, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco3B(:,2),eco6A(:,2),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-7 10^-1])
ylim([10^-7 10^-1])
xticks([10^-7 10^-5 10^-3 10^-1])
yticks([10^-7 10^-5 10^-3 10^-1])
title("Nutrient Export (\muM P d^{-1}) ")
set(gcf, 'Position', [800, 800, 450, 450])

% Carbon sink

subplot(2,4,7)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-13, 10^4];
y = [10^-13, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco3B(:,3),eco6A(:,3),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-5 10^3])
ylim([10^-5 10^3])
xticks([10^-7 10^-5 10^-3 10^-1 10^1 10^3])
yticks([10^-7 10^-5 10^-3 10^-1 10^1 10^3])
title("Carbon sink (\muM C d^{-1})")

set(gcf, 'Position', [800, 800, 950, 450])

saveas(gcf,'model6A_3B.png')