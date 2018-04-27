% Generate random sets of parameters to evaluate equilibrium of Model 8A and
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
results8A = zeros(1,7);
results5B = zeros(1,5);

% Store positive results when coexistence
results8APos = zeros(1,7);
results5BPos = zeros(1,5);

results8Aopti = zeros(1,7);
results5Bopti = zeros(1,5);

eco5B = zeros(1,4);
eco8A = zeros(1,4);

eco5Bopti = zeros(1,4);
eco8Aopti = zeros(1,4);


coexist8A = 0; % Number of all positive equilibria for 7A
coexist5B = 0; % Number of all positive equilibria for 4B

comp5B = 0; % Number equilibria 4B (!= reps because quadratic formula)
comp = 0; % Keep track of # coexistence 
comp_pos = 0; % For optimized sets

coexistV = 0;
coexistNV = 0; 

for j = 1:reps
    
    p = randomScaled(j,:);
    
    % Evaluate 8A equilibria
    equiV = Model8AEqui(p);
    
    % Store values
    results8A(j,:) = equiV;
    
    %Check for coexistence
    if (length(equiV) == sum(equiV > 0))
        coexist8A = coexist8A + 1;
    end
    
    % Evaluate 5B equilibria
    equiNVMult = zeros(1,5);
    equiNVMult = Model5BEqui(p);
    
    
    % Store results and check for coexistence
    for i = 1:length(equiNVMult(:,1))
        
        equiNV = equiNVMult(i,:);
        comp5B = comp5B + 1;
        results5B(comp5B,:) = equiNV;
        
        if (length(equiNV) == sum(equiNV > 0))
 
            coexist5B = coexist5B + 1;
            
            if (length(equiV) == sum(equiV > 0))
            
                comp = comp + 1;
                
                results8APos(comp,:) = equiV;
                results5BPos(comp,:) = equiNV;
                
                % DOM recycling
                eco8A(comp,1) = equiV(3)*p(15)*(p(9)*p(5)*equiV(1)+p(28)*p(23)*equiV(2)) + (p(9)-p(17)*p(18))*p(6)*equiV(1)*equiV(6) + (p(28) - p(17)*p(25))*p(24)*equiV(2)*equiV(7);
                eco5B(comp,1) = equiNV(3)*p(15)*(p(9)*p(5)*equiNV(1)+p(28)*p(23)*equiNV(2));
                
                % Export
                eco8A(comp,2) = p(10)*p(12)*equiV(3)^2;
                eco5B(comp,2) = p(10)*p(12)*equiNV(3)^2;
                
                % Carbon sink
                eco8A(comp,3) = p(33)*eco8A(comp,2)*p(35) + (p(9)*p(37)-p(17)*p(40)*p(18))*p(6)*equiV(1)*equiV(6)*p(38) + (p(28)*p(37) - p(17)*p(40)*p(25))*p(24)*equiV(2)*equiV(7)*p(38);
                eco5B(comp,3) = p(33)*eco5B(comp,2)*p(35);
                
                %Primary productivity
                eco8A(comp,4) = p(19)*equiV(2)*equiV(4)*p(28)/p(20);
                eco5B(comp,4) = p(19)*equiNV(2)*equiNV(4)*p(28)/p(20);
                
                parameters(comp,:) = p;
                
            end
  
        end
        
    end
    
    % Optimization procedure
    
    if coexist8A
               
        x0 = p(1:32);
        [pOptimized, fval] = fminunc(@Model8AEquiOpti,x0);
        
        equiVopti = Model8AEqui(pOptimized);
        
        if(sum(equiVopti > 0) == length(equiVopti))
           
            equiNVopti = Model5BEqui(pOptimized);
            
            for g = 1:length(equiNVopti(:,1))
                
                if (sum(equiNVopti(g,:) > 0) == length(equiNVopti(g,:)))                    
                  
                    p_bigMin = pOptimized >= pMin(1:32);
                    p_lowMax = pOptimized <= pMax(1:32);                  
                    
                    if (sum(p_bigMin) == length(p_bigMin))
                        
                        if (sum(p_lowMax) == length(p_lowMax))
                            
                            
                            pNew(1:32) = pOptimized; 
                            pNew(33:40) = p(33:40);
                            comp_pos = comp_pos + 1;
                            p_positive(:,comp_pos) = pOptimized;
                            results8Aopti(comp_pos,:) = equiVopti;
                            results5Bopti(comp_pos,:) = equiNVopti(g,:);
                            valueF(comp_pos) = fval;
                            
                            % DOM recycling
                            eco8Aopti(comp,1) = equiVopti(3)*p(15)*(p(9)*p(5)*equiVopti(1)+p(28)*p(23)*equiVopti(2)) + (p(9)-p(17)*p(18))*p(6)*equiVopti(1)*equiVopti(6) + (p(28) - p(17)*p(25))*p(24)*equiVopti(2)*equiVopti(7);
                            eco5Bopti(comp,1) = equiNVopti(g,3)*p(15)*(p(9)*p(5)*equiNVopti(g,1)+p(28)*p(23)*equiNVopti(g,2));
                
                            % Export
                            eco8Aopti(comp,2) = p(10)*p(12)*equiVopti(3)^2;
                            eco5Bopti(comp,2) = p(10)*p(12)*equiNVopti(g,3)^2;
                
                            % Carbon sink
                            eco8Aopti(comp,3) = p(33)*eco8Aopti(comp,2)*p(35) + (p(9)*p(37)-p(17)*p(40)*p(18))*p(6)*equiVopti(1)*equiVopti(6)*p(38) + (p(28)*p(37) - p(17)*p(40)*p(25))*p(24)*equiVopti(2)*equiVopti(7)*p(38);
                            eco5Bopti(comp,3) = p(33)*eco5Bopti(comp,2)*p(35);
                
                            %Primary productivity
                            eco8Aopti(comp,4) = p(19)*equiVopti(2)*equiVopti(4)*p(28)/p(20);
                            eco5Bopti(comp,4) = p(19)*equiNVopti(g,2)*equiNVopti(g,4)*p(28)/p(20);
                            
                        end
                    end
                end
            end
        end
    end              
end

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
scatter(results5BPos(:,1),results8APos(:,1),sz,'black' )
scatter(results5Bopti(:,1),results8Aopti(:,1),sz,'blue' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^5 10^12])
ylim([10^5 10^12])
xticks([10^5 10^7 10^9 10^11])
yticks([10^5 10^7 10^9 10^11])
title("Heterotrophic bacteria (L^{-1})")

% Cyanobacteria

subplot(3,3,2)
z = 3*10^8;
plot(z,z,'g^')
hold on
x = [1, 10^13];
y = [1, 10^13];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,2),results8APos(:,2),sz,'black' )
scatter(results5Bopti(:,2),results8Aopti(:,2),sz,'blue' )
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
scatter(results5BPos(:,3),results8APos(:,3),sz,'black' )
scatter(results5Bopti(:,3),results8Aopti(:,3),sz,'blue' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([1 10^6])
ylim([1 10^6])
xticks([1 10^2 10^4 10^6])
yticks([1 10^2 10^4 10^6])
title("Zooplankton (L^{-1})")

% Inorganic phosphorus
subplot(3,3,4)
z = 0.007;
plot(z,z,'g^')
hold on
x = [10^-6, 10^2];
y = [10^-6, 10^2];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,4),results8APos(:,4),sz,'black' )
scatter(results5Bopti(:,4),results8Aopti(:,4),sz,'blue' )
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
scatter(results5BPos(:,5),results8APos(:,5),sz,'black' )
scatter(results5Bopti(:,5),results8Aopti(:,5),sz,'blue' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-7 10^2])
ylim([10^-7 10^2])
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
scatter(eco5B(:,4),eco8A(:,4),sz,'black')
scatter(eco5Bopti(:,4),eco8Aopti(:,4),sz,'blue')
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-7 10^1])
ylim([10^-7 10^1])
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
scatter(eco5B(:,1),eco8A(:,1),sz,'black')
scatter(eco5Bopti(:,1),eco8Aopti(:,1),sz,'blue')
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
scatter(eco5B(:,2),eco8A(:,2),sz,'black' )
scatter(eco5Bopti(:,2),eco8Aopti(:,2),sz,'blue')
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
scatter(eco5B(:,3),eco8A(:,3),sz,'black' )
scatter(eco5Bopti(:,3),eco8Aopti(:,3),sz,'blue')
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

saveas(gcf,'model8A_5B.png')