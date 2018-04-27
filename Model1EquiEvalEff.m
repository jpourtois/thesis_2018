% Generate random sets of parameters to evaluate equilibrium of Model 1A,
% 1B and 1C (with growth efficiency). Create graphs to visualize and compare results.

% Define parameters' minimum
pMin = (1:9);

pMin(1) = 0.5; %rHmax 
pMin(2) = 0.005; %kInH 
pMin(3) = 0.005; %kOrgH
pMin(4) = 0.005; %lambdaH
pMin(5) = 0.0015; %sigmaH
pMin(6) = 10^(-11); %qH 
pMin(7) = 2*10^(-13); %psiH
pMin(8) = 0.08; %N
pMin(9) = 0.05; %epsiH


% Define parameters' maximum
pMax = (1:9);

pMax(1) = 5;%rHmax 
pMax(2) = 1.5;%kInH
pMax(3) = 2.5; %kOrgH
pMax(4) = 0.05;%lambdaH
pMax(5) = 0.003;%sigmaH
pMax(6) = 5*10^(-10);%qH
pMax(7) = 200*10^(-12); %psiH
pMax(8) = 0.18; %N
pMax(9) = 0.2; %epsiH

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
randomUnscaled = lhsdesign(reps,9);

% Scale the random number with the range spread in log space
randomdiffLog = zeros(reps,9);
 
for n = 1:9
    randomdiffLog(:,n) = randomUnscaled(:,n)*diffLog(n);    
end

randomLogScaled = logMin + randomdiffLog;

%Get out of log space
randomScaled = exp(randomLogScaled);

% Store positive results
resultsProp1Pos = zeros(1,3);
resultsProp1bPos = zeros(1,3);
resultsProp1cPos = zeros(1,3);

comp = 0; % Keep track of # coexistence 

% For each set of parameter, evaluate equilibrium for Model 1A, 1B and 1C
for j = 1:reps
    
    p = randomScaled(j,:);
    
    % Call Model 1A
    equi1a = Model1AEqui(p);
    
    % Check for coexistence
    if sum(equi1a > 0) == 3
        
        % Call Model 1B
        equi1b = Model1bEquiEff(p);
               
        %Check for coexistence
        if sum(equi1b > 0) == 3
                      
            % Call Model 1C
            equi1c = Model1cEquiEff(p);
         
            % Check for coexistence
            if sum(equi1c > 0) == 3
                
            comp = comp + 1;
            
            % Store results 
            resultsProp1Pos(comp,:) = equi1a;
            resultsProp1bPos(comp,:) = equi1b;
            resultsProp1cPos(comp,:) = equi1c;
            
            
            end            
        end       
    end
   
end 

% Define point size for graphs
sz = 7;

% Create graphs comparing Model 1A and 1B

% Heterotrophic bacteria
figure;
subplot(1,3,1)
z = 6*10^8;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [1, 10^11];
y = [1, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(resultsProp1Pos(:,1),resultsProp1bPos(:,1),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("Phosphatase")
ylabel("DOP uptake")
xlim([10^6 10^11])
ylim([10^6 10^11])
xticks([10^6 10^8 10^10])
yticks([10^6 10^8 10^10])
title("Heterotrophic bacteria (L^{-1})")

% Inorganic phosphorus
subplot(1,3,2)
z = 7*10^-3;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-6, 10^1];
y = [10^-6, 10^1];
line(x,y,'Color','red','LineWidth',1.5)
scatter(resultsProp1Pos(:,2),resultsProp1bPos(:,2),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("Phosphatase")
ylabel("DOP uptake")
xlim([10^-6 10])
ylim([10^-6 10])
xticks([10^-6 10^-4 10^-2 1])
yticks([10^-6 10^-4 10^-2 1])
title("Inorganic phosphorus (\muM)")

% Organic phosphorus
subplot(1,3,3)
z = 0.1;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-7, 10^1];
y = [10^-7, 10^1];
line(x,y,'Color','red','LineWidth',1.5)
scatter(resultsProp1Pos(:,3),resultsProp1bPos(:,3),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("Phosphatase")
ylabel("DOP uptake")
xlim([10^-7 10])
ylim([10^-7 10])
xticks([10^-7 10^-5 10^-3 10^-1])
yticks([10^-7 10^-5 10^-3 10^-1])
title("Organic phosphorus (\muM)")

set(gcf, 'Position', [800, 800, 700, 250])

saveas(gcf,'model1A_1B.png')

% Create graphs comparing Model 1B and 1C

% Heterotrophic bacteria
figure;
subplot(1,3,1)
z = 6*10^8;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
pbaspect([1 1 1])
x = [1, 10^11];
y = [1, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(resultsProp1bPos(:,1),resultsProp1cPos(:,1),sz,'black')
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("Type II")
ylabel("Type I")
xlim([10^6 10^11])
ylim([10^6 10^11])
xticks([10^6 10^8 10^10])
yticks([10^6 10^8 10^10])
title("Heterotrophic bacteria (L^{-1})")

% Inorganic phosphorus
subplot(1,3,2)
z = 7*10^-3;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-6, 10^1];
y = [10^-6, 10^1];
line(x,y,'Color','red','LineWidth',1.5)
scatter(resultsProp1bPos(:,2),resultsProp1cPos(:,2),sz,'black')
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("Type II")
ylabel("Type I")
xlim([10^-6 10])
ylim([10^-6 10])
xticks([10^-6 10^-4 10^-2 1])
yticks([10^-6 10^-4 10^-2 1])
title("Inorganic phosphorus (\muM)")

% Organic phosphorus
subplot(1,3,3)
z = 0.1;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-7, 10^1];
y = [10^-7, 10^1];
line(x,y,'Color','red','LineWidth',1.5)
scatter(resultsProp1bPos(:,3),resultsProp1cPos(:,3),sz,'black')
pbaspect([1 1 1])
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("Type II")
ylabel("Type I")
xlim([10^-7 1])
ylim([10^-7 1])
xticks([10^-7 10^-5 10^-3 10^-1])
yticks([10^-7 10^-5 10^-3 10^-1])
title("Organic phosphorus (\muM)")

set(gcf, 'Position', [800, 800, 700, 250])

saveas(gcf,'model1B_1C.png')
