% Make graph comparing Model 7 and 8
% FIRST RUN 'Comp7A_4B.m' and 'Comp8A_5B.m'

%%%% 1st Graph %%%%

%%% Nitrogen-limited system

sz = 7;

% Heterotrophic bacteria
figure;
subplot(3,2,1)
z = 6*10^8;
plot(z,z,'g^')
hold on
x = [1, 10^11];
y = [1, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results4BPos(:,1),results7APos(:,1),sz,'black' )
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

subplot(3,2,3)
z = 1*10^8;
plot(z,z,'g^')
hold on
x = [1, 10^13];
y = [1, 10^13];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results4BPos(:,2),results7APos(:,2),sz,'black' )
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

subplot(3,2,5)
z = 4*10^4;
plot(z,z,'g^')
hold on
x = [1, 10^11];
y = [1, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results4BPos(:,3),results7APos(:,3),sz,'black' )
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

%%% Phosphorus-limited system

% Heterotrophic bacteria

subplot(3,2,2)
z = 3*10^8;
plot(z,z,'g^')
hold on
x = [1, 10^12];
y = [1, 10^12];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,1),results8APos(:,1),sz,'black' )
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

subplot(3,2,4)
z = 3*10^8;
plot(z,z,'g^')
hold on
x = [1, 10^13];
y = [1, 10^13];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,2),results8APos(:,2),sz,'black' )
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

subplot(3,2,6)
z = 5*10^4;
plot(z,z,'g^')
hold on
x = [1, 10^11];
y = [1, 10^11];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,3),results8APos(:,3),sz,'black' )
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

set(gcf, 'Position', [800, 800, 600, 1000])
saveas(gcf,'poster1.png')

%%% 2nd graph 

%%% Nitrogen-limited system

% Inorganic nitrogen
figure;
subplot(3,2,1)
z = 0.1;
plot(z,z,'g^')
hold on
x = [10^-6, 10^2];
y = [10^-6, 10^2];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results4BPos(:,4),results7APos(:,4),sz,'black' )
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
subplot(3,2,3)
z = 5;
plot(z,z,'g^')
hold on
x = [10^-6, 10^4];
y = [10^-6, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results4BPos(:,5),results7APos(:,5),sz,'black' )
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
title("Organic nitrogen (\muM)")

% Primary Productivity

subplot(3,2,5)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-6, 10^4];
y = [10^-6, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco4B(:,4),eco7A(:,4),sz,'black')
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-3 10^3])
ylim([10^-3 10^3])
xticks([10^-5 10^-3 10^-1 10^1 10^3])
yticks([10^-5 10^-3 10^-1 10^1 10^3])
title("Primary Productivity (\muM N d^{-1})")
ax = gca;
ax.TitleFontSizeMultiplier = 1;

%%% Phosphorus-limited system

% Inorganic phosphorus
subplot(3,2,2)
z = 0.007;
plot(z,z,'g^')
hold on
x = [10^-6, 10^2];
y = [10^-6, 10^2];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,4),results8APos(:,4),sz,'black' )
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
subplot(3,2,4)
z = 0.1;
plot(z,z,'g^')
hold on
x = [10^-7, 10^4];
y = [10^-7, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(results5BPos(:,5),results8APos(:,5),sz,'black' )
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

subplot(3,2,6)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-7, 10^4];
y = [10^-7, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco5B(:,4),eco8A(:,4),sz,'black')
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-7 10^1])
ylim([10^-7 10^1])
xticks([10^-7 10^-5 10^-3 10^-1 10^1])
yticks([10^-7 10^-5 10^-3 10^-1 10^1])
title("Primary Productivity (\muM P d^{-1})")
ax = gca;
ax.TitleFontSizeMultiplier = 1;

set(gcf, 'Position', [800, 800, 600, 1000])
saveas(gcf,'poster2.png')

%%% 3rd graph

%%% Nitrogen-limited system

% DON recycling

figure;
subplot(3,2,1)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-6, 10^4];
y = [10^-6, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco4B(:,1),eco7A(:,1),sz,'black')
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-5 10^3])
ylim([10^-5 10^3])
xticks([10^-5 10^-3 10^-1 10^1 10^3])
yticks([10^-5 10^-3 10^-1 10^1 10^3])
title("DON recycling (\muM d^{-1})")

% Nutrient export

subplot(3,2,3)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-7, 10^4];
y = [10^-7, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco4B(:,2),eco7A(:,2),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-5 10])
ylim([10^-5 10])
xticks([10^-5 10^-3 10^-1 10])
yticks([10^-5 10^-3 10^-1 10])
title("Nutrient Export (\muM N d^{-1}) ")
set(gcf, 'Position', [800, 800, 450, 450])

% Carbon sink

subplot(3,2,5)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-13, 10^4];
y = [10^-13, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco4B(:,3),eco7A(:,3),sz,'black' )
pbaspect([1 1 1])
plot(z,z,'g^')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("without viruses")
ylabel("with viruses")
xlim([10^-3 10^3])
ylim([10^-3 10^3])
xticks([10^-7 10^-5 10^-3 10^-1 10^1 10^3])
yticks([10^-7 10^-5 10^-3 10^-1 10^1 10^3])
title("Carbon sink (\muM C d^{-1})")

%%% Phosphorus-limited system

% DOP recycling

subplot(3,2,2)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-6, 10^4];
y = [10^-6, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco5B(:,1),eco8A(:,1),sz,'black')
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

subplot(3,2,4)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-7, 10^4];
y = [10^-7, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco5B(:,2),eco8A(:,2),sz,'black' )
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

subplot(3,2,6)
z = 10^10;
plot(z,z,'g^','MarkerSize',7,'LineWidth', 2)
hold on
x = [10^-13, 10^4];
y = [10^-13, 10^4];
line(x,y,'Color','red','LineWidth',1.5)
scatter(eco5B(:,3),eco8A(:,3),sz,'black' )
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

set(gcf, 'Position', [800, 800, 600, 1000])
saveas(gcf,'poster3.png')

%%% Regression graph

sz = 10

%8A
V8A = results8APos(:,6) + results8APos(:,7);
logV8A = log10(V8A);
logCS8A = log10(eco8A(:,3));
tblV8A = table(logV8A, logCS8A,'VariableNames',{'V','CS'});
mdlV8A = fitlm(tblV8A,'CS ~ V'); 

x8A = [10^7 10^16];
y8A = 10.^(log10(x8A)*table2array(mdlV8A.Coefficients(2,1)) + table2array(mdlV8A.Coefficients(1,1)));

%7A

V7A = results7APos(:,6) + results7APos(:,7);
logV7A = log10(V7A);
logCS7A = log10(eco7A(:,3));
tblV7A = table(logV7A, logCS7A,'VariableNames',{'V','CS'});
mdlV7A = fitlm(tblV7A,'CS ~ V'); 

x7A = [10^7 10^16];
y7A = 10.^(log10(x7A)*table2array(mdlV7A.Coefficients(2,1)) + table2array(mdlV7A.Coefficients(1,1)));


figure;
scatter(V7A,eco7A(:,3),sz)
hold on
scatter(V8A,eco8A(:,3),sz,'red')
line(x8A,y8A,'LineWidth',1.5,'Color','red');
line(x7A,y7A,'LineWidth',1.5);
pbaspect([1 1 1])
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("Total virus abundance (L^{-1})")
ylabel("Carbon sink (\muM C d^{-1})")
xlim([10^9 10^16])
ylim([10^-3 10^4])
xticks([10^9 10^11 10^13 10^15])
yticks([10^-3 10^-1 10^1 10^3])
hold off
set(gcf, 'Position', [800, 800, 300, 300])
saveas(gcf,'poster4.png')
