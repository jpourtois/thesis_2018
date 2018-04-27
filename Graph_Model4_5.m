% Make alternate graph comparing P and N-limited conditions for Model 4 and
% Model 5. 
% FIRST RUN 'Comp4A_4B.m' and 'Comp5A_5B.m'

%%% First Graph %%%
figure;

%%% Nitrogen Limited

% Heterotrophic bacteria

subplot(2,3,1)
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

subplot(2,3,2)
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

subplot(2,3,3)
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

%%% Phosphorus Limited

% Heterotrophic bacteria

subplot(2,3,4)
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

subplot(2,3,5)
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

subplot(2,3,6)
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

set(gcf, 'Position', [800, 800, 1000, 650])
saveas(gcf,'4A_5A_1.png')

%%% Second figure %%%

figure;

%%% Nitrogen limited

% Inorganic nitrogen
subplot(2,3,1)
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
subplot(2,3,2)
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

subplot(2,3,3)
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

%%% Phosphorus-limited

% Inorganic phosphorus
subplot(2,3,4)
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
subplot(2,3,5)
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

subplot(2,3,6)
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

set(gcf, 'Position', [800, 800, 1000, 650])
saveas(gcf,'4A_5A_2.png')

%%% Third Graph %%%
figure;
%%% Nitrogen-limited

% DOP recycling

subplot(2,3,1)
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

subplot(2,3,2)
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

% Carbon sink

subplot(2,3,3)
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

%%% Phosphorus-limited

% DOP recycling

subplot(2,3,4)
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

subplot(2,3,5)
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

subplot(2,3,6)
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

set(gcf, 'Position', [800, 800, 1000, 650])
saveas(gcf,'4A_5A_3.png')
