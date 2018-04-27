% Compare Diff7A_4B and Diff8A_5B. 
% FIRST RUN 'Comp7A_4B.m' and 'Comp8A_5B.m'

% This does:
% 1. Run Wilcoxon rank sum test --> Compare effect of viruses between
% N-limited and P-limited system
% 2. Run ranksum for viruses (includes graph) --> Compare virus abundance in 7A and 8A
% 3. Linear Regression Viruses and bacteria (includes graph)
% 4. Linear Regression Viruses and carbon sink
% 5. Linear Regression Zooplankton and carbon sink
% 6. Linear Regression Viruses to zooplankton ratio and carbon sink


%%% 1. Run Wilcoxon rank sum test
for i = 1:length(diffLogResults7A_4B(1,:))
    
    pCompResults(i) = ranksum(diffLogResults7A_4B(:,i),diffLogResults8A_5B(:,i));
    
end

for j = 1:length(diffLogEco7A_4B(1,:))
    
    pCompEco(j) = ranksum(diffLogEco7A_4B(:,j), diffLogEco8A_5B(:,j));
end

%%% 2. Run ranksum for viruses

logVh7A = log10(results7APos(:,7));
logVc7A = log10(results7APos(:,6));
logVtot7A = log10(results7APos(:,7)+results7APos(:,6));

logVh8A = log10(results8APos(:,6));
logVc8A = log10(results8APos(:,7));
logVtot8A = log10(results8APos(:,6)+results8APos(:,7));

pRankVh = ranksum(logVh7A, logVh8A);
pRankVc = ranksum(logVc7A, logVc8A);
pRankVtot = ranksum(logVtot7A,logVtot8A);

% Graph 

figure;

% N-limited

vector1_5 = transpose(repelem(1.5,171));
vector3 = transpose(repelem(3,171));
vector4_5 = transpose(repelem(4.5,171));

n_vectorNum = vertcat(vector1_5, vector3, vector4_5);

n_vc = results7APos(:,6);
n_vh = results7APos(:,7);
n_vtot = n_vc + n_vh;

med_n_vc = median(n_vc);
med_n_vh = median(n_vh);
med_n_vtot = median(n_vtot);

n_vectorVal = vertcat(n_vc, n_vh, n_vtot);

scatter(n_vectorVal,n_vectorNum);
%scatter([med_n_vc med_n_vh med_n_vtot], [1.5 3 4.5]);
hold on

% P-limited

vector1 = transpose(repelem(1,125));
vector2_5 = transpose(repelem(2.5,125));
vector4 = transpose(repelem(4,125));

p_vectorNum = vertcat(vector1, vector2_5, vector4);

p_vc = results8APos(:,7);
p_vh = results8APos(:,6);
p_vtot = p_vc + p_vh;

med_p_vc = median(p_vc);
med_p_vh = median(p_vh);
med_p_vtot = median(p_vtot);

p_vectorVal = vertcat(p_vc, p_vh, p_vtot);

scatter(p_vectorVal,p_vectorNum,[],[0.9 0.5 0]);
%scatter([med_p_vc med_p_vh med_p_vtot], [1 2.5 4]);
legend('N-limited', 'P-limited')
xlabel('Abundance (L^{-1})')
ylim([0 5.5])
yticks([1.25 2.75 4.25])
yticklabels({'V_{C}', 'V_{H}', 'V_{tot}'})
set(gca,'xscale','log')
set(gca, 'FontSize', 14)
LEG.FontSize = 16;
set(gcf, 'Position', [800, 800, 500, 400])
saveas(gcf,'virus.png')

%%% 3. Linear Regression Viruses and bacteria

x = [10^8 10^16];

% Vh 

tbl_Vh_H_7A = table(log10(results7APos(:,7)),log10(results7APos(:,1)),'VariableNames',{'Vh','H'});
mdl_Vh_H_7A = fitlm(tbl_Vh_H_7A, 'H ~ Vh');

tbl_Vh_H_8A = table(log10(results8APos(:,6)),log10(results8APos(:,1)),'VariableNames',{'Vh','H'});
mdl_Vh_H_8A = fitlm(tbl_Vh_H_8A, 'H ~ Vh');

y_Vh_H_8A = 10.^(log10(x)*table2array(mdl_Vh_H_8A.Coefficients(2,1)) + table2array(mdl_Vh_H_8A.Coefficients(1,1)));

% Vc

tbl_Vc_C_7A = table(log10(results7APos(:,6)),log10(results7APos(:,2)),'VariableNames',{'Vc','C'});
mdl_Vc_C_7A = fitlm(tbl_Vc_C_7A, 'C ~ Vc');

y_Vc_C_7A = 10.^(log10(x)*table2array(mdl_Vc_C_7A.Coefficients(2,1)) + table2array(mdl_Vc_C_7A.Coefficients(1,1)));

tbl_Vc_C_8A = table(log10(results8APos(:,7)),log10(results8APos(:,2)),'VariableNames',{'Vc','C'});
mdl_Vc_C_8A = fitlm(tbl_Vc_C_8A, 'C ~ Vc');

y_Vc_C_8A = 10.^(log10(x)*table2array(mdl_Vc_C_8A.Coefficients(2,1)) + table2array(mdl_Vc_C_8A.Coefficients(1,1)));

% Vtot

tbl_Vtot_Btot_7A = table(log10(results7APos(:,6)+results7APos(:,7)),log10(results7APos(:,2)+results7APos(:,1)),'VariableNames',{'Vtot','Btot'});
mdl_Vtot_Btot_7A = fitlm(tbl_Vtot_Btot_7A, 'Btot ~ Vtot');

y_Vtot_Btot_7A = 10.^(log10(x)*table2array(mdl_Vtot_Btot_7A.Coefficients(2,1)) + table2array(mdl_Vtot_Btot_7A.Coefficients(1,1)));

tbl_Vtot_Btot_8A = table(log10(results8APos(:,7)+results8APos(:,6)),log10(results8APos(:,2)+results8APos(:,1)),'VariableNames',{'Vtot','Btot'});
mdl_Vtot_Btot_8A = fitlm(tbl_Vtot_Btot_8A, 'Btot ~ Vtot');

y_Vtot_Btot_8A = 10.^(log10(x)*table2array(mdl_Vtot_Btot_8A.Coefficients(2,1)) + table2array(mdl_Vtot_Btot_8A.Coefficients(1,1)));

% Graph
figure;
subplot(1,3,1)
scatter(results7APos(:,7), results7APos(:,1))
hold on
scatter(results8APos(:,6), results8APos(:,1),[],[0.9 0.5 0])
line(x,y_Vh_H_8A,'LineWidth',1.5,'Color', [0.9 0.5 0]);
set(gca,'xscale','log')
set(gca,'yscale','log')
pbaspect([1 1 1])
xlabel("Heterotrophic bacteria viruses (L^{-1})",'FontSize',14)
ylabel("Heterotrophic bacteria (L^{-1})",'FontSize',14)
legend({'N-limited', 'P-limited'},'FontSize',12)

subplot(1,3,2)
scatter(results7APos(:,6), results7APos(:,2))
hold on
scatter(results8APos(:,7), results8APos(:,2),[],[0.9 0.5 0])
line(x,y_Vc_C_7A,'LineWidth',1.5,'Color', [0 0.5 0.7]);
line(x,y_Vc_C_8A,'LineWidth',1.5,'Color', [0.9 0.5 0]);
set(gca,'xscale','log')
set(gca,'yscale','log')
pbaspect([1 1 1])
xlabel("Cyanobacteria viruses (L^{-1})",'FontSize',14)
ylabel("Cyanobacteria (L^{-1})",'FontSize',14)
legend({'N-limited', 'P-limited'},'FontSize',12)

subplot(1,3,3)
scatter(results7APos(:,6)+results7APos(:,7), results7APos(:,2)+results7APos(:,1))
hold on
scatter(results8APos(:,7) + results8APos(:,6), results8APos(:,2)+ results8APos(:,1),[],[0.9 0.5 0])
line(x,y_Vtot_Btot_7A,'LineWidth',1.5,'Color', [0 0.5 0.7]);
line(x,y_Vtot_Btot_8A,'LineWidth',1.5,'Color', [0.9 0.5 0]);
set(gca,'xscale','log')
set(gca,'yscale','log')
pbaspect([1 1 1])
xlabel("Total viruses (L^{-1})",'FontSize',14)
ylabel("Total bacteria (L^{-1})",'FontSize',14)
legend({'N-limited', 'P-limited'},'FontSize',12)

set(gcf, 'Position', [800, 800, 1200, 400])
saveas(gcf,'virus_bacteria.png')

% Stats 

ranksum(log10(results7APos(:,7)./results7APos(:,1)), log10(results8APos(:,6)./results8APos(:,1)))
ranksum(log10(results7APos(:,2)./results7APos(:,6)), log10(results8APos(:,2)./results8APos(:,7)))
ranksum(log10((results7APos(:,6)+results7APos(:,7))./(results7APos(:,2)+results7APos(:,1))), log10((results8APos(:,7) + results8APos(:,6))./(results8APos(:,2)+ results8APos(:,1))))

%%% 4. Linear Regression Viruses and carbon sink

sz = 7

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
subplot(1,3,1)
scatter(V7A,eco7A(:,3),sz)
hold on
scatter(V8A,eco8A(:,3),sz,[0.9 0.5 0])
line(x8A,y8A,'LineWidth',1.5,'Color',[0.9 0.5 0]);
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
legend('N-limited', 'P-limited')
hold off

%%% 5. Linear Regression Zooplankton and carbon sink

%8A

Z8A = results8APos(:,3);
logZ8A = log10(Z8A);
logCS8A = log10(eco8A(:,3));
tblZ8A = table(logZ8A, logCS8A,'VariableNames',{'Z','CS'});
mdlZ8A = fitlm(tblZ8A,'CS ~ Z'); 

x8A = [10^0 10^6];
y8A = 10.^(log10(x8A)*table2array(mdlZ8A.Coefficients(2,1)) + table2array(mdlZ8A.Coefficients(1,1)));

%7A

Z7A = results7APos(:,3);
logZ7A = log10(Z7A);
logCS7A = log10(eco7A(:,3));
tblZ7A = table(logZ7A, logCS7A,'VariableNames',{'Z','CS'});
mdlZ7A = fitlm(tblZ7A,'CS ~ Z'); 

x7A = [10^0 10^6];
y7A = 10.^(log10(x7A)*table2array(mdlZ7A.Coefficients(2,1)) + table2array(mdlZ7A.Coefficients(1,1)));


subplot(1,3,2)
scatter(Z7A,eco7A(:,3),sz)
hold on
scatter(Z8A,eco8A(:,3),sz,[0.9 0.5 0])
% No relationship so no line
%line(x8A,y8A,'LineWidth',1.5,'Color','red');
%line(x7A,y7A,'LineWidth',1.5);
pbaspect([1 1 1])
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("Zooplankton abundance (L^{-1})")
ylabel("Carbon sink (\muM C d^{-1})")
legend('N-limited', 'P-limited')
hold off

%%% 6. Linear Regression Viruses to zooplankton ratio and carbon sink

%8A

VZratio8A = (results8APos(:,6) + results8APos(:,7))./results8APos(:,3);
logVZratio8A = log10(VZratio8A);
logCS8A = log10(eco8A(:,3));
tbl8A = table(logVZratio8A, logCS8A,'VariableNames',{'VZratio','CS'});
mdl8A = fitlm(tbl8A,'CS ~ VZratio'); 

x8A = [10^5 10^14];
y8A = 10.^(log10(x8A)*table2array(mdl8A.Coefficients(2,1)) + table2array(mdl8A.Coefficients(1,1)));

%7A

VZratio7A = (results7APos(:,6) + results7APos(:,7))./results7APos(:,3);

logVZratio7A = log10(VZratio7A);
logCS7A = log10(eco7A(:,3));
tbl7A = table(logVZratio7A, logCS7A,'VariableNames',{'VZratio','CS'});
mdl7A = fitlm(tbl7A,'CS ~ VZratio'); 

x7A = [10^5 10^14];
y7A = 10.^(log10(x7A)*table2array(mdl7A.Coefficients(2,1)) + table2array(mdl7A.Coefficients(1,1)));


subplot(1,3,3)
scatter(VZratio7A,eco7A(:,3),sz)
hold on
scatter(VZratio8A,eco8A(:,3),sz,[0.9 0.5 0])
line(x8A,y8A,'LineWidth',1.5,'Color',[0.9 0.5 0]);
line(x7A,y7A,'LineWidth',1.5);
pbaspect([1 1 1])
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("Virus to zooplankton ratio")
ylabel("Carbon sink (\muM C d^{-1})")
xlim([10^5 10^14])
ylim([10^-3 10^6])
xticks([10^5 10^7 10^9 10^11 10^13])
yticks([10^-3 10^-1 10^1 10^3 10^5])
legend('N-limited', 'P-limited')
hold off

set(gcf, 'Position', [800, 800, 800, 300])

saveas(gcf,'regressions.png')


