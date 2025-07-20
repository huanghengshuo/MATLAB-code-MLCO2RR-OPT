clc
clear all
close all
rng('default') % For reproducibility
%   Seperated analysis for ALL
%   Seperated analysis for ALL
%   Seperated analysis for ALL
tic
%% Input original data
% M0 = readmatrix("jCOER multi test.xlsx"); % GDL CL  phi
% M1 = readmatrix("jHER multi test.xlsx");
% M2 = readmatrix("jTotal multi test.xlsx");
M0 = readmatrix("jCOER test.xlsx"); % GDL CL  phi
M1 = readmatrix("jHER test.xlsx");
M2 = readmatrix("jTotal test.xlsx");
MCOER = rmmissing(M0);
MHER  = rmmissing(M1);
MTotal= rmmissing(M2);

MFE = MCOER;
% MFE(:,end) = MCOER(:,end)./(MCOER(:,end)+MTotal(:,end));
MFE(:,end) = MCOER(:,end)./MTotal(:,end);

phix = unique(MFE(:,1))';


%%
pH = 7.23;
% phi_ref = [-1.138 -1.02 -0.901 -0.756 -0.615 -0.479]; % vs RHE
% % phi_ref = phi_ref - 0.05916*pH; % vs NHE
% jCO = [69.305 45.997 28.481 12.989 4.546 0.162]; % mA/cm^2
% jH2 = [1.902 1.455 1.565 0.948 0.385 0.159]; % mA/cm^2

phi_ref = [-1.15 -1.03 -0.903 -0.764 -0.624 -0.474]; % vs RHE
% phi_ref = phi_ref - 0.05916*pH; % vs NHE
jCO = [70.1 47.3 29.3 13.6 3.8 0.01]; % mA/cm^2
jH2 = [2.0 1.5 1.6 1.0 0.5 0.3]; % mA/cm^2
iR = (jCO+jH2)*10*0.75e-4/0.5; % iR correction
% phi_ref = phi_ref + iR;
FE_ref = jCO./(jCO+jH2);
% phi = -0.2:-0.02:-2.4; % vs NHE
phix = phix + 0.05916*pH;     % vs RHE

%%
jCOa = readmatrix("0 data CO2.xlsx");
jH2a = readmatrix("0 data H2.xlsx");
phi_refa = jCOa(:,1);
jCOa = spline(jCOa(:,1),jCOa(:,2),phi_refa);
jH2a = spline(jH2a(:,1),jH2a(:,2),phi_refa);
FE_refa = jCOa./(jCOa+jH2a);
%% AM 2022--Synergistic Cr2O3@Ag Heterostructure Enhanced
load("ExpAM2022Data001.mat")
AMx = -1.4:0.1:-0.4;
AMy = ExpAM2022Data001(:,2)';
AMy = AMy./100;

load("ExpAMI2022Data002.mat")
AMIx = -0.5:-0.1:-1.1;
AMIy1 = 1-ExpAMI2022Data002(1:7,2)';
AMIy2 = 1-ExpAMI2022Data002(8:end,2)';

load("ExpNANO2022Data001.mat")
NANOx = -0.8:-0.2:-1.6;
NANOy = ExpNANO2022Data001(:,2)';

load("ExpAngew2024Data004.mat")
Angewx = -0.6:-0.1:-0.9;
Angewy = ExpAngew2024Data004(:,2)';

load("ExpJMCA2019Data005.mat")
JMCAx = -1.2:0.1:-0.4;
JMCAy = ExpJMCA2019Data005(:,2)';

fig = figure;
% fig.Color = 'none';
hold on

picture1 = plot(MFE(:,1)+0.05916*pH,MFE(:,2),'k-','LineWidth',2.5);

picturea = scatter(phi_ref,FE_ref,'SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','m',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);

pictureNANO = scatter(NANOx,NANOy,'SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);

pictureAngew = scatter(Angewx,Angewy,'SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','g',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);

pictureJMCA = scatter(JMCAx,JMCAy,'SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','r',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
% picture0 = scatter(phi_ref,FE_ref,'SizeData',50,...
%     'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','r',...
%     'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);

ylim([0 1])
xlabel("$$ \phi \ \mathrm{vs\ RHE}$$","FontSize",15,"Interpreter","latex");
ylabel("$$ FE_\mathrm{CO} $$","FontSize",15,"Interpreter","latex");
legend([picture1 picturea pictureNANO pictureAngew pictureJMCA],[...
    "$$ \mathrm{Theoretical\ data } $$",...
    "$$ \mathrm{Experimental\ data\ Ag} $$",...
    "$$ \mathrm{Experimental\ data\ Ag900D} $$",...
    "$$ \mathrm{Experimental\ data\ AgMOCNO_3} $$",...
    "$$ \mathrm{Experimental\ data\ AgCLNP} $$"],"FontSize",10,"Interpreter","latex",...
    "Location","southeast");

ax = gca;
set(ax,'Color','none');
print('0Fig01 Test good comsol model','-djpeg','-r1200') 
% exportgraphics(fig,'0Fig01 Test good comsol model export.jpg')
% saveas(fig,'0Fig01 Test good comsol model saveas.jpg');
figure
hold on
scatter(phi_ref,jCO,'*');
plot(MFE(:,1)+0.05916*pH,-MCOER(:,2),'k-');


%% FEexp VS FEsim

load("ExpNANO2022Data001.mat")
NANOx = -0.8:-0.2:-1.6;
NANOy = ExpNANO2022Data001(:,2)';

load("ExpAngew2024Data004.mat")
Angewx = -0.6:-0.1:-0.9;
Angewy = ExpAngew2024Data004(:,2)';

load("ExpJMCA2019Data005.mat")
JMCAx = -1.2:0.1:-0.4;
JMCAy = ExpJMCA2019Data005(:,2)';


xq = [phi_ref,NANOx,Angewx,JMCAx];
FEsim = spline(MFE(:,1)+0.05916*pH,MFE(:,2),xq);
FEexp = [FE_ref,NANOy,Angewy,JMCAy];

figure3 = figure;
hold on
grid minor
box on
scatter(FEexp,FEsim,'SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','c',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7)
plot([0.75 1],[0.75 1],'k-','LineWidth',2)
xlim([0.75 1])
ylim([0.75 1])
xlabel("$$ \mathrm{Experimental\ FE} $$","FontSize",13,"Interpreter","latex");
ylabel("$$ \mathrm{Predicted\ FE} $$","FontSize",13,"Interpreter","latex");



figure3 = figure;
hold on
grid minor
box on
scatter(FEexp,FEsim,'SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','c',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7)
plot([min(FEsim) max(FEsim)],[min(FEsim) max(FEsim)],'k-','LineWidth',2)

xlim([min(FEsim) max(FEsim)])
ylim([min(FEsim) max(FEsim)])

% title("$$ \mathrm{Decision\ Tree\ Model} $$","FontSize",13,"Interpreter","latex");
xlabel("$$ \mathrm{Experiment\ FE} $$","FontSize",13,"Interpreter","latex");
ylabel("$$ \mathrm{Simulation\ FE} $$","FontSize",13,"Interpreter","latex");



