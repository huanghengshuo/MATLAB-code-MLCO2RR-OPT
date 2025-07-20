% Pearson & Spearman ori
% Coorrelation coefficients


clc
clear all
close all
rng('default') % For reproducibility
%   Seperated analysis for ALL
%   Seperated analysis for ALL
%   Seperated analysis for ALL
load('c.mat')
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

phi_ref = [-1.15 -1.03 -0.903 -0.764 -0.624]; % vs RHE
% phi_ref = phi_ref - 0.05916*pH; % vs NHE
jCO = [70.1 47.3 29.3 13.6 3.8 ]; % mA/cm^2
jH2 = [2.0 1.5 1.6 1.0 0.5]; % mA/cm^2
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


load("ExpAngew2024Data004.mat")
Angewx = -0.6:-0.1:-0.9;
Angewy = ExpAngew2024Data004(:,2)';

load("ExpJMCA2019Data005.mat")
JMCAx = -1.2:0.1:-0.4;
JMCAy = ExpJMCA2019Data005(:,2)';

JMCAx(end-1:end) = [];
JMCAy(end-1:end) = [];

%%
fig = figure;
set(gca,'FontSize',18)
% fig.Color = 'none';
hold on
box on
picture1 = plot(MFE(:,1)+0.05916*pH,MFE(:,2),'k-','LineWidth',2.5);

picturea = scatter(phi_ref,FE_ref,'SizeData',100,...
    'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','m',...
    'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.5);
pict1 = plot(phi_ref,FE_ref,':m','LineWidth',1.5);




pictureAngew = scatter(Angewx,Angewy,'SizeData',100,...
    'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','g',...
    'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.5);
pict2 = plot(Angewx,Angewy,':g','LineWidth',1.5);

pictureJMCA = scatter(JMCAx,JMCAy,'SizeData',100,...
    'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','r',...
    'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.5);
pict3 = plot(JMCAx,JMCAy,':r','LineWidth',1.5);

ylim([0.6 1])
xlim([-1.58 -0.4])
xlabel("Potential (V \it vs. \rm RHE)","FontSize",18,"FontName","Arial");
ylabel("FE","FontSize",18,"FontName","Arial");


legend([picture1 picturea pictureAngew pictureJMCA],[...
    "Simulation",...
    "refence a",...
    "refence b",...
    "refence c"],"FontSize",15,"FontName","Arial",...
    "Location","southeast");

ax = gca;
set(ax,'Color','none');
% print('0Fig00 Test good comsol model 01','-djpeg','-r1200') 

%%
fig = figure;
set(gca,'FontSize',18)
% fig.Color = 'none';
hold on
box on
picture1 = plot(MFE(:,1)+0.05916*pH,MFE(:,2),'k-','LineWidth',2.5);


phi_refa(end-4:end) = [];
FE_refa(end-4:end) = [];

picturea = scatter(phi_refa,FE_refa,'SizeData',100,...
    'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','m',...
    'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.5);

pict1 = plot(phi_refa,FE_refa,':m','LineWidth',1.5);




pictureAngew = scatter(Angewx,Angewy,'SizeData',100,...
    'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','g',...
    'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.5);
pict2 = plot(Angewx,Angewy,':g','LineWidth',1.5);

pictureJMCA = scatter(JMCAx,JMCAy,'SizeData',100,...
    'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','r',...
    'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.5);
pict3 = plot(JMCAx,JMCAy,':r','LineWidth',1.5);

ylim([0.6 1])
xlim([-1.58 -0.4])
xlabel("$$ \mathrm{Potential}\ \left( \mathrm{V}\ vs. \ \mathrm{RHE} \right)$$","FontSize",18,"Interpreter","latex");
ylabel("$$ FE_\mathrm{CO} $$","FontSize",18,"Interpreter","latex");
legend([picture1 picturea pictureAngew pictureJMCA],[...
    "$$ \mathrm{Simulation} $$",...
    "$$ \mathrm{refence\ a} $$",...
    "$$ \mathrm{refence\ b} $$",...
    "$$ \mathrm{refence\ c} $$"],"FontSize",15,"Interpreter","latex",...
    "Location","southeast");

ax = gca;
% set(ax,'Color','none');
% print('0Fig00 Test good comsol model 02','-djpeg','-r1200') 

%% FEexp VS FEsim
% xq = [phi_refa',Angewx,JMCAx];
% FEsim = spline(MFE(:,1)+0.05916*pH,MFE(:,2),xq);
% FEexp = [FE_refa',Angewy,JMCAy];

xq = [phi_refa'];
FEsim = spline(MFE(:,1)+0.05916*pH,MFE(:,2),xq);
FEexp = [FE_refa'];



SSE = sum((FEsim -       FEexp).^2);
SST = sum((FEexp - mean(FEexp)).^2);
rsquare = 1 - SSE./SST;
figure3 = figure;
set(gca,'FontSize',18)
hold on
% grid minor
box on

plot([0.85 1],[0.85 1],'k--','LineWidth',1.5)
scatter(FEexp,FEsim,'SizeData',100,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','c',...
    'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.5)

xlim([0.85 1])
ylim([0.85 1])
xlabel("Experiment FE","FontSize",18,"FontName","Arial");
% xlabel("$ \textsf{Experimental\ FE} $","FontSize",20,"Interpreter","latex");
ylabel("Predicted FE","FontSize",18,"FontName","Arial");
% print('0Fig00 Test good comsol model FEsim VS FEexp 01','-djpeg','-r1200') 


%% Input original data
M0 = readmatrix("jCOER multi test.xlsx"); % GDL CL  phi
M1 = readmatrix("jHER multi test.xlsx");
M2 = readmatrix("jTotal multi test.xlsx");
MCOER = rmmissing(M0);
MHER  = rmmissing(M1);
MTotal= rmmissing(M2);

MFE = MCOER;
% MFE(:,end) = MCOER(:,end)./(MCOER(:,end)+MTotal(:,end));
MFE(:,end) = MCOER(:,end)./MTotal(:,end);


% xCL =  unique(MFE(:,2))';
% 
% index = find(MFE(:,2)== xCL(6));
% MFE = MFE(index,:);
% MCOER = MCOER(index,:);
% 
index = find(MFE(:,1)<=MFE(:,2));
MFE(index,:)   = [];
MCOER(index,:) = [];


xGDL = unique(MFE(:,1))';
xCL =  unique(MFE(:,2))';
phix = unique(MFE(:,3))';



x = MFE(:,2);
y = MFE(:,3);
v = MFE(:,end);

[xq,yq] = meshgrid(linspace(min(xCL),max(xCL),500), linspace(min(phix),max(phix),500));
vq = griddata(x,y,v,xq,yq);

fig = figure;
fig.Position = [742.6,552.2,628,420];

box on
% splot = contourf(xq,yq,vq,[0.70:0.05:0.90,0.93:0.01:0.95,0.955],'ShowText','on');
splot = contourf(xq,yq,vq,[0.70:0.05:0.90,0.93:0.01:0.95,0.955]);
set(gca,'FontSize',18,'Position',[0.135031847133758,0.157619050888788,0.672611464968153,0.767380949111212]);
% xlim([0.45 0.7])

view([0 90])

c = colormap(hot);
c = rot90(c,2);
c = flipud(c);
colormap(flipud(c1));
c = colorbar;
c.Position = [0.851726071488318,0.157619050888788,0.033970276008493,0.609523809523811]
% c.Label.String = 'FE';
xlabel("Porosity of CL","FontSize",18,"FontName","Arial");

% xlabel("$$ \varepsilon_\mathrm{CL} $$","Interpreter","latex","FontSize",20)
ylabel("Potential (V \it vs. \rm NHE)","FontSize",18,"FontName","Arial");
% title("$$ \mathrm{residual\ FE} $$","Interpreter","latex","FontSize",20)
% zlabel("$$ residual\ FE $$","Interpreter","latex")
annotation(fig,'textbox',...
    [0.867242038216561 0.822857142857145 0.0792547770700638 0.0761904761904765],...
    'String',{'FE'},...
    'LineStyle','none',...
    'FontSize',18,...
    'FontName','Arial',...
    'FitBoxToText','off');
% print('0Fig05 FE vs CL with phi','-djpeg','-r1200') 


toc
