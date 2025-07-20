% Pearson & Spearman ori-mean
% Coorrelation coefficients

clc
clear all
close all
rng('default') % For reproducibility
%   Seperated analysis for ALL
%   Seperated analysis for ALL
%   Seperated analysis for ALL
tic
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

index = find(MFE(:,1)<=MFE(:,2));
MFE(index,:)   = [];
MCOER(index,:) = [];

index = find(MFE(:,3)<=-2.001);
MFE(index,:)   = [];
MCOER(index,:) = [];
% phix = -0.600:-0.005:-2.400;

phix = unique(MFE(:,3))';
% To remove some phi value 
% for i = 2:2:length(phix)-1
%     
%     index = find(abs(MFE(:,3)-phix(i))<0.001);
%     MFE(index,:) = [];
%     MCOER(index,:) = [];
% end

phix = unique(MFE(:,3))';
X = MFE(:,1:end-1);
Y = MFE(:,end);


%% Prepare the data for Machine Learning
YmeanInd = zeros(size(phix));
Yminusmean = zeros(size(Y));

for i = 1:1:length(phix)
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    Y1 = Y(index,:);
    YmeanInd(1,i) = mean(Y1);
    Yminusmean(index,:) = Y1 - mean(Y1);
end

%%
% epiGDL = MFE(:,1);
% epiCL = MFE(:,2);
% phi = MFE(:,3);
MFE(:,4) = Yminusmean;
phix = unique(MFE(:,3))';
rhoP = partialcorr(MFE,'Type','Pearson');

rhoS = partialcorr(MFE,'Type','Spearman');

xvalues = {'GDL','CL','phi','FECO'};
yvalues = xvalues;
f1 = figure;
h1 = heatmap(f1,xvalues,yvalues,rhoP,'GridVisible','off','FontName',...
    'Times New Roman','FontSize',12,'CellLabelFormat','%0.4f','ColorLimits',[-0.3 1],'Interpreter','latex');
h1.XDisplayLabels = {'$$\varepsilon_\mathrm{GDL} $$','$$\varepsilon_\mathrm{CL} $$','$$ \varphi $$','$$\widetilde{\mathrm{FE}} $$'};
h1.YDisplayLabels = h1.XDisplayLabels;

title('Pearson partial correlation');
print('0Fig01bb partialcorrPearson residual FE} ','-djpeg','-r1200') 

f2 = figure;
h2 = heatmap(f2,xvalues,yvalues,rhoS,'GridVisible','off','FontName',...
    'Times New Roman','FontSize',12,'CellLabelFormat','%0.4f','ColorLimits',[-0.3 1],'Interpreter','latex');
h2.XDisplayLabels = {'$$\varepsilon_\mathrm{GDL} $$','$$\varepsilon_\mathrm{CL} $$','$$ \varphi $$','$$\widetilde{\mathrm{FE}} $$'};
h2.YDisplayLabels = h2.XDisplayLabels;

title('Spearman partial correlation');
print('0Fig01bb partialcorrSpearman residual FE','-djpeg','-r1200') 
 

%%
rhoP = corr(MFE,'Type','Pearson');
rhoS = corr(MFE,'Type','Spearman');
f3 = figure;
h3 = heatmap(f3,xvalues,yvalues,rhoP,'GridVisible','off','FontName',...
    'Times New Roman','FontSize',12,'CellLabelFormat','%.4f','ColorLimits',[-0.3 1],'Interpreter','latex');
h3.XDisplayLabels = {'$$\varepsilon_\mathrm{GDL} $$','$$\varepsilon_\mathrm{CL} $$','$$ \varphi $$','$$\widetilde{\mathrm{FE}}  $$'};
h3.YDisplayLabels = h3.XDisplayLabels;

title('Pearson correlation');
print('0Fig01bb corrPearson residual FE','-djpeg','-r1200') 


f4 = figure;
h4 = heatmap(f4,xvalues,yvalues,rhoS,'GridVisible','off','FontName',...
    'Times New Roman','FontSize',12,'CellLabelFormat','%0.4f','ColorLimits',[-0.3 1],'Interpreter','latex');
h4.XDisplayLabels = {'$$\varepsilon_\mathrm{GDL} $$','$$\varepsilon_\mathrm{CL} $$','$$ \varphi $$','$$\widetilde{\mathrm{FE}} $$'};
h4.YDisplayLabels = h4.XDisplayLabels;

title('Spearman partial correlation');
print('0Fig01bb corrSpearman residual FE','-djpeg','-r1200') 



toc
