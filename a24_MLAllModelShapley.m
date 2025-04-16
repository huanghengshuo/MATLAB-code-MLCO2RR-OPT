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

% phix = unique(MFE(:,3))';
% for i = 2:2:length(phix)-1
%     
%     index = find(abs(MFE(:,3)-phix(i))<0.001);
%     MFE(index,:) = [];
%     MCOER(index,:) = [];
% end

phix = unique(MFE(:,3))';

%% plot swarmchart
figure
hold on
yyaxis left
set(gca,'YColor','k')
X = MFE(:,1:end-1);
Y = MFE(:,end);
ystd = zeros(1,length(phix));
cmap = colormap(autumn(length(phix)));
for i = 1:5:length(phix)
    
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    
    X1 = X(index,1:2);
    Y1 = Y(index,:);
    Y1 = Y1 - mean(Y1);
    
    x = ones(size(Y1,1),1)*phix(i)*1;
    % x = ones(size(Y1,1),1)*i;
    s = swarmchart(x,Y1,5,cmap(i,:),'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    s.XJitterWidth = 0.05;
end

%% Prepare the data for Machine Learning
YmeanInd = zeros(size(phix));
Yminusmean = zeros(size(Y));

for i = 1:1:length(phix)
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    
    X1 = X(index,1:2);
    Y1 = Y(index,:);
    ystd(1,i) = std(Y1);
    YmeanInd(1,i) = mean(Y1);
    Yminusmean(index,:) = Y1 - mean(Y1);
end
yyaxis right
set(gca,'YColor','k')
plot(phix,ystd,'k-','LineWidth',2)
ylim([-max(ystd)*1.0 max(ystd)*1.05])


figure
hold on

plot(X(:,3),Y,'k.')
plot(phix,YmeanInd,'c-','LineWidth',2)


figure
hold on
plot(X(:,3),Yminusmean,'k.')
plot(phix,YmeanInd-YmeanInd,'c-','LineWidth',2)



%% Determine the Performance of ALL Machine Learning
KernelFunctionName = {'ardsquaredexponential';...
    'ardexponential';...
    'ardmatern32';...
    'ardmatern52';...
    'ardrationalquadratic'};
% 1 done
% 2 done
% 3 done
% 4
% 5
for j = 1:1
    ModelCVName= "RGPCV"+string(KernelFunctionName{j})+".mat";
    ModelName= "RGP"+string(KernelFunctionName{j})+".mat";
    load(ModelName)
    load(ModelCVName)
ShapleyValue = zeros(size(X));
warning('off')
% for i = 1:10
for i = 1:size(X,1)

explainer1 = shapley(GaussModel,'QueryPoint',X(i,:),'UseParallel','on');

% figure
% plot(explainer1)
ShapleyValue(i,:) = table2array(explainer1.ShapleyValues(:,2))';
end
varnames = ["epi_GDL","epi_CL","phi"];
writematrix(varnames,"Shapley value All Gauss.xlsx","Sheet",string(KernelFunctionName{j}),"Range","B1")

writematrix(ShapleyValue,"Shapley value All Gauss.xlsx","Sheet",string(KernelFunctionName{j}),"Range","B2")
end

%% Tree Model Shalpey
CVname = "Tree"+"All.mat";
ModelCVName= "TreeCV"+".mat";
ModelName= "Tree"+".mat";
load(ModelName)
load(ModelCVName)
ShapleyValue = zeros(size(X));
warning('off')
% for i = 1:10
for i = 1:size(X,1)

explainer1 = shapley(TreeModel,'QueryPoint',X(i,:),'UseParallel','on');

% figure
% plot(explainer1)
ShapleyValue(i,:) = table2array(explainer1.ShapleyValues(:,2))';
end
varnames = ["epi_GDL","epi_CL","phi"];
writematrix(varnames,"Shapley value All.xlsx","Sheet","Tree","Range","B1")

writematrix(ShapleyValue,"Shapley value All.xlsx","Sheet","Tree","Range","B2")




toc