clc
clear all
close all
rng("default")
load("c.mat")
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
%%
KernelFunctionName = {'ardsquaredexponential';...
    'ardexponential';...
    'ardmatern32';...
    'ardmatern52';...
    'ardrationalquadratic'};

%
% for j = 1:2
% ShapleyValue = readmatrix("Shapley value All.xlsx","Sheet",string(KernelFunctionName{j}),"Range","B2");
% %
% %% Determine each feature
%
%
% x = ones(size(ShapleyValue,1),1)*1;
% yShapleyValue = ShapleyValue(:,1);
% c = 1:size(ShapleyValue,1);
% % c(1:4990) = 1;
% s = swarmchart(x,yShapleyValue,[],MCOER(:,1),'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
% colormap parula
% % freezeColors
% % freezeColors(jicolorbar)
%
% hold on
%
% x = ones(size(ShapleyValue,1),1)*2;
% yShapleyValue = ShapleyValue(:,2);
% % c = 1:size(ShapleyValue,1);
% % c(1:4990) = 1;
% s = swarmchart(x,yShapleyValue,[],MCOER(:,2),'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
% % caxis([0.0 99])
% colormap gray
% % freezeColors
% % freezeColors(colorbar)
%
%
% x = ones(size(ShapleyValue,1),1)*3;
% yShapleyValue = ShapleyValue(:,3);
% % c = 1:size(ShapleyValue,1);
% % c(1:4990) = 1;
% s = swarmchart(x,yShapleyValue,[],MCOER(:,3),'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
% % caxis([0.0 99])
% colormap autumn
% % freezeColors
% % freezeColors(colorbar)


%% Normalize
KernelFunctionName = "Tree";
Nc = normalize(abs(MCOER),"range");
ShapleyValue = readmatrix("Shapley value All.xlsx","Sheet",string(KernelFunctionName),"Range","B2");

r = randperm(length(ShapleyValue),floor(length(ShapleyValue)*0.01));
ShapleyValue = ShapleyValue(r,:);
Nc = Nc(r,:);
figShapley = figure;
set(figShapley,"Position",[1000,891,646,447]);
% plot([0.5 3.5],[0 0],'-k')
hold on

for i = 1:size(Nc,2)-1
    
    x = ones(size(ShapleyValue,1),1)*i;
    yShapleyValue = ShapleyValue(:,i);
    
        s = swarmchart(x,yShapleyValue,95,Nc(:,i),'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.05,'MarkerEdgeColor',[0.15,0.15,0.15]);
%     s = swarmchart(x,yShapleyValue,95,Nc(:,i),'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.9);
    
    s(1,1).XJitterWidth = 0.75;
%     s(1,1).XJitter = 'randn'
    
end
% colormap spring
% colormap
c = colormap(hot);
c = rot90(c,2);
colormap(c(24:130,:));
% colormap(c);
colormap(c1)

% xticklabels({'$$ \varepsilon_\mathrm{GDL} $$','$$ \varepsilon_\mathrm{GDL} $$','$$ \phi $$'},'TickLabelInterpreter','latex')
axes1 = gca;
set(axes1,'XTick',1:3,...
    'XTickLabel',...
    {'\epsilon_{GDL}','\epsilon_{CL}','\phi'},...
    'YGrid','on','FontSize',18);
% title("$$ \mathrm{Shapley\ summary\ plot} $$","FontSize",15,"Interpreter","latex");

box on

ylabel("Shapley value","FontSize",18,"Interpreter","tex","FontName","Arial");
xlabel("Predictor","FontSize",18,"Interpreter","tex","FontName","Arial");
c = colorbar('Ticks',[0,1],...
    'TickLabels',{'low','high'})
% c.Label.String = 'Predictor value ';

set(c.Label,"String","Predictor value","FontSize",18,...
    "Position",[1.288888738268898,0.502915928732202,0],"Interpreter","tex")
print('0Fig25 Shapleyplot1','-djpeg','-r1200') 

%% ver 2
figShapley = figure;
plot([0.5 3.5],[0 0],'-k')
hold on

for i = 1:size(Nc,2)-1
    
    x = ones(size(ShapleyValue,1),1)*i;
    yShapleyValue = ShapleyValue(:,i);
    
%     s = swarmchart(x,yShapleyValue,95,Nc(:,i),'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.05,'MarkerEdgeColor',[0.15,0.15,0.15]);
    s = swarmchart(x,yShapleyValue,95,Nc(:,i),'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.9);
    
    s(1,1).XJitterWidth = 0.7;
%     s(1,1).XJitter = 'randn'
    
end
% colormap spring
% colormap
c = colormap(spring);
% c = rot90(c,2);
% colormap(c(24:130,:));
colormap(c);


% xticklabels({'$$ \varepsilon_\mathrm{GDL} $$','$$ \varepsilon_\mathrm{GDL} $$','$$ \phi $$'},'TickLabelInterpreter','latex')
axes1 = gca;
set(axes1,'TickLabelInterpreter','latex','XTick',1:3,...
    'XTickLabel',...
    {'$$ \varepsilon_\mathrm{GDL} $$','$$ \varepsilon_\mathrm{CL} $$','$$ \phi $$'},...
    'YGrid','on');
title("$$ \mathrm{Shapley\ summary\ plot} $$","FontSize",15,"Interpreter","latex");
box on

ylabel("$$ \mathrm{Shapley\ value} $$","FontSize",13,"Interpreter","latex");
xlabel("$$ \mathrm{Predictor} $$","FontSize",13,"Interpreter","latex");
c = colorbar('Ticks',[0,1],...
    'TickLabels',{'low','high'})
% c.Label.String = 'Predictor value ';

set(c.Label,"String","$$ \mathrm{Predictor\ value} $$","FontSize",13,...
    "Position",[1.288888738268898,0.502915928732202,0],"Interpreter","latex")
print('0Fig25 Shapleyplot2','-djpeg','-r1200') 

% figure
% surfc(peaks(50))
% c = colormap(hot);
% c = rot90(c,2)
% colormap(c(60:160,:))
% colorbar