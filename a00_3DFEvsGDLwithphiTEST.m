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
% M0 = readmatrix("jCOER test.xlsx"); % GDL CL  phi
% M1 = readmatrix("jHER test.xlsx");
% M2 = readmatrix("jTotal test.xlsx");
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

xGDL = unique(MFE(:,1))';
xCL =  unique(MFE(:,2))';
phix = unique(MFE(:,3))';

% figure
% hold on


for i = 1:1:length(xGDL)
    index = find(abs(MFE(:,1)-xGDL(i))<0.001);
    
x = MFE(index,2);
y = MFE(index,3);
v = MFE(index,end);

[xq,yq] = meshgrid(linspace(min(xCL),max(xCL),500), linspace(min(phix),max(phix),500));
vq = griddata(x,y,v,xq,yq);

fig = figure;
fig.Position = [742.6,552.2,628,420];

box on
splot = contourf(xq,yq,vq,'ShowText','on');
set(gca,'FontSize',18);
% xlim([0.45 0.7])

% view([0 90])

c = colormap(hot);
c = rot90(c,2);
c = flipud(c);
colormap(c(130:end,:));
c = colorbar;
% c.Label.String = 'FE';
xlabel("Porosity of CL","FontSize",18,"FontName","Arial");
ylabel("Potential (V \it vs. \rm NHE)","FontSize",18,"FontName","Arial");
title(['xGDL = ',num2str(xGDL(i))],"FontSize",18,"FontName","Arial")
end


for i = 1:1:length(xGDL)
    index = find(abs(MFE(:,1)-xGDL(i))<0.001);
    
x = MFE(index,2);
y = MFE(index,3);
v = MFE(index,end);

[xq,yq] = meshgrid(linspace(min(xCL),max(xCL),500), linspace(min(phix),max(phix),500));
vq = griddata(x,y,v,xq,yq);

fig = figure;
fig.Position = [742.6,552.2,628,420];

box on
splot = surf(xq,yq,vq,'EdgeColor','none');
set(gca,'FontSize',18);
% xlim([0.45 0.7])

% view([0 90])

c = colormap(hot);
c = rot90(c,2);
c = flipud(c);
colormap(c(130:end,:));
c = colorbar;
% c.Label.String = 'FE';
xlabel("Porosity of CL","FontSize",18,"FontName","Arial");
ylabel("Potential (V \it vs. \rm NHE)","FontSize",18,"FontName","Arial");
title(['xGDL = ',num2str(xGDL(i))],"FontSize",18,"FontName","Arial")
end




for i = 1:1:length(xGDL)
    index = find(abs(MFE(:,1)-xGDL(i))<0.001);
    figure
x = MFE(index,2);
y = MFE(index,3);
v = MFE(index,end);

picturea = scatter(y,v,'SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','m',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
ylabel("FE","FontSize",18,"FontName","Arial");
xlabel("Potential (V \it vs. \rm NHE)","FontSize",18,"FontName","Arial");
title(['xGDL = ',num2str(xGDL(i))],"FontSize",18,"FontName","Arial")
end












