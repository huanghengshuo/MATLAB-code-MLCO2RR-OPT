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
% index = find(MFE(:,2)== xCL(3));
% MFE = MFE(index,:);
% MCOER = MCOER(index,:);
% 
index = find(MFE(:,1)<=MFE(:,2));
MFE(index,:)   = [];
MCOER(index,:) = [];
% 

xGDL = unique(MFE(:,1))';
xCL =  unique(MFE(:,2))';
phix = unique(MFE(:,3))';



x = MFE(:,1);
y = MFE(:,3);
v = MFE(:,end);

[xq,yq] = meshgrid(linspace(min(xGDL),max(xGDL),500), linspace(min(phix),max(phix),500));
vq = griddata(x,y,v,xq,yq);

fig = figure;
fig.Position = [742.6,552.2,628,420];

splot = contourf(xq,yq,vq,[0.70:0.05:0.95,0.955]);
set(gca,'FontSize',18,'Position',[0.135031847133758,0.157619050888788,0.672611464968153,0.767380949111212]);

view([0 90])

c = colormap(hot);
c = rot90(c,2);
c = flipud(c);
colormap(flipud(c1));
c = colorbar;
c.Position = [0.851726071488318,0.157619050888788,0.033970276008493,0.609523809523811]

xlabel("Porosity of GDL","FontSize",18,"FontName","Arial");

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

print('0Fig05 FE vs GDL with phi','-djpeg','-r1200') 

toc
