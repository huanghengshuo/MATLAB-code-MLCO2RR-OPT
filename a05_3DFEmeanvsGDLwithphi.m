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
index = find(MFE(:,1)<=MFE(:,2));
MFE(index,:)   = [];
MCOER(index,:) = [];


X = MFE(:,1:end-1);
Y = MFE(:,end);
xGDL = unique(MFE(:,1))';
xCL =  unique(MFE(:,2))';
phix = unique(MFE(:,3))';
%% Prepare the data for Machine Learning
YmeanInd = zeros(size(phix));
Yminusmean = zeros(size(Y));
index = find(MFE(:,1)<=MFE(:,2));
MFE(index,:)   = [];
MCOER(index,:) = [];

for i = 1:1:length(phix)
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    Y1 = Y(index,:);
    YmeanInd(1,i) = mean(Y1);
    Yminusmean(index,:) = Y1 - mean(Y1);
end


% xCL =  unique(MFE(:,2))';
% 
% index = find(MFE(:,2)== xCL(6));
% MFE = MFE(index,:);
% MCOER = MCOER(index,:);
% 





x = MFE(:,1);
y = MFE(:,3);
v = Yminusmean;

[xq,yq] = meshgrid(linspace(min(xGDL),max(xGDL),500), linspace(min(phix),max(phix),500));
vq = griddata(x,y,v,xq,yq);

fig = figure;
fig.Position = [742.6,552.2,628,420];
splot = contourf(xq,yq,vq,[-0.07:0.01:-0.01,-0.01:0.003:0,0:0.003:0.016,0.019]);
% splot = contourf(xq,yq,vq,'ShowText','on');

view([0 90])
set(gca,'FontSize',18,'Position',[0.135031847133758,0.157619050888788,0.672611464968153,0.767380949111212]);

c = colormap(hot);
c = rot90(c,2);
c = flipud(c);
colormap(flipud(c1));
c = colorbar
c.Position = [0.851726071488318,0.157619050888788,0.033970276008493,0.609523809523811]

% c.Label.String = "$$ \widetilde{FE} $$";
% c.Label.Interpreter = "latex";
% c.Label.FontSize = 20;

xlabel("Porosity of GDL","FontSize",18,"FontName","Arial");

% xlabel("$$ \varepsilon_\mathrm{CL} $$","Interpreter","latex","FontSize",20)
ylabel("Potential (V \it vs. \rm NHE)","FontSize",18,"FontName","Arial");
% title("$$ \mathrm{residual\ FE} $$","Interpreter","latex","FontSize",20)
% zlabel("$$ residual\ FE $$","Interpreter","latex")
annotation(fig,'textbox',...
    [0.821656050955414 0.803809523809527 0.138853503184714 0.0761904761904765],...
    'String',{'res FE'},...
    'LineStyle','none',...
    'FontSize',18,...
    'FontName','Arial',...
    'FitBoxToText','off');
print('0Fig05 dFE vs GDL with phi','-djpeg','-r1200') 


toc
