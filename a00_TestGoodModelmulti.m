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

% To remove some phi value 
% for i = 2:2:length(phix)-1
%     
%     index = find(abs(MFE(:,3)-phix(i))<0.001);
%     MFE(index,:) = [];
%     MCOER(index,:) = [];
% end
%% plot swarmchart
figure
hold on
yyaxis left
set(gca,'YColor','k')
X = MFE(:,1:end-1);
Y = MFE(:,end);
ystd = zeros(1,length(phix));
cmap = colormap(autumn(111));
% for i = 1:5:length(phix)
%     
%     index = find(abs(MFE(:,1)-phix(i))<0.001);
%     
%     X1 = X(index,1:1);
%     Y1 = Y(index,:);
%     Y1 = Y1 - mean(Y1);
%     
%     x = ones(size(Y1,1),1)*phix(i)*1;
%     x = ones(size(Y1,1),1)*i;
%     s = swarmchart(x,Y1,5,cmap(i,:),'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
%     s.XJitterWidth = 0.05;
% end

%% Prepare the data for Machine Learning
YmeanInd = zeros(size(phix));
Yminusmean = zeros(size(Y));

for i = 1:1:length(phix)
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    
    X1 = X(index,1:1);
    Y1 = Y(index,:);
    ystd(1,i) = std(Y1);
    YmeanInd(1,i) = mean(Y1);
    Yminusmean(index,:) = Y1 - mean(Y1);
end
yyaxis right
set(gca,'YColor','k')
plot(phix,ystd,'k-','LineWidth',2)
% ylim([-max(ystd)*1.0 max(ystd)*1.05])

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


fig = figure;
% fig.Color = 'none';
hold on

picture1 = plot(X(:,end)+0.05916*pH,Y,'k.');
picture2 = plot(phix,YmeanInd,'r-','LineWidth',1.2);
picturea = scatter(phi_refa,FE_refa,'SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','m',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
% picture0 = scatter(phi_ref,FE_ref,'SizeData',50,...
%     'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','r',...
%     'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);

ylim([0 1])
xlabel("$$ \phi \ \mathrm{vs\ RHE}$$","FontSize",15,"Interpreter","latex");
ylabel("$$ FE_\mathrm{CO} $$","FontSize",15,"Interpreter","latex");
legend([picture1 picture2 picturea],[...
    "$$ \mathrm{Theoretical\ data } $$",...
    "$$ \mathrm{Theoretical\ mean}$$",...
    "$$ \mathrm{Experimental\ data} $$"],"FontSize",10,"Interpreter","latex",...
    "Location","southeast");
print('0Fig01 Test good comsol model','-djpeg','-r1200') 
ax = gca;
% exportgraphics(fig,'0Fig01 Test good comsol model export.jpg','BackgroundColor','none')


fig = figure;
% fig.Color = 'none';
hold on

picture1 = plot(X(:,end)+0.05916*pH,Y,'k.');
picture2 = plot(phix,YmeanInd,'r-','LineWidth',2);
picturea = scatter(phi_refa,FE_refa,'SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','m',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
picture0 = scatter(phi_ref,FE_ref,'SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','r',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
xlabel("$$ \phi \ \mathrm{vs\ RHE}$$","FontSize",15,"Interpreter","latex");
ylabel("$$ FE_\mathrm{CO} $$","FontSize",15,"Interpreter","latex");
legend([picture1 picture2 picturea picture0],[...
    "$$ \mathrm{Theoretical\ data } $$",...
    "$$ \mathrm{Theoretical\ mean}$$",...
    "$$ \mathrm{Experimental\ data\ ref more} $$",...
    "$$ \mathrm{Experimental\ data\ ref} $$"],"FontSize",10,"Interpreter","latex",...
    "Location","southeast");
% 
figure
hold on
scatter(phi_ref,jCO,'*');
plot(X(:,end)+0.05916*pH,-MCOER(:,end),'k-');




