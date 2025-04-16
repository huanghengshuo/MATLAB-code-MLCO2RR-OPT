% Pearson & Spearman ori
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

phix = unique(MFE(:,3))';
rhoP = partialcorr(MFE,'Type','Pearson');

rhoS = partialcorr(MFE,'Type','Spearman');

xvalues = {'GDL','CL','phi','FECO'};
yvalues = xvalues;
f1 = figure;
h1 = heatmap(f1,xvalues,yvalues,rhoP,'GridVisible','off','FontName',...
    'Arial','FontSize',12,'CellLabelFormat','%0.2f','ColorLimits',[-0.3 1]);
h1.XDisplayLabels = {'x_{GDL}','x_{CL}','\phi','FE'}


title('Pearson partial correlation');
print('0Fig01b partialcorrPearson','-djpeg','-r1200') 

f2 = figure;
h2 = heatmap(f2,xvalues,yvalues,rhoS,'GridVisible','off','FontName',...
    'Times New Roman','FontSize',12,'CellLabelFormat','%0.2f','ColorLimits',[-0.3 1]);
title('Spearman partial correlation');
print('0Fig01b partialcorrSpearman','-djpeg','-r1200') 
h2.XDisplayLabels  


%%
rhoP = corr(MFE,'Type','Pearson');
rhoS = corr(MFE,'Type','Spearman');
f3 = figure;
h3 = heatmap(f3,xvalues,yvalues,rhoP,'GridVisible','off','FontName',...
    'Times New Roman','FontSize',12,'CellLabelFormat','%0.2f','ColorLimits',[-0.3 1]);
title('Pearson correlation');
print('0Fig01b corrPearson','-djpeg','-r1200') 


f4 = figure;
h4 = heatmap(f4,xvalues,yvalues,rhoS,'GridVisible','off','FontName',...
    'Times New Roman','FontSize',12,'CellLabelFormat','%0.2f','ColorLimits',[-0.3 1]);
title('Spearman partial correlation');
print('0Fig01b corrSpearman','-djpeg','-r1200') 






%% Train tree model
% X1 = MFE(:,1:3);
% Y1 = MFE(:,4);
% 
% % NoKFold = 10;
% TreeModelCV = fitrtree(X1,Y1);
% % % GaussModelCV = fitrgp(X1,Y1,'CrossVal','on','KFold',NoKFold);
% %
% fig = figure;
% fig.Position = [185,51,1271,912];
% p = plotPartialDependence(TreeModelCV,[1,3]);
% xlabel('x');
% ylabel('y');
% view([89.9 90]);
% xlim([min(MFE(:,1)),max(MFE(:,1))]);
% ylim([min(MFE(:,3)),max(MFE(:,3))]);
% colormap(c)
% p.Children.EdgeColor = 'none';
% p.Children.FaceColor = 'none';
% 
% title([],'Visible','off');
% 
% figure
% p = plotPartialDependence(TreeModelCV,3,"Conditional","absolute");
% p.Children(2,1).MarkerFaceColor = 'none';
% p.Children(2,1).MarkerEdgeColor = 'none';
% f2 = figure;
% %
% p = partialDependence(TreeModelCV,[1,3]);
% h = heatmap(f2,p,'GridVisible','off');

%% Train Linear Regression Model
epiGDL = MFE(:,1);
epiCL = MFE(:,2);
phi = MFE(:,3);
FE = MFE(:,4);
T = table(epiGDL,epiCL,phi,FE,'VariableNames',{'epiGDL','epiCL','phi','FE'});
lm = fitlm(T,'linear');
figure
p = plotPartialDependence(lm,3,"Conditional","absolute");
p.Children(2,1).MarkerFaceColor = 'none';
p.Children(2,1).MarkerEdgeColor = 'none';
[p.Children(3:end,1).Color] = deal([0.8, 0.8, 0.8]);
[ypred,yci] = predict(lm,MFE(:,1:3));

figure
plotmatrix(MFE)


figure
corrplot(MFE,"Type","Spearman",TestR="on")




% figure
% hold on
% plot(phi,FE,'.');
% plot(phi',ypred','r-');
% fill(phi',yci(:,1)');
% fill(phi',yci(:,2)');


% fig = figure;
% fig.Position = [185,51,1271,912];
% t = tiledlayout(3,3);
% for i = 1:3
%     for j = 1:3
%         if i <= j
%             if i == j
%                 nexttile
%                 hold on
%                 PDP = plotPartialDependence(lm,i,"Conditional","absolute");
%                 PDP.Children(2,1).MarkerFaceColor = 'none';
%                 PDP.Children(2,1).MarkerEdgeColor = 'none';
%                 [PDP.Children(3:end,1).Color] = deal([0.8, 0.8, 0.8]);
%                 scatter(MFE(:,i),MFE(:,4),'o',...
%                     'MarkerFaceColor',[0.3010 0.7450 0.9330],...
%                     'MarkerEdgeColor',[1 1 1],...
%                     'MarkerFaceAlpha',0.5,...
%                     'MarkerEdgeAlpha',0.5);
%             else
%                 nexttile
%                 p = plotPartialDependence(lm,[i,j]);
%                 %                 xlabel('x');
%                 %                 ylabel('y');
%                 view([89.9 90]);
%                 xlim([min(MFE(:,i)),max(MFE(:,i))]);
%                 ylim([min(MFE(:,j)),max(MFE(:,j))]);
%                 colormap(c)
%                 p.Children.EdgeColor = 'none';
%             end
%         else
%             nexttile
%             plot(MFE(:,j),MFE(:,i),'o','Color',[0.3010 0.7450 0.9330]);
%         end
%     end
% end








toc
