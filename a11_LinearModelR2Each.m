% Step 1
clc
clear all
close all
rng('default') % For reproducibility
%   Seperated analysis for ALL
%   Seperated analysis for ALL
%   Seperated analysis for ALL
tic
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

% index = find(MFE(:,3)>-0.6);
% MFE(index,:)   = [];
% MCOER(index,:) = [];

phix = unique(MFE(:,3))';
for i = 2:2:length(phix)-1
    
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    MFE(index,:) = [];
    MCOER(index,:) = [];
end
phix = unique(MFE(:,3))';
for i = 2:2:length(phix)-1
    
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    MFE(index,:) = [];
    MCOER(index,:) = [];
end
phix = unique(MFE(:,3))';
for i = 2:2:length(phix)-1
    
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    MFE(index,:) = [];
    MCOER(index,:) = [];
end
phix = unique(MFE(:,3))';
for i = 2:2:length(phix)-1
    
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    MFE(index,:) = [];
    MCOER(index,:) = [];
end

xGDL = unique(MFE(:,1))';
xCL  = unique(MFE(:,2))';
phix = unique(MFE(:,3))';
%% plot swarmchart
figure
hold on
yyaxis left

set(gca,'YColor','k')
X = MFE(:,1:end-1);
Y = MFE(:,end);
axes1 = gca;
ystd = zeros(1,length(phix));
cmap = colormap(cool(length(phix)));
% b = boxplot(MFE(:,4),num2str(MFE(:,3)),'ColorGroup',num2str(MFE(:,3)),'Jitter',0.5);
% axes1 = gca;
% axes1.XDir = "Reverse";

b = boxchart(MFE(:,3),MFE(:,4),'BoxWidth',0.07);
b.JitterOutliers = 'off';
b.MarkerStyle = 'none';
meanFE = groupsummary(MFE(:,4),MFE(:,3),'mean');
for i = 1:1:length(phix)
    
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    
    X1 = X(index,1:2);
    Y2 = Y(index,:);
    Y1 = Y2 - mean(Y2);
    
    x = ones(size(Y2,1),1)*phix(i)*1;
    % x = ones(size(Y1,1),1)*i;
    s = swarmchart(axes1,x,Y2,20,cmap(i,:),'o','filled',...
        'MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8);
    s.XJitterWidth = 0.02;
end
hold on
pmean = plot(phix,meanFE,'k-p');
xlim([-2.05 -0.55])

xlabel("$$ \phi \ \mathrm{vs\ NHE}$$","FontSize",15,"Interpreter","latex");
ylabel("$$ FE_\mathrm{CO} $$","FontSize",15,"Interpreter","latex")
%% Prepare the data for Machine Learning
YmeanInd = zeros(size(phix));
Yminusmean = zeros(size(Y));
slopeCL = YmeanInd;
slopeGDL= YmeanInd;
R2      = YmeanInd;
rsquare = YmeanInd;
for i = 1:1:length(phix)
    % for i = 1:1:1
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    
    X1 = X(index,1:2);
    Y2 = Y(index,:);
    Y1 = Y2 - mean(Y2);
    ystd(1,i) = std(Y1);
    YmeanInd(1,i) = mean(Y1);
    Yminusmean(index,:) = Y1 - mean(Y1);
    
    figure
    %     scatter3(X1(:,1),X1(:,2),Y1,'filled');
    stem3(X1(:,1),X1(:,2),Y1,'filled');
    hold on
    
    
    XX = [ones(size(X1,1),1) X1];
    lm = fitlm(X1,Y1);
    [b,bint,r,rint,stats] = regress(Y1,XX);
    slopeGDL(i) = b(2);
    slopeCL(i)  = b(3);
    R2(i)       = stats(1);
    
    xGDLfit = linspace(min(xGDL),max(xGDL),100);
    xCLfit  = linspace(min(xCL),max(xCL),100);
    [X1FIT,X2FIT] = meshgrid(xGDLfit,xCLfit);
    YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
    surf(X1FIT,X2FIT,YFIT,'FaceAlpha',0.8,'EdgeColor','none')
    set(gca,'FontSize',16)
    title(['\phi \rm = ',num2str(phix(i)),' (V \it vs. \rm NHE)'],"FontSize",18,"Interpreter","tex");
    xlabel("Porosity of GDL","FontSize",16,"FontName","Arial");
    ylabel("Porosity of CL","FontSize",16,"FontName","Arial");
    zlabel("res FE","FontSize",16,"FontName","Arial");
    view([-40 10])
    zlim([min(YFIT,[],"all") max(YFIT,[],"all")])
    %     zlim([min(Y1,[],"all")*0.99 max(Y1,[],"all")*1.01])
    colormap(c1)
    %     legend(["Data","Linear model"],'FontSize',18,...
    %         'Position',[0.666544727892742,0.8361904783476,0.228373319076158,0.086904759634109])
    %     print(['0Fig11 scatter & LinearModel',num2str(i)],'-djpeg','-r1200')
    
    
    %     TreeModel = fitrtree(X1,Y1,'OptimizeHyperparameters','auto',...
    %         'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    %         'expected-improvement-plus','CVPartition',cvpartition(length(X1),"KFold",NoKFold)));
    %     imp_each = predictorImportance(TreeModel);
    NoKFold = 10;
    TreeModelCV = fitrtree(X1,Y1,'CrossVal','on','KFold',NoKFold);
    
    f = @(CMP,Xtrain,Ytrain,Wtrain,Xtest,Ytest,Wtest)...
        predict(CMP,X1);
    ypredx = kfoldfun(TreeModelCV,f);
    ypredCV = mean(ypredx)';
    SSE = sum((ypredCV      -       Y1).^2);
    SST = sum((Y1 - mean(Y1)).^2);
    rsquare(i) = 1 - SSE./SST;
      
end

%% Bar plot
fig = figure
fig.Position = [745,409.8,585.6,420];
y = [slopeGDL; slopeCL];
barplot = bar(phix(1:end-1),y(:,1:end-1),'stacked','BarWidth',0.8);

barplot(1).FaceColor = [0.30,0.75,0.93];
barplot(2).FaceColor = [0.00,0.45,0.74];
set(gca,'FontSize',18)
% xlabel("$$ \phi \ \mathrm{vs\ NHE}$$","FontSize",20,"Interpreter","latex");
% ylabel("$$ \mathrm{slope\ of\ linear\ model}$$","FontSize",20,"Interpreter","latex");
xlabel("Potential (V \it vs. \rm NHE)","FontSize",18,"FontName","Arial");
ylabel("slope of linear model","FontSize",18,"FontName","Arial");


yyaxis right
hold on

plot(phix,R2,'o-','LineWidth',2.5,'MarkerSize',8,'Color',[0.36,0.56,1.00],'MarkerEdgeColor',[0.36,0.56,1.00],'MarkerFaceColor',[0.36,0.56,1.00])

% plot(phix,R2,'o-','LineWidth',1.5,'Color',[0.36,0.56,1.00])
% plot(phix,rsquare,'o-','LineWidth',1.5,'Color',[0.00,0.56,1.00])
set(gca,'YColor','k','YTick',0.3:0.1:1,'YLim',[0.3 1.01])

ylabel("r^2","FontSize",18,"FontName","Arial");


legend(["\beta_{GDL} ","\beta_{CL}","r^2"],'FontSize',18,'Position',[0.661309508425035,0.606507961333744,0.172500001941408,0.250000005676633],'EdgeColor','none')
print('0Fig11 slope r2','-djpeg','-r1200')



%% bar plot for linear model and DTr model

fig = figure
box on
hold on
fig.Position = [745,409.8,585.6,420];
y = [R2; rsquare];
% barplot = bar(phix(1:end-1),y(:,1:end-1),'BarWidth',0.8);

% barplot(1).FaceColor = [0.30,0.75,0.93];
% barplot(2).FaceColor = [0.00,0.45,0.74];

plot(phix,R2,'o-','LineWidth',2.5,'MarkerSize',8,'Color',[0.30,0.75,0.93],'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93])
plot(phix,rsquare,'^-','LineWidth',2.5,'MarkerSize',8,'Color',[0.00,0.45,0.74],'MarkerEdgeColor',[0.00,0.45,0.74],'MarkerFaceColor',[0.00,0.45,0.74])

ylim([0 1])
set(gca,'FontSize',18)
xlabel("Potential (V \it vs. \rm NHE)","FontSize",18,"FontName","Arial");
ylabel("r^2","FontSize",18,"FontName","Arial");
legend(["Linear","ML"],'FontSize',18,'Position',[0.663093183584991 0.686228919185047 0.182593853807286 0.134523805833998],'EdgeColor','none')
print('0Fig11 slope r2 Linear vs ML ','-djpeg','-r1200')



% Step 2

rng('default') % For reproducibility
%   Seperated analysis for Each
%   Seperated analysis for Each
%   Seperated analysis for Each
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


XX = [ones(size(X,1),1) X];
lm = fitlm(X,Yminusmean);
[b,bint,r,rint,stats] = regress(Yminusmean,XX);
slopeGDL = b(2);
slopeCL  = b(3);
R2      = stats(1);



figure2 = figure
YFIT = XX*b;
plot(Yminusmean,YFIT,'*')
box on
hold on
%     grid minor
scatter(Yminusmean,YFIT,'SizeData',70,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor',[0.07,0.62,1.00],...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5)
%     plot(Yminusmean,ypred,'ko','MarkerSize',5,'MarkerFaceColor','c',"MarkerIndices",1:1:length(ypred))
plot([min(Yminusmean) max(Yminusmean)],[min(Yminusmean) max(Yminusmean)],'k-','LineWidth',2)
xlim([min(Yminusmean) max(Yminusmean)])
ylim([min(Yminusmean) max(Yminusmean)])
set(gca,'FontSize',18)
annotation(figure2,'textbox',...
    [0.52778571428571,0.235714289333143,0.250198466324851,0.145238091619242],...
    'String',...
    ['r^2=',char(num2str(R2,'%.6f')),''],...
    'FontSize',18,...
    'FitBoxToText','off',...
    'FontName','Arial',...
    'EdgeColor','none');
xlabel("True response ","FontSize",18,"FontName","Arial");
ylabel("Predicted response","FontSize",18,"FontName","Arial");

print(['0Fig11 scatter True vs Predicted ','Linear model'],'-djpeg','-r1200')









