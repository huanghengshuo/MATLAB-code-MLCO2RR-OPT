clc
clear all
close all
rng('default') % For reproducibility
%   Seperated analysis for ALL
%   Seperated analysis for ALL
%   Seperated analysis for ALL
tic
global s
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

% phix = -0.600:-0.005:-2.400;

% phix = unique(MFE(:,3))';
% for i = 2:2:length(phix)-1
%
%     index = find(abs(MFE(:,3)-phix(i))<0.001);
%     MFE(index,:) = [];
%     MCOER(index,:) = [];
% end

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
for i = 1:2:length(phix)
    
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    
    X1 = X(index,1:3);
    Y1 = Y(index,:);
    si = size(X1,1);
    index = 1:floor(si/30):si;
    bubblechart(X1(index,3),Y1(index),X1(index,1),'LineWidth',0.1,...
        'MarkerEdgeColor',[0.3010 0.7450 0.9330],...
        'MarkerFaceColor',[0.3010 0.7450 0.9330],...
        'MarkerEdgeAlpha',0.9,'MarkerFaceAlpha',0.2);
    bubblesize([3 7])
    scatter(X(:,3),Y,[],X(:,3),'LineWidth',0.2,...
        'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5)
    scatter(X1(:,3),Y1,rescale(X1(:,1),0.1,1).*150,...
        'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.3010 0.7450 0.9330],...
        'MarkerEdgeAlpha',0.7,'MarkerFaceAlpha',0.6)
    % plot(X(:,3),Y,'k.')
end
pmean = plot(phix,YmeanInd,'k-','LineWidth',2.5)

blgd = bubblelegend('$$ \varepsilon_\mathrm{GDL} $$','Interpreter','latex');
blgd.Location = 'east';
blgd.NumBubbles = 3;
blgd.FontSize = 15;
legend(pmean,'$$ \overline{FE} $$','location','northwest',...
    'FontSize',15,'Interpreter','latex');

%%
figure
hold on
for i = 1:2:length(phix)
    
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    
    X1 = X(index,1:3);
    Y1 = Yminusmean(index,:);
    si = size(X1,1);
    index = 1:floor(si/30):si;
    bubblechart(X1(index,3),Y1(index),X1(index,1),'LineWidth',0.1,...
        'MarkerEdgeColor',[0.3010 0.7450 0.9330],...
        'MarkerFaceColor',[0.3010 0.7450 0.9330],...
        'MarkerEdgeAlpha',0.9,'MarkerFaceAlpha',0.2);
    bubblesize([3 7])
    % scatter(X(:,3),Y,[],X(:,3),'LineWidth',0.2,...
    %     'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5)
    %     scatter(X1(:,3),Y1,rescale(X1(:,1),0.1,1).*150,...
    %         'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.3010 0.7450 0.9330],...
    %         'MarkerEdgeAlpha',0.7,'MarkerFaceAlpha',0.6)
    % plot(X(:,3),Y,'k.')
end
pmean = plot(phix,YmeanInd-YmeanInd,'k-','LineWidth',2.5);
blgd = bubblelegend('$$ \varepsilon_\mathrm{GDL} $$','Interpreter','latex');
blgd.Location = 'east';
blgd.NumBubbles = 3;
blgd.FontSize = 15;
legend(pmean,'$$ FE-\overline{FE} $$','location','northwest',...
    'FontSize',15,'Interpreter','latex');






figure
hold on
scatter(X(:,3),Yminusmean,X(:,1).*80,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.4470 0.7410],...
    'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.2)

plot(phix,YmeanInd-YmeanInd,'k-','LineWidth',2.5);

% figure
% plotmatrix(X(:,3),Yminusmean)


%% Determine the Performance of ALL Machine Learning
KernelFunctionName = {'ardsquaredexponential';...
    'ardexponential';...
    'ardmatern32';...
    'ardmatern52';...
    'ardrationalquadratic'};
Performance = zeros(6,6); %% r2 r2CV RMSE RMSECV MAE MAECV




for j = 1:length(KernelFunctionName)
% for j = 1:1
    ModelCVName= "RGPCV"+string(KernelFunctionName{j})+".mat";
    ModelName= "RGP"+string(KernelFunctionName{j})+".mat";
    load(ModelName)
    load(ModelCVName)
    % ------------  R2 RMSE MAE
    ypred = resubPredict(GaussModel);
    SSE = sum((ypred      -       Yminusmean).^2);
    SST = sum((Yminusmean - mean(Yminusmean)).^2);
    rsquare = 1 - SSE./SST;
    RMSE = rms(ypred - Yminusmean);
    MAE  = mean(abs(ypred - Yminusmean));
    Performance(j,[1 3 5]) = [rsquare RMSE MAE];
    
    figure1 = figure;
    box on
    hold on
%     grid minor
    scatter(Yminusmean,ypred,'SizeData',70,...
        'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor',[0.07,0.62,1.00],...
        'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5)
    %     plot(Yminusmean,ypred,'ko','MarkerSize',5,'MarkerFaceColor','c',"MarkerIndices",1:1:length(ypred))
    plot([min(Yminusmean) max(Yminusmean)],[min(Yminusmean) max(Yminusmean)],'k-','LineWidth',2)
    xlim([min(Yminusmean) max(Yminusmean)])
    ylim([min(Yminusmean) max(Yminusmean)])
    set(gca,'FontSize',18)
xlabel("True response ","FontSize",18,"FontName","Arial");
ylabel("Predicted response","FontSize",18,"FontName","Arial");

annotation(figure1,'textbox',...
    [0.52778571428571,0.235714289333143,0.250198466324851,0.145238091619242],...
    'String',...
    {['r^2=',char(num2str(rsquare,'%.6f')),''],...
    ['MAE=',char(num2str(MAE,'%.10f')),''],...
    ['RMSE=',char(num2str(RMSE,'%.10f')),'']},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'FontName','Arial',...
    'EdgeColor','none');
    print(['0Fig23 scatter True vs Predicted ',char("RGP "+string(KernelFunctionName{j}))],'-djpeg','-r1200')

    % ----------- residual
    figres = figure;
    figres.Position = [2,341,1908,420];
    Res = Yminusmean - ypred;
    plot(Res,'o','MarkerFaceColor',[0 0.447058823529412 0.741176470588235],...
        'MarkerEdgeColor',[0 0.447058823529412 0.741176470588235])
    title(string(KernelFunctionName{j})+" residual")
    
    
    
    f = @(CMP,Xtrain,Ytrain,Wtrain,Xtest,Ytest,Wtest)...
        predict(CMP,X);
    
    ypredx = kfoldfun(GaussModelCV,f);
    ypredCV = mean(ypredx)';
    SSECV = sum((ypredCV    -       Yminusmean).^2);
    SSTCV = sum((Yminusmean - mean(Yminusmean)).^2);
    rsquareCV = 1 - SSECV./SSTCV;
    RMSECV = rms(ypredCV - Yminusmean);
    MAECV= mean(abs(ypredCV - Yminusmean));
    Performance(j,[2 4 6]) = [rsquareCV RMSECV MAECV];
    figure2 = figure;
    box on
    hold on
%     grid minor
    %     plot(Yminusmean,ypred,'kp','MarkerSize',12,'MarkerFaceColor','c',"MarkerIndices",1:1:length(ypred))
    scatter(Yminusmean,ypredCV,'SizeData',70,...
        'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor',[0.07,0.62,1.00],...
        'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5)
    plot([min(Yminusmean) max(Yminusmean)],[min(Yminusmean) max(Yminusmean)],'k-','LineWidth',2)
    xlim([min(Yminusmean) max(Yminusmean)])
    ylim([min(Yminusmean) max(Yminusmean)])
    set(gca,'FontSize',18)
xlabel("True response ","FontSize",18,"FontName","Arial");
ylabel("Predicted response","FontSize",18,"FontName","Arial");

annotation(figure2,'textbox',...
    [0.52778571428571,0.235714289333143,0.250198466324851,0.145238091619242],...
    'String',...
    {['r^2=',char(num2str(rsquareCV,'%.6f')),''],...
    ['MAE=',char(num2str(MAECV,'%.10f')),''],...
    ['RMSE=',char(num2str(RMSECV,'%.10f')),'']},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'FontName','Arial',...
    'EdgeColor','none');
    print(['0Fig23 scatter True vs Predicted ',char("RGP "+string(KernelFunctionName{j})+"CV")],'-djpeg','-r1200')

    % ----------- residual CV
    figres = figure;
    figres.Position = [2,341,1908,420];
    Res = Yminusmean - ypred;
    plot(Res,'o')
    title(string(KernelFunctionName{j})+" CV residual")
    
    
    %     figure
    %     hold on
    %     plot(X(:,3),Yminusmean,'.')
    %     plot(X(:,3),ypred,'.')
end


%% Analyse Tree Regression Model
CVname = "Tree"+"All.mat";
ModelCVName= "TreeCV"+".mat";
ModelName= "Tree"+".mat";
load(ModelName)
load(ModelCVName)
ypred = resubPredict(TreeModel);
% Tree
SSE = sum((ypred      -      Yminusmean ).^2);
SST = sum((Yminusmean - mean(Yminusmean)).^2);
rsquare = 1 - SSE./SST;
RMSE = rms(ypred - Yminusmean);
MAE  = mean(abs(ypred - Yminusmean));
Performance(6,[1 3 5]) = [rsquare RMSE MAE];
figure1 = figure;
hold on
% grid minor
box on
scatter(Yminusmean,ypred,'SizeData',70,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor',[0.07,0.62,1.00],...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7)
%     plot(Yminusmean,ypred,'ko','MarkerSize',5,'MarkerFaceColor','c',"MarkerIndices",1:1:length(ypred))
plot([min(Yminusmean) max(Yminusmean)],[min(Yminusmean) max(Yminusmean)],'k-','LineWidth',2)
xlim([min(Yminusmean) max(Yminusmean)])
ylim([min(Yminusmean) max(Yminusmean)])
set(gca,'FontSize',18)
% title("Decision Tree Model","FontSize",18,"FontName","Arial");
xlabel("True response ","FontSize",18,"FontName","Arial");
ylabel("Predicted response","FontSize",18,"FontName","Arial");

annotation(figure1,'textbox',...
    [0.52778571428571,0.235714289333143,0.250198466324851,0.145238091619242],...
    'String',...
    {['r^2=',char(num2str(rsquare,'%.6f')),''],...
    ['MAE=',char(num2str(MAE,'%.10f')),''],...
    ['RMSE=',char(num2str(RMSE,'%.10f')),'']},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'FontName','Arial',...
    'EdgeColor','none');
    print(['0Fig23 scatter True vs Predicted ',char("DTr")],'-djpeg','-r1200')

% ----------- residual
figres = figure;
figres.Position = [2,341,1908,420];
Res = Yminusmean - ypred;
plot(Res,'o','MarkerFaceColor',[0 0.447058823529412 0.741176470588235],...
    'MarkerEdgeColor',[0 0.447058823529412 0.741176470588235])
title("Tree Model")
xlabel("No of reponse")
ylabel("residual")

%% Analyse Tree Regression Model CV
f = @(CMP,Xtrain,Ytrain,Wtrain,Xtest,Ytest,Wtest)...
    predict(CMP,X);

ypredx = kfoldfun(TreeModelCV,f);
ypredCV = mean(ypredx)';
% Tree CV
SSECV = sum((ypredCV    -      Yminusmean ).^2);
SSTCV = sum((Yminusmean - mean(Yminusmean)).^2);
rsquareCV = 1 - SSECV./SSTCV;
RMSECV = rms(ypredCV - Yminusmean);
MAECV= mean(abs(ypredCV - Yminusmean));
Performance(6,[2 4 6]) = [rsquareCV RMSECV MAECV];
figure1 = figure;
hold on
%     plot(Yminusmean,ypred,'kp','MarkerSize',12,'MarkerFaceColor','c',"MarkerIndices",1:1:length(ypred))
scatter(Yminusmean,ypredCV,'SizeData',70,...
    'LineWidth',0.3,'MarkerEdgeColor','k','MarkerFaceColor',[0.07,0.62,1.00],...
    'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.7)
plot([min(Yminusmean) max(Yminusmean)],[min(Yminusmean) max(Yminusmean)],'k-','LineWidth',2)
xlim([min(Yminusmean) max(Yminusmean)])
ylim([min(Yminusmean) max(Yminusmean)])
box on
hold on

set(gca,'FontSize',18)
% title("Decision Tree Model","FontSize",18,"FontName","Arial");
xlabel("True response ","FontSize",18,"FontName","Arial");
ylabel("Predicted response","FontSize",18,"FontName","Arial");

annotation(figure1,'textbox',...
    [0.52778571428571,0.235714289333143,0.250198466324851,0.145238091619242],...
    'String',...
    {['r^2=',char(num2str(rsquareCV,'%.6f')),''],...
    ['MAE=',char(num2str(MAECV,'%.10f')),''],...
    ['RMSE=',char(num2str(RMSECV,'%.10f')),'']},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'FontName','Arial',...
    'EdgeColor','none');

print(['0Fig23 scatter True vs Predicted ',char("DTr CV")],'-djpeg','-r1200')
% ----------- residual CV
figres = figure;
figres.Position = [2,341,1908,420];
ResCV = Yminusmean - ypredCV;
% plot(ResCV,'o','MarkerFaceColor',[0 0.447058823529412 0.741176470588235],...
%     'MarkerEdgeColor',[0 0.447058823529412 0.741176470588235])

scatter(1:length(ResCV),ResCV,'SizeData',30,...
    'LineWidth',0.3,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.447058823529412 0.741176470588235],...
    'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.3)
title("Tree Model CV residual")
xlabel("No of reponse")
ylabel("residual")


figres = figure;
figres.Position = [2,341,1908,420];
scatter(MFE(:,3),ResCV,'SizeData',30,...
    'LineWidth',0.3,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.447058823529412 0.741176470588235],...
    'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.3)
% title("Tree Model CV residual")

xlabel("phi")
ylabel("residual")


%% ICE
% 
% figure
% subplot(2,2,1)
% plotPartialDependence(TreeModel,1,"Conditional","absolute")
% subplot(2,2,2)
% plotPartialDependence(TreeModel,2,"Conditional","absolute")
% subplot(2,2,3)
% plotPartialDependence(TreeModel,3,"Conditional","absolute")
% 
% 
% 
% fig18 = figure
% hold on
% pt = linspace(min(xGDL),max(xGDL),50)';
% % s = plotPartialDependence(TreeModelCV.Trained{1,1},1)
% PDPdata = cell(TreeModelCV.KFold,34000);
% for i = 1:20
%     s = plotPartialDependence(TreeModelCV.Trained{i,1},2,TreeModelCV.X,"Conditional","absolute");
%     
%     % p = s.Children(1); % PDP
%     pdata = [s.Children(1).XData;s.Children(1).YData];
%     PDPdata{i,1} = pdata;
%     
%     % ps= s.Children(2); % Scatter
%     psdata = [s.Children(2).XData;s.Children(2).YData];
%     PDPdata{i,2} = psdata;
%     for j = 3:size(s.Children)
%         % pr= s.Children(3); % ICE...
%         prdata = [s.Children(j).XData;s.Children(j).YData];
%         PDPdata{i,j} = prdata;
%     end
% end
% save PDPdata2.mat PDPdata
% toc
% 
% % xval = s.Children.XData;
% % yval = s.Children.YData;
% % zval = s.Children.ZData;
% figure
% p = s.Children(1); % PDP
% ps= s.Children(2); % Scatter
% pr= s.Children(3); % ICE...
% 
% plot(s.Children(33000).XData,s.Children(33000).YData,'o-')
%% Test
% ypred = resubPredict(trainedModel.RegressionTree);
% SSRCV = sum((ypredCV    - mean(Yminusmean)).^2);
% SSTCV = sum((Yminusmean - mean(Yminusmean)).^2);
% rsquareCV = SSRCV./SSTCV;

% function y = Myfun(j)
% global s
% y = [s.Children(j).XData;s.Children(j).YData];
% end