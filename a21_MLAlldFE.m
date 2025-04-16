%% Step 1

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

%% Gauss process & Tree regression Machine Learning for All phi

Y = Yminusmean;
KernelFunctionName = {'ardsquaredexponential';...
    'ardexponential';...
    'ardmatern32';...
    'ardmatern52';...
    'ardrationalquadratic'};
NoKFold = 20;
% for j = 1:length(KernelFunctionName)+1
for j  = 6:6
% for j  = 1:1
    disp(j)
    if j <= length(KernelFunctionName)
        WeightGauss = zeros(NoKFold+1,3);
        %         WeightTree     = WeightGauss;
        WeightGaussCV10All= zeros(NoKFold+1,3);
        
        X1 = X(:,1:3);
        Y1 = Y(:,:);
        
        %% Train Gauss model with different Kernel functions
        
        GaussModel = fitrgp(X1,Y1, ...
            'KernelFunction',KernelFunctionName{j});
        GaussModelCV = fitrgp(X1,Y1,...
            'KernelFunction',KernelFunctionName{j},'CrossVal','on','KFold',NoKFold);
        
        %% KFold Cross-Validation
        for k = 1:NoKFold
            WeightGaussCVrecent = GaussModelCV.Trained{k,1}.KernelInformation.KernelParameters(1:3);
            WeightGaussCVrecent = exp(-WeightGaussCVrecent);
            WeightGaussCVrecent = WeightGaussCVrecent/sum(WeightGaussCVrecent);
            WeightGaussCV10All(k,:) = WeightGaussCVrecent';
        end
        CVname = string(KernelFunctionName{j})+"All.mat";
        ModelCVName= "RGPCV"+string(KernelFunctionName{j})+".mat";
        ModelName= "RGP"+string(KernelFunctionName{j})+".mat";
        save(CVname,"WeightGaussCV10All")
        save(ModelCVName,"GaussModelCV")
        save(ModelName,"GaussModel")
        %% Gauss weight
                    sigmaL = GaussModel.KernelInformation.KernelParameters(1:3); % Learned length scales
                    weights = exp(-sigmaL); % Predictor weights
                    weights = weights/sum(weights); % Normalized predictor weights
        
                    tbl_weight = table(weights,'VariableNames',{'Predictor Weight'}, ...
                        'RowNames',GaussModel.ExpandedPredictorNames);
                    % tbl_weight = sortrows(tbl_weight,'Predictor Weight');
                    WeightGaussCV10All(k+1,:) = weights';
        % plot
        % b = barh(categorical(tbl_weight.Row,tbl_weight.Row),tbl_weight.('Predictor Weight'));
        % b.Parent.TickLabelInterpreter = 'none';
        % xlabel('Predictor Weight')
        % ylabel('Predictor')
        %
        
        
        varnames = ["epi_GDL","epi_CL","phi"];
        KernelFunctionNamej = string(KernelFunctionName{j});
        writematrix(varnames,"Weights ALL.xlsx","Sheet",KernelFunctionNamej,"Range","B1")
        writematrix(WeightGaussCV10All,"Weights ALL.xlsx","Sheet",KernelFunctionNamej,"Range","B2")
    else
        WeightTreeCV10ALL= zeros(NoKFold+1,3);
        
        X1 = X(:,1:3);
        Y1 = Y(:,:);
        
        %% Train Tree model and CV
%         TreeModel = fitrtree(X1,Y1);
        TreeModel = fitrtree(X1,Y1,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus','KFold',NoKFold));

%         TreeModel = fitrtree(X1,Y1,'OptimizeHyperparameters','auto');
        imp_each = predictorImportance(TreeModel);
        
        TreeModelCV = fitrtree(X1,Y1,'CrossVal','on','KFold',NoKFold);
        
        for k = 1:NoKFold
            WeightTreeCVrecent = predictorImportance(TreeModelCV.Trained{k,1});
            
            WeightTreeCV10ALL(k,:) = WeightTreeCVrecent;
        end
%         save("WeightTreeCV10All.mat","WeightTreeCV10All")
        %% Tree Importance
        varnames = ["epi_GDL","epi_CL","phi"]
        WeightTreeCV10ALL(k+1,:) = imp_each;
        writematrix(WeightTreeCV10ALL,"Weights ALL.xlsx","Sheet","Tree","Range","B2")
        writematrix(varnames,"Weights ALL.xlsx","Sheet","Tree","Range","B1")
        
        CVname = "Tree"+"All.mat";
        ModelCVName= "TreeCV"+".mat";
        ModelName= "Tree"+".mat";
        save(CVname,"WeightTreeCV10ALL")
        save(ModelCVName,"TreeModelCV")
        save(ModelName,"TreeModel")
        
        
        
        
        
        
    end
end
fig4 = figure(4)
fig4.Position = [1000,918,560,420];
fig4.Children(2,1).Children(1,1).LineWidth = 2.5;
fig4.Children(2,1).Children(2,1).LineWidth = 2.5;
ax1 = gca;
ax1.Title.Visible = 'off';
ax1.FontSize = 20;
box on
ylabel('Estimated value')
xlabel('Min Leaf Size')
legend('Position',[0.280892865897289 0.70896825513906 0.544642843625375 0.158333328933943]);

print('0Fig21 Min obj','-djpeg','-r1200') 





fig5 = figure(5)
fig5.Position = [1000,918,560,420];
ylabel('Min objective')
xlabel('Function evaluations')

fig5.Children(2,1).Children(1,1).Marker = 'p';
fig5.Children(2,1).Children(1,1).MarkerSize = 10;
fig5.Children(2,1).Children(1,1).MarkerEdgeColor = 'r';
fig5.Children(2,1).Children(1,1).MarkerFaceColor = 'r';


fig5.Children(2,1).Children(3,1).LineWidth = 1.5;
fig5.Children(2,1).Children(4,1).LineWidth = 1.5;
fig5.Children(2,1).Children(5,1).LineWidth = 1.5;
fig5.Children(2,1).Children(6,1).LineWidth = 1.5;
fig5.Children(2,1).Children(7,1).LineWidth = 1.5;




ax1 = gca;
ax1.Title.Visible = 'off';
ax1.FontSize = 18;
box on

print('0Fig21 Estimated obj 01','-djpeg','-r1200') 
xlim([1 100])
ylim([-2e-6 10e-6])
print('0Fig21 Estimated obj 02','-djpeg','-r1200') 


load c.mat
c = c1;
fig = figure;
fig.Position = [185,51,1271,912];
t = tiledlayout(3,3);
% for i = 1:3
%     for j = 1:3
%         if i <= j
%             if i == j
%                 nexttile
%                 hold on
%                 PDP = plotPartialDependence(TreeModel,i,"Conditional","absolute");
%                 PDP.Children(2,1).MarkerFaceColor = 'none';
%                 PDP.Children(2,1).MarkerEdgeColor = 'none';
%                 [PDP.Children(3:end,1).Color] = deal([0.8, 0.8, 0.8]);
%                 scatter(MFE(:,i),Y,'o',...
%                     'MarkerFaceColor',[0.3010 0.7450 0.9330],...
%                     'MarkerEdgeColor',[1 1 1],...
%                     'MarkerFaceAlpha',0.5,...
%                     'MarkerEdgeAlpha',0.5);
%             else
%                 nexttile
%                 p = plotPartialDependence(TreeModel,[i,j]);
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

%%

figPDP = figure;
xname = ["porosity of GDL";"porosity of CL";"Potential (V \it vs. \rm NHE)"]
PDP = plotPartialDependence(TreeModel,1,"Conditional","absolute");
PDP.Title.Visible = 'off'
ylabel({'$$ FE-\overline{FE}$$'},'Interpreter','latex');

xlabel({xname(2,1)},'Interpreter','latex');

figPDP = figure;
hold on
figPDP.Position = [330.6,167.4,1372,761.6];
% scatter(MFE(:,1),Y,'o',...
%     'MarkerFaceColor',[0.3010 0.7450 0.9330],...
%     'MarkerEdgeColor',[0 0 0],...
%     'MarkerFaceAlpha',0.6,...
%     'MarkerEdgeAlpha',0.2);
for i = 1:1
    %     subplot(2,3,i)
    figure
    
    PDP = plotPartialDependence(TreeModel,i,"Conditional","absolute");
    
    PDP.Children(2,1).MarkerFaceColor = [0.3010 0.7450 0.9330];
    PDP.Children(2,1).MarkerEdgeColor = [0.00,0.45,0.74];
    PDP.Children(2,1).SizeData = 50;
    PDP.Title.Visible = 'off';
    set(gca,'FontSize',18)
    ylabel('res FE','FontSize',18,'FontName','Arial');
    
    xlabel({xname(i,1)},'FontSize',18,'FontName','Arial','Interpreter','tex');
    box on
    set(gca,"XMinorTick","off")
    % PDP.Children(2,1).MarkerEdgeColor = 'none';
    % PDP.Children(2,1).MarkerFaceAlpha = 0.6;
    % PDP.Children(2,1).MarkerEdgeAlpha = 0.2;
    print(['0Fig21 PDP 1D ',num2str(i)],'-djpeg','-r1200')
    
end


%%
% pfigure = figure
% pfigure.Position = [];
% subplot(1,3,1)
fig = figure;
fig.Position = [645,285,560,543.2];
p = plotPartialDependence(TreeModel,[1 2]);
colormap(c)
% p.Children.EdgeColor = 'none';
p.Children.EdgeAlpha = 0.2;
p.Title.Visible = 'off';
xlabel({xname(1,1)},'FontSize',18,"FontName","Arial",'Interpreter','tex');
ylabel({xname(2,1)},'FontSize',18,"FontName","Arial",'Interpreter','tex');        
zlabel('res FE','FontSize',18,"FontName","Arial");
print('0Fig21 PDP 2D 1','-djpeg','-r1200') 


% subplot(1,3,2)
fig = figure;
fig.Position = [645,285,560,543.2];
p = plotPartialDependence(TreeModel,[1 3]);
colormap(c)
% p.Children.EdgeColor = 'none';
p.Children.EdgeAlpha = 0.2;
p.Title.Visible = 'off';
xlabel({xname(1,1)},'FontSize',18,"FontName","Arial",'Interpreter','tex');
ylabel({xname(3,1)},'FontSize',18,"FontName","Arial",'Interpreter','tex');        
zlabel('res FE','FontSize',18,"FontName","Arial");
print('0Fig21 PDP 2D 2','-djpeg','-r1200') 



% subplot(1,3,3)
fig = figure;
fig.Position = [645,285,560,543.2];
p = plotPartialDependence(TreeModel,[2 3]);
colormap(c)
% p.Children.EdgeColor = 'none';
p.Children.EdgeAlpha = 0.2;
p.Title.Visible = 'off';
xlabel({xname(2,1)},'FontSize',18,"FontName","Arial",'Interpreter','tex');
ylabel({xname(3,1)},'FontSize',18,"FontName","Arial",'Interpreter','tex');        
zlabel('res FE','FontSize',18,"FontName","Arial");

print('0Fig21 PDP 2D 3','-djpeg','-r1200') 




toc















