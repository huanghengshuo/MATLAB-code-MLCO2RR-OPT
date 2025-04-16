% Step 2
clc
clear all
close all
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

%% Gauss process & Tree regression Machine Learning for each phi

Y = Yminusmean;
KernelFunctionName = {'ardsquaredexponential';...
    'ardexponential';...
    'ardmatern32';...
    'ardmatern52';...
    'ardrationalquadratic'};
NoKFold = 20;
% for j = 1:length(KernelFunctionName)+1
for j = 6:6
    % for j = 1:1
    if j <= length(KernelFunctionName)
        WeightGauss = zeros(length(phix),2);
        WeightTree     = WeightGauss;
        WeightGaussCV10= zeros(length(phix),2,NoKFold);
        for i = 1:length(phix)
            
            % for i = 1:10
            index = find(abs(MFE(:,3)-phix(i))<0.001);
            
            X1 = X(index,1:2);
            Y1 = Y(index,:);
            
            %% Train Gauss model with different Kernel functions
            
            GaussModel = fitrgp(X1,Y1, ...
                'KernelFunction',KernelFunctionName{j});
            GaussModelCV = fitrgp(X1,Y1,...
                'KernelFunction',KernelFunctionName{j},'CrossVal','on','KFold',NoKFold);
            
            %% KFold Cross-Validation
            for k = 1:NoKFold
                WeightGaussCVrecent = GaussModelCV.Trained{k,1}.KernelInformation.KernelParameters(1:2);
                WeightGaussCVrecent = exp(-WeightGaussCVrecent);
                %                 WeightGaussCVrecent = WeightGaussCVrecent/sum(WeightGaussCVrecent);
                WeightGaussCV10(i,:,k) = WeightGaussCVrecent';
            end
            CVname = string(KernelFunctionName{j})+".mat";
            save(CVname,"WeightGaussCV10")
            %% Gauss weight
            sigmaL = GaussModel.KernelInformation.KernelParameters(1:2); % Learned length scales
            weights = exp(-sigmaL); % Predictor weights
%             weights = weights/sum(weights); % Normalized predictor weights
            
            tbl_weight = table(weights,'VariableNames',{'Predictor Weight'}, ...
                'RowNames',GaussModel.ExpandedPredictorNames);
            % tbl_weight = sortrows(tbl_weight,'Predictor Weight');
            WeightGauss(i,:) = weights';
            % plot
            % b = barh(categorical(tbl_weight.Row,tbl_weight.Row),tbl_weight.('Predictor Weight'));
            % b.Parent.TickLabelInterpreter = 'none';
            % xlabel('Predictor Weight')
            % ylabel('Predictor')
            %
            
        end
        varnames = ["epi_GDL","epi_CL"];
        KernelFunctionNamej = string(KernelFunctionName{j});
        writematrix(varnames,"Weights CL.xlsx","Sheet",KernelFunctionNamej,"Range","B1")
        writematrix(WeightGauss,"Weights CL.xlsx","Sheet",KernelFunctionNamej,"Range","B2")
        toc
    else
        WeightTreeCV10= zeros(length(phix),2,NoKFold);
        rsquare = zeros(length(phix),1);
        for i = 1:length(phix)
            index = find(abs(MFE(:,3)-phix(i))<0.001);
            
            X1 = X(index,1:2);
            Y1 = Y(index,:);
            % TreeModel = fitrtree(X,Y,'MaxNumSplits',10);
            %% Train Tree model and CV
            TreeModel = fitrtree(X1,Y1,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus','CVPartition',cvpartition(length(X1),"KFold",NoKFold)));
            imp_each = predictorImportance(TreeModel);
            
            TreeModelCV = fitrtree(X1,Y1,'CrossVal','on','KFold',NoKFold);  
            
            % TreeModelCV = fitrtree(X1,Y1,'OptimizeHyperparameters','auto',...
            %     'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
            %     'expected-improvement-plus','CVPartition',cvpartition(length(X1),"KFold",NoKFold)));
            
%             ypred = kfoldPredict(TreeModelCV);
%             SSE = sum((ypred      -       Y1).^2);
%             SST = sum((Y1 - mean(Y1)).^2);
%             rsquare(i) = 1 - SSE./SST
            
            
            
            f = @(CMP,Xtrain,Ytrain,Wtrain,Xtest,Ytest,Wtest)...
                predict(CMP,X1);            
            ypredx = kfoldfun(TreeModelCV,f);
            ypredCV = mean(ypredx)';
            SSE = sum((ypredCV      -       Y1).^2);
            SST = sum((Y1 - mean(Y1)).^2);
            rsquare(i) = 1 - SSE./SST;
            close all
            
            for k = 1:NoKFold
                WeightTreeCVrecent = predictorImportance(TreeModelCV.Trained{k,1});
                
                WeightTreeCV10(i,:,k) = WeightTreeCVrecent;
            end
            save("WeightTreeCV10.mat","WeightTreeCV10")
            save("rsquare10CV.mat","rsquare")
            %% Tree Importance
            varnames = ["epi_GDL","epi_CL"];
            WeightTree(i,:) = imp_each;
            writematrix(WeightTree,"Weights CL.xlsx","Sheet","Tree","Range","B2")
            writematrix(varnames,"Weights CL.xlsx","Sheet","Tree","Range","B1")
        end
    end
end


toc















