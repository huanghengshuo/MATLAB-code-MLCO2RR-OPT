clc
clear all
close all
global  TreeModel TreeModeljCO  GaussModel
global  TreeModelCV TreeModelCVjCO GaussModeljCO
rng('default') % For reproducibility
tic
%% Input original data

% -- 5 points
% M0 = readmatrix("jCOER multi test 01.xlsx"); % GDL CL  phi
% M1 = readmatrix("jHER multi test 01.xlsx");
% M2 = readmatrix("jTotal multi test 01.xlsx");

M0 = readmatrix("jCOER multi test.xlsx"); % GDL CL  phi
M1 = readmatrix("jHER multi test.xlsx");
M2 = readmatrix("jTotal multi test.xlsx");

% -- 3 points
% M0 = readmatrix("jCOER multi test 03.xlsx"); % GDL CL  phi
% M1 = readmatrix("jHER multi test 03.xlsx");
% M2 = readmatrix("jTotal multi test 03.xlsx");
MCOER = rmmissing(M0);
MHER  = rmmissing(M1);
MTotal= rmmissing(M2);

MFE = MCOER;
% MFE(:,end) = MCOER(:,end)./(MCOER(:,end)+MTotal(:,end));
MFE(:,end) = MCOER(:,end)./MTotal(:,end);
MCOER(:,end) = -MCOER(:,end);

index = find(MFE(:,1)<=MFE(:,2));
MFE(index,:)   = [];
MCOER(index,:) = [];

index = find(MFE(:,end-1)<=-2.001);
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


% j0   = unique(MFE(:,3))';
% alphac=unique(MFE(:,4))';
phix = unique(MFE(:,3))';
%% plot ori

% plot(MFE(:,3),MFE(:,4),'.k')


%% Prepare the data for Machine Learning
X = MFE(:,1:end-1);
Y = MFE(:,end);
ystd = zeros(1,length(phix));
YmeanInd = zeros(size(phix));
Yminusmean = zeros(size(Y));

for i = 1:1:length(phix)
    index = find(abs(MFE(:,end-1)-phix(i))<0.001);
    
    X1 = X(index,1:end-1);
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

plot(X(:,end),Y,'k.')
plot(phix,YmeanInd,'c-','LineWidth',2)


figure
hold on
plot(X(:,end),Yminusmean,'k.')
plot(phix,YmeanInd-YmeanInd,'c-','LineWidth',2)

%% Prepare the data for Machine Learning for jCOER
X = MFE(:,1:end-1);
YjCOER = MCOER(:,end);
ystdjCOER = zeros(1,length(phix));
YmeanIndjCOER = zeros(size(phix));
YminusmeanjCOER = zeros(size(Y));

for i = 1:1:length(phix)
    index = find(abs(MFE(:,end-1)-phix(i))<0.001);
    
    X1 = X(index,1:end);
    Y1jCOER = YjCOER(index,:);
    ystdjCOER(1,i) = std(Y1jCOER);
    YmeanIndjCOER(1,i) = mean(Y1jCOER);
    YminusmeanjCOER(index,:) = Y1jCOER - mean(Y1jCOER);
end
yyaxis right
set(gca,'YColor','k')
plot(phix,ystdjCOER,'k-','LineWidth',2)
ylim([-max(ystdjCOER)*1.0 max(ystdjCOER)*1.05])


% figure
% hold on
% 
% plot(X(:,3),YjCOER,'k.')
% plot(phix,YmeanIndjCOER,'c-','LineWidth',2)
% 
% 
% figure
% hold on
% plot(X(:,3),YminusmeanjCOER,'k.')
% plot(phix,YmeanIndjCOER-YmeanIndjCOER,'c-','LineWidth',2)
%% X-x1 x2 phi  Y-dFE
X = MFE(:,1:end-1);

Y = Yminusmean;
YCO = YminusmeanjCOER;
% 
% X = Nnew;
% Y = Mnew(:,end);

NoKFold = 20;

TreeModel = fitrtree(X,Y);
TreeModelCV = fitrtree(X,Y,'CrossVal','on','KFold',NoKFold);

TreeModeljCO = fitrtree(X,YCO);
TreeModelCVjCO = fitrtree(X,YCO,'CrossVal','on','KFold',NoKFold);

GaussModel = fitrgp(phix',YmeanInd');
GaussModeljCO = fitrgp(phix',YmeanIndjCOER');
figure
hold on
plot(phix,YmeanInd,'*');
plot(phix,resubPredict(GaussModel),'*');

figure
hold on
plot(phix,YmeanIndjCOER,'*');
plot(phix,resubPredict(GaussModeljCO),'*');
% GaussModel= fitrgp(X,Y,'CategoricalPredictors',[1 2], ...
%     'KernelFunction','ardsquaredexponential');

% GaussModel= fitrgp(X,Y,'CategoricalPredictors',[1 2]);
% GaussModel= fitrgp(X,Y);
% GaussModeljCO = fitrgp(X,YCO);




% view(TreeModel);
% imp = predictorImportance(TreeModelCV.Trained{1,1});
% figure
% bar(imp)

% figure
% view(TreeModel,'Mode','graph')

% fig = figure;
% fig.Position = [86,78,1808,415];
% subplot(2,3,1)
% plotPartialDependence(TreeModel,1)
% subplot(2,3,2)
% plotPartialDependence(TreeModel,2)
% subplot(2,3,3)
% plotPartialDependence(TreeModel,3)
% subplot(2,3,4)
% plotPartialDependence(TreeModel,4)
% subplot(2,3,5)
% plotPartialDependence(TreeModel,5)
% figure
% subplot(2,2,1)
% plotPartialDependence(GaussModel,1)
% subplot(2,2,2)
% plotPartialDependence(GaussModel,2)
% subplot(2,2,3)
% plotPartialDependence(GaussModel,3)
%% figure
% figure
% plot(X(:,end),Y,'.');
% 
% Ypred = kfoldPredict(TreeModelCV);
% hold on
% plot(X(:,end),Ypred,'.k');

%% figure
f = @(CMP,Xtrain,Ytrain,Wtrain,Xtest,Ytest,Wtest)...
    predict(CMP,X);

ypredx = kfoldfun(TreeModelCV,f);
ypred = mean(ypredx)';
% figure
% hold on
% plot(Ypred,ypred,'*');
% plot([min(Ypred),max(Ypred)],[min(Ypred),max(Ypred)],'-');
% xlabel('kfoldPredict')
% ylabel('meanPredictkfold')
% 
% figure
% hold on
% plot(Y,ypred,'*');
% plot([min(Ypred),max(Ypred)],[min(Ypred),max(Ypred)],'-');
% xlabel('true y')
% ylabel('meanPredictkfold')


%% Gauss weight
% sigmaL = GaussModel.KernelInformation.KernelParameters(1:end-1); % Learned length scales
% weights = exp(-sigmaL); % Predictor weights
% weights = weights/sum(weights); % Normalized predictor weights
% 
% tbl_weight = table(weights,'VariableNames',{'Predictor Weight'}, ...
%     'RowNames',GaussModel.ExpandedPredictorNames);
% tbl_weight = sortrows(tbl_weight,'Predictor Weight');
% figure
% b = barh(categorical(tbl_weight.Row,tbl_weight.Row),tbl_weight.('Predictor Weight'));
% b.Parent.TickLabelInterpreter = 'none'; 
% xlabel('Predictor Weight')
% ylabel('Predictor')


%% plotPartialDependence
fig = figure;
fig.Position = [81,51,1808,864];
for i = 1:3
    subplot(1,3,i)
    plotPartialDependence(TreeModel,i,"Conditional","absolute")
end

%% plotPartialDependence
fig = figure;
fig.Position = [81,51,1808,864];
for i = 1:3
    subplot(1,3,i)
    plotPartialDependence(TreeModeljCO,i,"Conditional","absolute")
end
%% plotPartialDependence
% figure
% subplot(2,2,1)
% plotPartialDependence(GaussModel,1,"Conditional","absolute")
% subplot(2,2,2)
% plotPartialDependence(GaussModel,2,"Conditional","absolute")
% subplot(2,2,3)
% plotPartialDependence(GaussModel,3,"Conditional","absolute")
% subplot(2,2,4)
% plotPartialDependence(GaussModel,4,"Conditional","absolute")

%% Compare with Experimental data

% phi = -0.2:-0.0002:-2.4; % vs NHE
% % phi = -0.2:-0.02:-2.4; % vs NHE
% % phi = phi + 0.05916*pH;     % vs RHE
% 
% pH = 7.23;
% % phi_ref = [-1.138 -1.02 -0.901 -0.756 -0.615 -0.479]; % vs RHE
% % % phi_ref = phi_ref - 0.05916*pH; % vs NHE
% % jCO = [69.305 45.997 28.481 12.989 4.546 0.162]; % mA/cm^2
% % jH2 = [1.902 1.455 1.565 0.948 0.385 0.159]; % mA/cm^2
% 
% phi_ref = [-1.15 -1.03 -0.903 -0.764 -0.624 -0.474]; % vs RHE
% % phi_ref = phi_ref - 0.05916*pH; % vs NHE
% jCO = [70.1 47.3 29.3 13.6 3.8 0.01]; % mA/cm^2
% jH2 = [2.0 1.5 1.6 1.0 0.5 0.3]; % mA/cm^2
% 
% 
% %
% jCO = (jCO-min(jCO))./min(jCO);
% 
% 
% 
% jTotal = jCO + jH2;
% FECO = jCO./jTotal;
% 
% 
% Xpred = zeros(length(phi),3);
% Xpred(:,1) = 0.7;
% Xpred(:,2) = 0.3;
% Xpred(:,3) = phi;
% Xpred(:,4) = 1.8e-5;
% % Xpred(:,5) = phi;
% YpredTree = predict(TreeModel,Xpred);
% YpredGauss= predict(GaussModel,Xpred);
% index = find((abs(MFE(:,1)-0.7)<0.001) & (abs(MFE(:,2)-0.3)<0.001));
% figure
% hold on
% % plot(MFE(index,3)+0.05916*pH,Mnew(index,1),'*',"MarkerSize",12)
% 
% % plot(phi_ref,jCO,"s","MarkerSize",12,"Color",[0.9290 0.6940 0.1250],"MarkerFaceColor",[0.9290 0.6940 0.1250])
% 
% plot(phi+0.05916*pH,YpredTree,"-","LineWidth",3)
% plot(phi+0.05916*pH,YpredGauss,"r--","LineWidth",3)
% % legend({'COMSOL';'exp';'Tree';'Gauss'})

%% Pareto
figure
xGDL = unique(MFE(:,1));
xCL = unique(MFE(:,2));
% j0   = unique(MFE(:,3))';
% alphac=unique(MFE(:,4))';

lb = [min(xGDL);min(xCL);-1.7];
ub = [max(xGDL);max(xCL);-1.4];
% options = optimoptions('paretosearch','ParetoSetSize',100,...
%     'PlotFcn',{'psplotparetof' 'psplotparetox'});

%% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',300);

% [xx,fval,exitflag,output,residuals] = paretosearch(@HHtarget,3,[],[],[],[],lb,ub,@HHcirclecons,options);
[xx,fval,exitflag,output,residuals] = paretosearch(@HHtarget,3,[-1,1,0],0,[],[],lb,ub,[],options);

%% gamultiobj

% options = optimoptions("gamultiobj","PopulationSize",500);
% [xx,fval,exitflag,output,population,scores] = gamultiobj(@HHtarget,3,[-1,1,0],0,[],[],lb,ub,[],options);

% index = find(xx(:,1)<xx(:,2))
% xx(index,:) = [];
% fval(index,:) = [];


%%
[Xori0,Yori0,Zori0] = meshgrid(linspace(min(xGDL),max(xGDL),20),...
    linspace(min(xCL),max(xCL),20),...
    linspace(-1.7,-1.4,300));
Xori1 = reshape(Xori0,[],1);
Yori1 = reshape(Yori0,[],1);
Zori1 = reshape(Zori0,[],1);
% XX = [Xori1';Yori1';Zori1'];
XX = cat(2,Xori1,Yori1,Zori1);
index = find(XX(:,1)<XX(:,2));
XX(index,:) = [];
% Xori0 = meshgrid(linspace(min(xGDL),max(xGDL),20),linspace(min(xCL),max(xCL),20));
% Xori1 = reshape(Xori0,[],2);

% XX = X;

Fval = HHtargetori(XX);
Fval(:,1) = -Fval(:,1);
Fval(:,3) = -Fval(:,3);

fval(:,1) = -fval(:,1);
fval(:,3) = -fval(:,3);


% save xx and fval
save Pareto3varsPS.mat xx fval XX Fval

% X = XX;

figure
hold on
plot(X(:,end),Y,'.')
plot(xx(:,end),predict(TreeModel,xx),'r*')
% plot(xx(:,3),kfoldpredict(GaussModel,xx),'r*')



figure
hold on
plot(X(:,end),MFE(:,end),'o')
plot(xx(:,end),predict(TreeModel,xx)+predict(GaussModel,xx(:,end)),'r*')
plot(XX(:,end),Fval(:,1),'.')

figure
hold on
plot(X(:,end),MFE(:,end),'o')
plot(xx(:,end),fval(:,1),'r*')
plot(XX(:,end),Fval(:,1),'.')

%% Pareto Search Plot All
f = @(CMP,Xtrain,Ytrain,Wtrain,Xtest,Ytest,Wtest)...
    predict(CMP,xx);

figure
% -------- diag --------begin
subplot(3,3,1)
hold on
plot(X(:,1),MFE(:,end),'.');
plot(xx(:,1),mean(kfoldfun(TreeModelCV,f))'+predict(GaussModel,xx(:,end)),'r*');
xlabel("$$ x_{\mathrm{GDL}} $$","Interpreter","latex")
ylabel("$$ \mathrm{FE} $$","Interpreter","latex")

subplot(3,3,5)
hold on
plot(X(:,2),MFE(:,end),'.');
plot(xx(:,2),mean(kfoldfun(TreeModelCV,f))'+predict(GaussModel,xx(:,end)),'r*');
xlabel("$$ x_{\mathrm{CL}} $$","Interpreter","latex")
ylabel("$$ \mathrm{FE} $$","Interpreter","latex")

subplot(3,3,9)
hold on
plot(X(:,3),MFE(:,end),'.');
plot(xx(:,3),mean(kfoldfun(TreeModelCV,f))'+predict(GaussModel,xx(:,end)),'r*');
xlabel("$$ \phi $$","Interpreter","latex")
ylabel("$$ \mathrm{FE} $$","Interpreter","latex")
% -------- diag --------end


% -------- upper--------begin
subplot(3,3,2)
hold on
plot(X(:,1),X(:,2),'.');
plot(xx(:,1),xx(:,2),'r*');
xlabel("$$ x_{\mathrm{GDL}} $$","Interpreter","latex")
ylabel("$$ x_{\mathrm{CL}} $$","Interpreter","latex")

subplot(3,3,3)
hold on
plot(X(:,1),X(:,end),'.');
plot(xx(:,1),xx(:,end),'r*');
xlabel("$$ x_{\mathrm{GDL}} $$","Interpreter","latex")
ylabel("$$ \phi $$","Interpreter","latex")

subplot(3,3,6)
hold on
plot(X(:,2),X(:,end),'.');
plot(xx(:,2),xx(:,end),'r*');
xlabel("$$ x_{\mathrm{CL}} $$","Interpreter","latex")
ylabel("$$ \phi $$","Interpreter","latex")
% -------- upper--------end


% -------- lower--------begin
subplot(3,3,4)
hold on
plot(fval(:,1),fval(:,2),'cs');
xlabel("$$ obj_\mathrm{FE} $$","Interpreter","latex")
ylabel("$$ obj_\mathrm{EI} $$","Interpreter","latex")

subplot(3,3,7)
hold on
plot(fval(:,1),fval(:,3),'cs');
xlabel("$$ obj_\mathrm{FE} $$","Interpreter","latex")
ylabel("$$ obj_\mathrm{jCOER} $$","Interpreter","latex")

subplot(3,3,8)
hold on
plot(fval(:,2),fval(:,3),'cs');
xlabel("$$ obj_\mathrm{EI} $$","Interpreter","latex")
ylabel("$$ obj_\mathrm{jCOER} $$","Interpreter","latex")
% -------- lower--------end

%%

strnamex = ["x_{\mathrm{GDL}}";"x_{\mathrm{CL}}";"\phi"];
strnamex = '$$ ' + strnamex + ' $$';

strnamey = ["obj_\mathrm{FE}";"obj_\mathrm{EI}";"obj_\mathrm{jCOER}"];
strnamey = '$$ ' + strnamey + ' $$';
figPF = figure;
figPF.Position = [431.4,128.2,1143.2,852.8];
tiledlayout(3,3)
for i = 1:3
    for j = 1:3
        nexttile
        hold on
        scatter(XX(:,i),Fval(:,j),'o','SizeData',10,...
            'LineWidth',0.1,'MarkerEdgeColor','k','MarkerFaceColor','c',...
            'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.7)
        scatter(xx(:,i),fval(:,j),'o','SizeData',20,...
            'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','m',...
            'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7)
        xlabel(strnamex(i,1),"Interpreter","latex")
        ylabel(strnamey(j,1),"Interpreter","latex")
    end
end
print('0Fig31 Pareto set vs objective PS','-djpeg','-r1200') 




%% Pareto Search Plot All

% f = @(CMP,Xtrain,Ytrain,Wtrain,Xtest,Ytest,Wtest)...
%     predict(CMP,xx);
% strname = ["x_{\mathrm{GDL}}";"x_{\mathrm{CL}}";"j_{0}";"\alpha_c";"\phi"];
% strname = '$$ ' + strname + ' $$';
% 
% figure
% t = tiledlayout(5,5);
% for i = 1:5
%     for j = 1:5
%       if i == j
%           nexttile
%           plot(X(:,i),MFE(:,end),'k.');
%           hold on
%           plot(xx(:,i),mean(kfoldfun(TreeModelCV,f))'+predict(GaussModel,xx(:,end)),'r*');
%           xlabel(strname(i,1),"Interpreter","latex")
%           ylabel("$$ FE $$","Interpreter","latex")
%       elseif i > j
%           nexttile
%       else
%           nexttile
%           plot(X(:,i),X(:,j),'k.');
%           hold on
%           plot(xx(:,i),xx(:,j),'r*');
%           xlabel(strname(i,1),"Interpreter","latex")
%           ylabel(strname(j,1),"Interpreter","latex")
%       end
%     end
% end
          








toc






% [trainedModel, validationRMSE] = trainRegressionModel(MFE)
% T = MFE(:,1:3);
% yfit = trainedModel.predictFcn(T)
% figure
% plotPartialDependence(trainedModel,3)
% figure
% err = abs(yfit-MFE(:,4));
% plot(err)
% 
% figure
% hold on
% plot(yfit,'-s')
% plot(MFE(:,4),'*')
% 
% 
% 
% imp = predictorImportance(trainedModel)
% 

function y = HHtarget(x)
global  GaussModel
global TreeModelCV TreeModelCVjCO GaussModeljCO
% global TreeModel TreeModeljCO
f = @(CMP,Xtrain,Ytrain,Wtrain,Xtest,Ytest,Wtest)...
    predict(CMP,x);

% ypredx = kfoldfun(TreeModelCV,f);
% ypred = mean(ypredx)';


FE = mean(kfoldfun(TreeModelCV,f))';

FE = FE + predict(GaussModel,x(end));
jCOER = mean(kfoldfun(TreeModelCVjCO,f))'+predict(GaussModeljCO,x(end));

% FE = TreeModel(x);
% FE = predict(TreeModel,x)+predict(GaussModel,x(3));
% jCOER = predict(TreeModeljCO,x);



y(1) = -FE;
%%%%%%
% e0 = 1.6021766e-19;
Prate = 1.*1e+6/(28.01*24*3600);
Itotal = Prate.*2.*96485./FE;
PowerConsumed = Itotal.*abs(x(end));
Eelectrolyser = PowerConsumed./1000.*24.*0.0036;
y(2) = Eelectrolyser;
%%%%%%
% y(2) = -x(3)./FE;
y(3) = -jCOER;

% y = [-FE;Eelectrolyser;-jCOER];
end

function y = HHtargetori(x)
global GaussModel
global TreeModelCV TreeModelCVjCO GaussModeljCO
% global TreeModel TreeModeljCO
f = @(CMP,Xtrain,Ytrain,Wtrain,Xtest,Ytest,Wtest)...
    predict(CMP,x);

% ypredx = kfoldfun(TreeModelCV,f);
% ypred = mean(ypredx)';
y = zeros(size(x,1),3);

FE = mean(kfoldfun(TreeModelCV,f))'+predict(GaussModel,x(:,end));
jCOER = mean(kfoldfun(TreeModelCVjCO,f))'+predict(GaussModeljCO,x(:,end));

% FE = TreeModel(x);
% FE = predict(TreeModel,x)+predict(GaussModel,x(3));
% jCOER = predict(TreeModeljCO,x);


y(:,1) = -FE;

%%%%%% kW/ton CO2
% e0 = 1.6021766e-19;
Prate = 1.*1e+6/(28.01*24*3600);
Itotal = Prate.*2.*96485./FE;
PowerConsumed = Itotal.*abs(x(:,end));
Eelectrolyser = PowerConsumed./1000.*24.*0.0036;
y(:,2) = Eelectrolyser;
%%%%%%
% y(2) = -x(3)./FE;
y(:,3) = -jCOER;
end

function [c,ceq] = HHcirclecons(x)
ceq = [];
% ceq = (x(1)-0).^2+(x(2)-0).^2+(x(3)-0).^2-2^2;

% c(1) = x(end)+1.2;
% c(2) = -x(end)-1.8;
c = [];
end

