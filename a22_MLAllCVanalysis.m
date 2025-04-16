%% Step 2
clc
clear all
close all

%%
BarColors1 = [163 210 226; 134 192 203; 93 163 175; 16 137 148]./255;
BarColors2 = [168 191 145; 105 148  95; 77 128 72; 69 107 56]./255;

% BarColors1 = [163 210 226; 134 192 203; 16 137 148]./255;
% BarColors2 = [168 191 145; 105 148  95; 69 107 56]./255;

BarColors1 = [145 197 218; 71 146 196; 36 100 171]./255;
BarColors2 = [168 191 145; 105 148  95; 69 107 56]./255;


KernelFunctionName = {'ardsquaredexponential';...
    'ardexponential';...
    'ardmatern32';...
    'ardmatern52';...
    'ardrationalquadratic'};


%% plot Gauss
% for j = 1:length(KernelFunctionName)
for j = 1:1

WeightGauss = readmatrix("Weights ALL.xlsx","Sheet",KernelFunctionName{j},"Range","B2");
xset = 1:11;
figure

b1 = bar(xset,WeightGauss,'stacked');
for i = 1:3
b1(i).FaceColor = BarColors1(i,:);
b1(i).EdgeColor = BarColors1(i,:);
end

legend({'epiGDL';'epiCL';'phi'})
ylim([0 1])
end

for j = 1:5
WeightGauss = readmatrix("Weights ALL.xlsx","Sheet",KernelFunctionName{j},"Range","B2");

yNo = length(WeightGauss);
x = 1:yNo-1;
fig = figure;
hold on
% fig.Position = [389,81,929,785];
b2 = bar(x,WeightGauss(x,:),'stacked');
for i = 1:3
b2(i).FaceColor = BarColors1(i,:);
b2(i).EdgeColor = BarColors1(i,:);
% b2(i).BarWidth = 0.8;
end

b2 = bar(yNo,WeightGauss(end,:),'stacked');
for i = 1:3
b2(i).FaceColor = BarColors2(i,:);
b2(i).EdgeColor = BarColors2(i,:);
end

WeightGaussmean= mean(WeightGauss);
WeightGaussstd = std(WeightGauss,0,1);
% figure
hold on
b3 = bar(0,WeightGaussmean,'stacked');
for i = 1:3
b3(i).FaceColor = BarColors1(i,:);
b3(i).EdgeColor = BarColors1(i,:);
end

c = cumsum(WeightGaussmean);
for i = 1:3
  er = errorbar(0,c(i),WeightGaussstd(i),WeightGaussstd(i));
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.CapSize = 10;  
end
ylim([0 1])
set(gca,'FontSize',18)
title(string([char(string(KernelFunctionName{j}))]),"FontSize",18,"FontName","Arial");
xlabel("Cross-Validation","FontSize",18,"FontName","Arial");
ylabel("Feature importance","FontSize",18,"FontName","Arial");
legend(["\epsilon_{GDL}","\epsilon_{CL}","\phi"],...
    'FontSize',15,...
    "FontName","Arial");
print(['0Fig22 CVimpALL ',char(string(KernelFunctionName{j}))],'-djpeg','-r1200')
end

%% plot Tree
WeightTree = readmatrix("Weights ALL.xlsx","Sheet","Tree","Range","B2");
yNo = length(WeightTree);
x = 1:yNo-1;
fig = figure;
set(gca,'FontSize',18)
hold on
fig.Position = [185,316,863,424];
b2 = bar(x,WeightTree(x,:),'stacked');
for i = 1:3
b2(i).FaceColor = BarColors1(i,:);
b2(i).EdgeColor = BarColors1(i,:);
b2(i).BarWidth = 0.6;
end

b2 = bar(yNo,WeightTree(end,:),'stacked');
for i = 1:3
b2(i).FaceColor = BarColors2(i,:);
b2(i).EdgeColor = BarColors2(i,:);
b2(i).BarWidth = 0.6;
end

WeightTreemean= mean(WeightTree);
WeightTreestd = std(WeightTree,0,1);
% figure
hold on
b3 = bar(0,WeightTreemean,'stacked');
for i = 1:3
b3(i).FaceColor = BarColors1(i,:);
b3(i).EdgeColor = BarColors1(i,:);
b3(i).BarWidth = 0.6;
end

c = cumsum(WeightTreemean);
for i = 1:3
  er = errorbar(0,c(i),WeightTreestd(i),WeightTreestd(i));
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.CapSize = 8;  
end

title("Decision Tree","FontSize",18,"FontName","Arial");
xlabel("Cross-Validation","FontSize",18,"FontName","Arial");
ylabel("Feature importance","FontSize",18,"FontName","Arial");
legend(["\epsilon_{GDL}","\epsilon_{CL}","\phi"],...
    'FontSize',15,...
    "FontName","Arial","Position",[0.775589031964115,0.573899377324858,0.107763613845078,0.221698106951871]);
print('0Fig22 CVimpALL Tree','-djpeg','-r1200')


%%
WeightTreeFE = readmatrix("Weights ALL.xlsx","Sheet","Tree FE","Range","B2");
WeightTreeFE = WeightTreeFE./sum(WeightTreeFE,2);
WeightTreemeanFE= mean(WeightTreeFE);
WeightTreestdFE = std(WeightTree,0,1);
fig = figure;
set(fig,'Position',[993,447,560,423]);
hold on
b3 = bar(0,WeightTreemeanFE,'stacked');
for i = 1:3
b3(i).FaceColor = BarColors1(i,:);
b3(i).EdgeColor = BarColors1(i,:);
b3(i).BarWidth = 0.5;
end

c = cumsum(WeightTreemeanFE);
% for i = 1:3
% er = errorbar(0,c(i),WeightTreestdFE(i),WeightTreestdFE(i));
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% er.CapSize = 8;  
% end

%
WeightTree = readmatrix("Weights ALL.xlsx","Sheet","Tree","Range","B2");
WeightTree = WeightTree./sum(WeightTree,2);
WeightTreemean= mean(WeightTree);
WeightTreestd = std(WeightTree,0,1);
% figure
hold on
box on
b3 = bar(1.4,WeightTreemean,'stacked');
for i = 1:3
b3(i).FaceColor = BarColors1(i,:);
b3(i).EdgeColor = BarColors1(i,:);
b3(i).BarWidth = 0.5;
end

c = cumsum(WeightTreemean);
% for i = 1:3
%   er = errorbar(1,c(i),WeightTreestd(i),WeightTreestd(i));
% er.Color = [1 0 0];                            
% er.LineStyle = 'none';  
% er.CapSize = 8;  
% end
xlim([-1.0 2.5])
ylim([0 1.00])
set(gca,'FontSize',18,'XTick',[0 1.6])
set(gca,'XTickLabel',{'ML(FE)','ML(res FE)'}, 'XTickLabelRotation',15)


ylabel("Relative feature importance","FontSize",18,"FontName","Arial");
legend(["\epsilon_{GDL}","\epsilon_{CL}","\phi"],...
    'FontSize',12,'Position',[0.73056522954137 0.687043591742233 0.163851654063542 0.183333328224364],'EdgeColor','none',...
    'Color','none',"FontName","Arial");
print('0Fig22 CVimpALL FE vs resiFE','-djpeg','-r1200')

















