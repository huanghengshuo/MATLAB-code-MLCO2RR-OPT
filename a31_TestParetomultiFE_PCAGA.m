clc
clear all
close all

load Pareto3varsGA.mat
load c.mat
[wcoeff1,score1,latent1,tsquared1,explained1,mu1] = pca(xx,'VariableWeights','variance');
% Fval(:,3) = abs(Fval(:,3));

figure
plot(score1(:,1),score1(:,2),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

figure
pareto(explained1)
xlabel('Principal Component')
ylabel('Variance Explained (%)')

%%
X = XX;
% tt = score1*wcoeff1' + repmat(mu1,size(xx,1),1);
Sxx = (X - repmat(mu1,size(X,1),1))/wcoeff1';
% S = (xx - repmat(mu1,size(xx,1),1))/wcoeff1';
[xq,yq] = meshgrid(linspace(min(Sxx(:,1)),max(Sxx(:,1)),100),...
    linspace(min(Sxx(:,2)),max(Sxx(:,2)),100));
vq = griddata(Sxx(:,1),Sxx(:,2),Fval(:,1),xq,yq);
figure
hold on
scatter3(score1(:,1),score1(:,2),fval(:,1),'o','SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','c',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7)
% scatter3(Sxx(:,1),Sxx(:,2),-Fval(:,1),'.')
s = surf(xq,yq,vq,'FaceAlpha',0.8);
s.FaceColor = 'interp';
s.EdgeColor = 'none';
[~,con] = contour(xq,yq,vq,3);
% con.ShowText = 'on';
con.LineWidth = 2.5;
con.ZLocation = 0.85;
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('obj FE')
view([-25 28])
grid minor
colormap(summer)
colorbar
zlim([0.85 1])
%%
vq = griddata(Sxx(:,1),Sxx(:,2),Fval(:,2),xq,yq);
figure
hold on
scatter3(score1(:,1),score1(:,2),fval(:,2),'o','SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','c',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7)
s = surf(xq,yq,vq,'FaceAlpha',0.8);
s.FaceColor = 'interp';
s.EdgeColor = 'none';
[~,con] = contour(xq,yq,vq,3);
% con.ShowText = 'on';
con.LineWidth = 2.5;
con.ZLocation = 1.4e-15;
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('obj EI')
view([-25 28])
grid minor
colormap(spring)
colorbar
%%
vq = griddata(Sxx(:,1),Sxx(:,2),Fval(:,3),xq,yq);
figure
hold on

scatter3(score1(:,1),score1(:,2),fval(:,3),'o','SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','c',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7)
s = surf(xq,yq,vq,'FaceAlpha',0.8);
s.FaceColor = 'interp';
s.EdgeColor = 'none';
[~,con] = contour(xq,yq,vq,3);

% scatter3(score1(:,1),score1(:,2),log10(-fval(:,3)),'o','SizeData',50,...
%     'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','c',...
%     'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7)
% s = surf(xq,yq,log10(vq),'FaceAlpha',0.8);
% s.FaceColor = 'interp';
% % s.EdgeColor = 'none';
% [~,con] = contour(xq,yq,log10(vq),3);

% con.ShowText = 'on';
con.LineWidth = 2.5;
con.ZLocation = -1.5;
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('obj j_{COER}')
% zlabel('log10(j_{COER})')
view([-25 28])
grid minor
colormap(c1)
colorbar


%% PCA for fval
% [wcoeff,score,latent,tsquared,explained,mu] = pca(fval,'VariableWeights','variance');
% 
% figure
% plot(score(:,1),score(:,2),'+')
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% 
% figure
% pareto(explained)
% xlabel('Principal Component')
% ylabel('Variance Explained (%)')

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
phix = unique(MFE(:,3))';
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

tt = score1*wcoeff1' + repmat(mu1,size(xx,1),1);




