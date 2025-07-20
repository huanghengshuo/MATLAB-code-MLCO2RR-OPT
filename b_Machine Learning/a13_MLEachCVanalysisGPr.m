% Step 3
clc
clear all
close all
%%%%% Cross-validate machine learning model
KernelFunctionName = {'ardsquaredexponential';...
    'ardexponential';...
    'ardmatern32';...
    'ardmatern52';...
    'ardrationalquadratic'};
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


%% Junk plot
% for j = 1:length(KernelFunctionName)
% % for j = 1:1
%     figure(j)
%     hold on
%     load(string(KernelFunctionName{j})+".mat")
%     plot(phix,WeightGaussCV10(:,:,1),'-');
%     plot(phix,WeightGaussCV10(:,:,2),'-');
%     KernelFunctionNamej = KernelFunctionName{j,1};
%     WeightGauss = readmatrix("Weights CL.xlsx","Sheet",KernelFunctionNamej,"Range","B2");
%
%     plot(phix,WeightGauss(:,1),'-','LineWidth',3)
%     plot(phix,WeightGauss(:,2),'-','LineWidth',3)
%
%     title(string(KernelFunctionName{j}))
% end

% % % % % figure
% % % % % hold on
% % % % % for i = 1:10
% % % % %     plot(WeightGaussCV10(:,:,i),'-','LineWidth',2)
% % % % % end


%% plot 3D Gauss
for j = 1:length(KernelFunctionName)
% for j = 1:1
    
    load(string(KernelFunctionName{j})+".mat")
    KernelFunctionNamej = KernelFunctionName{j,1};
    WeightGauss = readmatrix("Weights CL.xlsx","Sheet",KernelFunctionNamej,"Range","B2");
    hold on
    
    c = colormap(cool(11));
    
    
    fig = figure;
    fig.Position = [1005,228,560*2,849];
    
    hold on
    t = tiledlayout(6,2);
    % t = tiledlayout("vertical")
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    for i = 1:1:10
        nexttile
        hold on
        box on
        x = [1, ones(size(phix)).*1,1]*i;
        y = [phix(1),phix,phix(end)];
        z = [0,WeightGaussCV10(:,1,i)',0];
        % fill3(x,y,z,'r','FaceColor',c(i,:),'FaceAlpha',0.2)
        fill(y,z,'r','FaceColor','k','FaceAlpha',0.2)
        
        y1 = [phix(1),phix,fliplr(phix)];
        z1 = [WeightGaussCV10(1,1,i),WeightGaussCV10(:,1,i)',fliplr(WeightGaussCV10(:,1,i)'+WeightGaussCV10(:,2,i)')];
        x1 = ones(size(z1))*i;
        fill(y1,z1,'r','FaceColor',c(i,:),'FaceAlpha',0.8)
        set(gca,'FontSize',15)
        xlim([-2.0 -0.6])
    end
    
    % WeightTree = readmatrix("Weights CL.xlsx","Sheet","Tree","Range","B2");
    % WeightTree = WeightTree./sum(WeightTree,2);
    nexttile
    box on
    hold on
    x = [1, ones(size(phix)).*1,1]*0;
    y = [phix(1),phix,phix(end)];
    z = [0,WeightGauss(:,1)',0];
    fill(y,z,'r','FaceColor','k','FaceAlpha',0.8,'LineWidth',2)
    
    y1 = [phix(1),phix,fliplr(phix)];
    z1 = [WeightGauss(1,1),WeightGauss(:,1)',fliplr(WeightGauss(:,1)'+WeightGauss(:,2)')];
    x1 = ones(size(z1))*0;
    fill(y1,z1,'r','FaceColor',c(1,:),'FaceAlpha',0.8)
    set(gca,'FontSize',15)
    xlim([-2.0 -0.6])
    % xlabel(t,"$$ \phi \ \mathrm{vs\ NHE}$$","FontSize",20,"Interpreter","latex");
    % ylabel(t,"$$ \mathrm{Feature\ Importance} $$","FontSize",20,"Interpreter","latex");
    
    xlabel(t,"Potential (V \it vs. \rm NHE)","FontSize",18,"FontName","Arial");
    ylabel(t,"Feature Importance","FontSize",18,"FontName","Arial");
    
    
    title(t,char(string(KernelFunctionName{j})),"FontSize",18,"FontName","Arial");
    legend(["\epsilon_{GDL}","\epsilon_{CL}"],'FontSize',13,"FontName","Arial",...
        "Position",[0.732987103187256 0.0800942285041224 0.0973700396698858 0.0699306658484859])
    print(['0Fig13 CVimp 2D',char(string(KernelFunctionName{j}))],'-djpeg','-r1200')
end

%% plot 3D Tree

figTree = figure;
figTree.Position = [241,168,999,810];
hold on
% hold on
% x = [1, ones(size(phix'))];
% y = [phix(1), phix'];
% z = [0, WeightGaussCV10(:,1,1)];
% fill3(x,y,z,'r')
c = colormap(cool(20));summer;
c1 = colormap(summer(10));
load("WeightTreeCV10.mat");

for i = 1:10
    
    x = [1, ones(size(phix)).*1,1]*i;
    y = [phix(1),phix,phix(end)];
    z = [0,WeightTreeCV10(:,1,i)',0];
    % fill3(x,y,z,'r','FaceColor',c(i,:),'FaceAlpha',0.2)
    fill3(x,y,z,'r','FaceColor','k','FaceAlpha',0.2)
    
    y1 = [phix(1),phix,fliplr(phix)];
    z1 = [WeightTreeCV10(1,1,i),WeightTreeCV10(:,1,i)',fliplr(WeightTreeCV10(:,1,i)'+WeightTreeCV10(:,2,i)')];
    x1 = ones(size(z1))*i;
    fill3(x1,y1,z1,'r','FaceColor',c(i,:),'FaceAlpha',0.8)
end

WeightTree = readmatrix("Weights CL.xlsx","Sheet","Tree","Range","B2");

x = [1, ones(size(phix)).*1,1]*0;
y = [phix(1),phix,phix(end)];
z = [0,WeightTree(:,1)',0];
fill3(x,y,z,'r','FaceColor','k','FaceAlpha',0.8,'LineWidth',2)

y1 = [phix(1),phix,fliplr(phix)];
z1 = [WeightTree(1,1),WeightTree(:,1)',fliplr(WeightTree(:,1)'+WeightTree(:,2)')];
x1 = ones(size(z1))*0;
fill3(x1,y1,z1,'r','FaceColor',c(1,:),'FaceAlpha',0.8)

view([-1.265975672200814e+02,58.077646538927603])
title("$$ \mathrm{Decision\ Tree} $$","FontSize",15,"Interpreter","latex");
xlabel("$$ \mathrm{Cross-Validation} $$","FontSize",15,"Interpreter","latex","Position",[5.866803796140631,-0.406437753674645,-6.3982719e-7],"Rotation",-39);
ylabel("$$ \phi \ \mathrm{vs\ NHE}$$","FontSize",15,"Interpreter","latex","Position",[-1.135933858773001,-1.343703293273506,-0.000003564697478],"Rotation",24);
zlabel("$$ \mathrm{Feature\ Importance} $$","FontSize",15,"Interpreter","latex");
grid on
grid minor
legend(["$$\varepsilon_\mathrm{GDL} $$","$$\varepsilon_\mathrm{CL} $$"],"Interpreter","latex",'FontSize',15,...
    "Position",[0.743344993380644 0.829320989420385 0.0918487750506336 0.0648148130487513])
print('0Fig13 CVimp Tree Accumulation 3D','-djpeg','-r1200')




%%  Accumulation
load("WeightTreeCV10.mat");
% WeightTreeCV10 = WeightTreeCV10./sum(WeightTreeCV10,2);
fig = figure;
fig.Position = [263,127,1120,1509];

hold on
t = tiledlayout(11,2);
% t = tiledlayout("vertical")
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:1:20
    nexttile
    hold on
    box on
    x = [1, ones(size(phix)).*1,1]*i;
    y = [phix(1),phix,phix(end)];
    z = [0,WeightTreeCV10(:,1,i)',0];
    % fill3(x,y,z,'r','FaceColor',c(i,:),'FaceAlpha',0.2)
    fill(y,z,'r','FaceColor','k','FaceAlpha',0.2)
    
    y1 = [phix(1),phix,fliplr(phix)];
    z1 = [WeightTreeCV10(1,1,i),WeightTreeCV10(:,1,i)',fliplr(WeightTreeCV10(:,1,i)'+WeightTreeCV10(:,2,i)')];
    x1 = ones(size(z1))*i;
    fill(y1,z1,'r','FaceColor',c(i,:),'FaceAlpha',0.8)
    set(gca,'FontSize',15)
    xlim([-2.0 -0.6])
end

WeightTree = readmatrix("Weights CL.xlsx","Sheet","Tree","Range","B2");
% WeightTree = WeightTree./sum(WeightTree,2);
nexttile
box on
hold on
x = [1, ones(size(phix)).*1,1]*0;
y = [phix(1),phix,phix(end)];
z = [0,WeightTree(:,1)',0];
fill(y,z,'r','FaceColor','k','FaceAlpha',0.8,'LineWidth',2)

y1 = [phix(1),phix,fliplr(phix)];
z1 = [WeightTree(1,1),WeightTree(:,1)',fliplr(WeightTree(:,1)'+WeightTree(:,2)')];
x1 = ones(size(z1))*0;
fill(y1,z1,'r','FaceColor',c(1,:),'FaceAlpha',0.8)
set(gca,'FontSize',15)
xlim([-2.0 -0.6])
% xlabel(t,"$$ \phi \ \mathrm{vs\ NHE}$$","FontSize",20,"Interpreter","latex");
% ylabel(t,"$$ \mathrm{Feature\ Importance} $$","FontSize",20,"Interpreter","latex");

xlabel(t,"Potential (V \it vs. \rm NHE)","FontSize",18,"FontName","Arial");
ylabel(t,"Feature Importance","FontSize",18,"FontName","Arial");

% legend(["\epsilon_{GDL}","\epsilon_{CL}"],'FontSize',13,"FontName","Arial")

legend(["\epsilon_{GDL}","\epsilon_{CL}"],'FontSize',15,"FontName","Arial",...
    'Position',[0.678511905746445 0.0516185448149628 0.0929166656821259 0.0515727856862083])
print('0Fig13 CVimp Tree Accumulation ALL 02','-djpeg','-r1200')









