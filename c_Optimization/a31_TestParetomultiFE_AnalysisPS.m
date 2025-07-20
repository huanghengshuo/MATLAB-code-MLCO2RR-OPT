clc
clear all
close all
load Pareto3varsPS.mat
load c.mat

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

% XX = MFE(:,1:3);


% scatter3(xx(:,1),xx(:,2),xx(:,3),'rp')
% scatter3(XX(:,1),XX(:,2),XX(:,3),'k.')

% scatter3(XX(:,1),XX(:,2),XX(:,3),'o','SizeData',20,...
%     'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
%     'MarkerEdgeAlpha',0.05,'MarkerFaceAlpha',0.7)


s = scatter3(xx(:,1),xx(:,2),xx(:,3),'o','SizeData',50,...
    'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','c',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7);
view([-39 20])
xlabel('xGDL')
ylabel('xCL')
zlabel('\phi')
%%
fig = figure
hold on

% scatter3(xx(:,1),xx(:,2),xx(:,3),'rp')
% scatter3(XX(:,1),XX(:,2),XX(:,3),'k.')

scatter3(XX(:,1),XX(:,2),XX(:,3),20,Fval(:,1),'filled',...
    'LineWidth',0.010,...
    'MarkerEdgeAlpha',0.05,'MarkerFaceAlpha',0.05)

s = scatter3(xx(:,1),xx(:,2),xx(:,3),50,fval(:,1),'filled',...
    'MarkerEdgeAlpha',0.5,'MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
box on
grid minor
colormap(c1)
c = colorbar;
c.Label.String = 'obj FE';
view([-39 20])
xlabel('xGDL')
ylabel('xCL')
zlabel('\phi')
%%
fig = figure
hold on

% scatter3(xx(:,1),xx(:,2),xx(:,3),'rp')
% scatter3(XX(:,1),XX(:,2),XX(:,3),'k.')
%
scatter3(XX(:,1),XX(:,2),Fval(:,1),20,Fval(:,1),'filled',...
    'LineWidth',0.010,...
    'MarkerEdgeAlpha',0.05,'MarkerFaceAlpha',0.05)

s = scatter3(xx(:,1),xx(:,2),fval(:,1),50,fval(:,1),'filled',...
    'MarkerEdgeAlpha',0.5,'MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
box on
grid minor
colormap(c1)
c = colorbar;
c.Label.String = 'obj FE';
xlabel('xGDL')
ylabel('xCL')
zlabel('obj FE')
view([-39 20])
%%

[xq,yq] = meshgrid(linspace(min(xx(:,1)),max(xx(:,1)),40),...
    linspace(min(xx(:,2)),max(xx(:,2)),40));
vq = griddata(XX(:,1),XX(:,2),Fval(:,1),xq,yq);
figure
hold on
surf(xq,yq,vq)
s = scatter3(xx(:,1),xx(:,2),fval(:,1),50,fval(:,1),'filled',...
    'MarkerEdgeAlpha',0.5,'MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
box on
grid minor
colormap(c1)
c = colorbar;
c.Label.String = 'obj FE';
xlabel('xGDL')
ylabel('xCL')
zlabel('obj FE')
view([-39 20])

%% Examples

index1 = find(fval(:,1) == max(fval(:,1)));
index1 = index1(xx(index1,1) == min(xx(index1,1)));

index2 = find(abs(xx(:,3)-(-1.6)) == min(abs(xx(:,3)-(-1.6))));

index3 = find(fval(:,3) == max(fval(:,3)));
index3 = index3(xx(index3,1) == min(xx(index3,1)));




%% 9 plot 3D
strnamex = ["Porosity of GDL";"Porosity of CL";"Potential (V \it vs. \rm NHE)"];
% strnamex = '$$ ' + strnamex + ' $$';

strnamey = ["obj FE";"obj EI (GJ/ton CO)";"obj j_{COER} (mA/cm^{2})"];
% strnamey = '$$ ' + strnamey + ' $$';
figPF = figure;
% figPF.Position = [431.4,128.2,1143.2,852.8];
% tiledlayout(3,3)
for i = 1:3
    for j = 1:3
        for k = 1:3
            if j<k
                nexttile
                hold on
                s = scatter3(xx(:,j),xx(:,k),fval(:,i),50,fval(:,i),'filled',...
                    'MarkerEdgeAlpha',0.5,'MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
                xlabel(strnamex(j,1),"Interpreter","tex")
                ylabel(strnamex(k,1),"Interpreter","tex")
                zlabel(strnamey(i,1),"Interpreter","tex")
                view([-39 20])
                box on
                grid minor
                colormap(c1)
            end
        end
    end
end
% print('0Fig31 Pareto set vs objective 3D GA','-djpeg','-r1200')

%% 9 plots 2D
% strnamex = ["x_{\mathrm{GDL}}";"x_{\mathrm{CL}}";"\phi"];
% strnamex = '$$ ' + strnamex + ' $$';
% 
% strnamey = ["obj_\mathrm{FE}";"obj_\mathrm{EI}";"obj_\mathrm{jCOER}"];
% strnamey = '$$ ' + strnamey + ' $$';
% figPF = figure;
% figPF.Position = [431.4,128.2,1143.2,852.8];
% tiledlayout(3,3)
for i = 1:3
    for j = 1:3
        figure
%         nexttile
        hold on
        box on
        oripic = scatter(XX(:,i),Fval(:,j),'o','SizeData',50,...
            'LineWidth',0.1,'MarkerEdgeColor','none','MarkerFaceColor',[0.444705882352941,0.692922352941176,0.856178823529412],...
            'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.2);
        
        PSpic  = scatter(xx(:,i),fval(:,j),'o','SizeData',60,...
            'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.4470 0.7410],...
            'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7);
       PSpicmaxFE  = scatter(xx(index1,i),fval(index1,j),'s','SizeData',65,...
            'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','r',...
            'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',1);
       PSpicsetphi  = scatter(xx(index2,i),fval(index2,j),'^','SizeData',65,...
            'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','m',...
            'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',1);
       PSpicmaxjCOER  = scatter(xx(index3,i),fval(index3,j),'d','SizeData',65,...
            'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','c',...
            'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',1);
  
        set(gca,'FontSize',18);
        xlabel(strnamex(i,1),"Interpreter","tex","FontSize",18)
        ylabel(strnamey(j,1),"Interpreter","tex","FontSize",18)
        legend([oripic,PSpic],["data set","PF by PS"],'location','southwest',...
            'FontSize',15,'Interpreter','tex');
        
        print(['0Fig31 Pareto set vs objective 2D PS ',num2str(i),num2str(j)],'-djpeg','-r1200')

    end
end



















