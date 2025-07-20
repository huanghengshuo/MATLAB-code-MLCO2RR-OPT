clc
clear all
close all
load Pareto3varsGA.mat
load c.mat

fig = figure;
hold on

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

%% 9 plot
strnamex = ["x_{\mathrm{GDL}}";"x_{\mathrm{CL}}";"\phi"];
strnamex = '$$ ' + strnamex + ' $$';

strnamey = ["obj_\mathrm{FE}";"obj_\mathrm{EI}";"obj_\mathrm{jCOER}"];
strnamey = '$$ ' + strnamey + ' $$';
figPF = figure;
figPF.Position = [431.4,128.2,1143.2,852.8];
tiledlayout(3,3)
for i = 1:3
    for j = 1:3
        for k = 1:3
            if j<k
                nexttile
                hold on
                s = scatter3(xx(:,j),xx(:,k),fval(:,i),50,fval(:,i),'filled',...
                    'MarkerEdgeAlpha',0.5,'MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
                xlabel(strnamex(j,1),"Interpreter","latex")
                ylabel(strnamex(k,1),"Interpreter","latex")
                zlabel(strnamey(i,1),"Interpreter","latex")
                view([-39 20])
                box on
                grid minor
                colormap(c1)
            end
        end
    end
end
print('0Fig31 Pareto set vs objective 3D GA','-djpeg','-r1200')

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
        oripic = scatter(XX(:,i),Fval(:,j),'o','SizeData',10,...
            'LineWidth',0.1,'MarkerEdgeColor','none','MarkerFaceColor',[0.4940 0.1840 0.5560],...
            'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.7);
        
        PSpic  = scatter(xx(:,i),fval(:,j),'o','SizeData',20,...
            'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor',[0.3010 0.7450 0.9330],...
            'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7);
      
        
        xlabel(strnamex(i,1),"Interpreter","latex")
        ylabel(strnamey(j,1),"Interpreter","latex")
        legend([oripic,PSpic],["$$ \mathrm{data\ set} $$","$$ \mathrm{GA} $$"],'location','southwest',...
            'FontSize',10,'Interpreter','latex');
    end
end
print('0Fig31 Pareto set vs objective 2D GA','-djpeg','-r1200')



















