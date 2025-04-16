clc
clear all
close all

load Pareto3varsPS.mat
load c.mat
xxPS = xx;
fvalPS = fval;
load Pareto3varsGA.mat
xxGA = xx;
fvalGA = fval;

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
            'LineWidth',0.1,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1);
        
        PSpic  = scatter(xxPS(:,i),fvalPS(:,j),'o','SizeData',20,...
            'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','c',...
            'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7);
        GApic  = scatter(xxGA(:,i),fvalGA(:,j),'o','SizeData',20,...
            'LineWidth',0.2,'MarkerEdgeColor','k','MarkerFaceColor','m',...
            'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.7);
        
        
        xlabel(strnamex(i,1),"Interpreter","latex")
        ylabel(strnamey(j,1),"Interpreter","latex")
        legend([oripic,GApic,PSpic],["$$ \mathrm{original\ data} $$","$$ \mathrm{GA} $$","$$ \mathrm{PS} $$"],'location','southwest',...
            'FontSize',10,'Interpreter','latex');
    end
end


%% 9 plot
strnamex = ["x_{\mathrm{GDL}}";"x_{\mathrm{CL}}";"\phi"];
strnamex = '$$ ' + strnamex + ' $$';

strnamey = ["obj_\mathrm{FE}";"obj_\mathrm{EI}";"obj_\mathrm{jCOER}"];
strnamey = '$$ ' + strnamey + ' $$';
figPF = figure;
figPF.Position = [431.4,128.2,1143.2,852.8];
T = tiledlayout(3,3)
for i = 1:3
    for j = 1:3
        for k = 1:3
            if j<k

                
                ax1 = nexttile(T);
                
                s1 = scatter3(ax1,xxGA(:,j),xxGA(:,k),fvalGA(:,i),50,fvalGA(:,i),'filled','Marker','o',...
                    'MarkerEdgeAlpha',0.5,'MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
                colormap(ax1,c1);
                view([-39 20])
%                 axis off
                hold on
                ax2 = axes('Position',ax1.Position);
                s2 = scatter3(ax2,xxPS(:,j),xxPS(:,k),fvalPS(:,i),50,fvalPS(:,i),'filled','Marker','s',...
                    'MarkerEdgeAlpha',0.5,'MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
                colormap(ax2,hot);
                ax2.Visible = 'off';
                xlabel(ax1,strnamex(j,1),"Interpreter","latex")
                ylabel(ax1,strnamex(k,1),"Interpreter","latex")
                zlabel(ax1,strnamey(i,1),"Interpreter","latex")
                view([-39 20])
                box on
                grid minor
                %                 colormap(c1)
                %                 legend([s1,s2],["$$ \mathrm{GA} $$","$$ \mathrm{PS} $$"],'location','southwest',...
                %                     'FontSize',10,'Interpreter','latex');
%                 hold off
            end
        end
    end
end

