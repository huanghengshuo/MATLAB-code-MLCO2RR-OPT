% Pearson & Spearman ori
% Coorrelation coefficients
% Plot more figures

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

phix = unique(MFE(:,3))';
for i = 2:2:length(phix)-1
    
    index = find(abs(MFE(:,3)-phix(i))<0.001);
    MFE(index,:) = [];
    MCOER(index,:) = [];
end

phix = unique(MFE(:,3))';
rhoP = partialcorr(MFE,'Type','Pearson');
rhoS = partialcorr(MFE,'Type','Spearman');



xvalues = {'epiGDL','epiCL','phi','FECO'};
yvalues = xvalues;
f1 = figure;
h1 = heatmap(f1,xvalues,yvalues,rhoP,'GridVisible','off');

f2 = figure;
h2 = heatmap(f2,xvalues,yvalues,rhoS,'GridVisible','off');


c = [0.9 0.9447 0.9741;0.896470588235294 0.942748235294118 0.973185882352941;0.892941176470588 0.940796470588235 0.972271764705882;0.889411764705882 0.938844705882353 0.971357647058824;0.885882352941176 0.936892941176471 0.970443529411765;0.882352941176471 0.934941176470588 0.969529411764706;0.878823529411765 0.932989411764706 0.968615294117647;0.875294117647059 0.931037647058824 0.967701176470588;0.871764705882353 0.929085882352941 0.966787058823529;0.868235294117647 0.927134117647059 0.96587294117647;0.864705882352941 0.925182352941177 0.964958823529412;0.861176470588235 0.923230588235294 0.964044705882353;0.857647058823529 0.921278823529412 0.963130588235294;0.854117647058824 0.919327058823529 0.962216470588235;0.850588235294118 0.917375294117647 0.961302352941176;0.847058823529412 0.915423529411765 0.960388235294118;0.843529411764706 0.913471764705882 0.959474117647059;0.84 0.91152 0.95856;0.836470588235294 0.909568235294118 0.957645882352941;0.832941176470588 0.907616470588235 0.956731764705882;0.829411764705882 0.905664705882353 0.955817647058823;0.825882352941177 0.903712941176471 0.954903529411765;0.822352941176471 0.901761176470588 0.953989411764706;0.818823529411765 0.899809411764706 0.953075294117647;0.815294117647059 0.897857647058824 0.952161176470588;0.811764705882353 0.895905882352941 0.951247058823529;0.808235294117647 0.893954117647059 0.950332941176471;0.804705882352941 0.892002352941176 0.949418823529412;0.801176470588235 0.890050588235294 0.948504705882353;0.797647058823529 0.888098823529412 0.947590588235294;0.794117647058824 0.88614705882353 0.946676470588235;0.790588235294118 0.884195294117647 0.945762352941176;0.787058823529412 0.882243529411765 0.944848235294118;0.783529411764706 0.880291764705882 0.943934117647059;0.78 0.87834 0.94302;0.776470588235294 0.876388235294118 0.942105882352941;0.772941176470588 0.874436470588235 0.941191764705882;0.769411764705882 0.872484705882353 0.940277647058824;0.765882352941176 0.870532941176471 0.939363529411765;0.762352941176471 0.868581176470588 0.938449411764706;0.758823529411765 0.866629411764706 0.937535294117647;0.755294117647059 0.864677647058824 0.936621176470588;0.751764705882353 0.862725882352941 0.935707058823529;0.748235294117647 0.860774117647059 0.934792941176471;0.744705882352941 0.858822352941176 0.933878823529412;0.741176470588235 0.856870588235294 0.932964705882353;0.737647058823529 0.854918823529412 0.932050588235294;0.734117647058824 0.852967058823529 0.931136470588235;0.730588235294118 0.851015294117647 0.930222352941177;0.727058823529412 0.849063529411765 0.929308235294118;0.723529411764706 0.847111764705882 0.928394117647059;0.72 0.84516 0.92748;0.716470588235294 0.843208235294118 0.926565882352941;0.712941176470588 0.841256470588235 0.925651764705882;0.709411764705882 0.839304705882353 0.924737647058824;0.705882352941176 0.837352941176471 0.923823529411765;0.702352941176471 0.835401176470588 0.922909411764706;0.698823529411765 0.833449411764706 0.921995294117647;0.695294117647059 0.831497647058824 0.921081176470588;0.691764705882353 0.829545882352941 0.920167058823529;0.688235294117647 0.827594117647059 0.919252941176471;0.684705882352941 0.825642352941176 0.918338823529412;0.681176470588235 0.823690588235294 0.917424705882353;0.677647058823529 0.821738823529412 0.916510588235294;0.674117647058823 0.819787058823529 0.915596470588235;0.670588235294118 0.817835294117647 0.914682352941177;0.667058823529412 0.815883529411765 0.913768235294118;0.663529411764706 0.813931764705882 0.912854117647059;0.66 0.81198 0.91194;0.656470588235294 0.810028235294118 0.911025882352941;0.652941176470588 0.808076470588235 0.910111764705882;0.649411764705882 0.806124705882353 0.909197647058824;0.645882352941176 0.804172941176471 0.908283529411765;0.642352941176471 0.802221176470588 0.907369411764706;0.638823529411765 0.800269411764706 0.906455294117647;0.635294117647059 0.798317647058824 0.905541176470588;0.631764705882353 0.796365882352941 0.904627058823529;0.628235294117647 0.794414117647059 0.903712941176471;0.624705882352941 0.792462352941177 0.902798823529412;0.621176470588235 0.790510588235294 0.901884705882353;0.617647058823529 0.788558823529412 0.900970588235294;0.614117647058824 0.786607058823529 0.900056470588235;0.610588235294118 0.784655294117647 0.899142352941176;0.607058823529412 0.782703529411765 0.898228235294118;0.603529411764706 0.780751764705882 0.897314117647059;0.6 0.7788 0.8964;0.596470588235294 0.776848235294118 0.895485882352941;0.592941176470588 0.774896470588235 0.894571764705882;0.589411764705882 0.772944705882353 0.893657647058824;0.585882352941176 0.770992941176471 0.892743529411765;0.582352941176471 0.769041176470588 0.891829411764706;0.578823529411765 0.767089411764706 0.890915294117647;0.575294117647059 0.765137647058824 0.890001176470588;0.571764705882353 0.763185882352941 0.889087058823529;0.568235294117647 0.761234117647059 0.888172941176471;0.564705882352941 0.759282352941177 0.887258823529412;0.561176470588235 0.757330588235294 0.886344705882353;0.557647058823529 0.755378823529412 0.885430588235294;0.554117647058823 0.753427058823529 0.884516470588235;0.550588235294118 0.751475294117647 0.883602352941177;0.547058823529412 0.749523529411765 0.882688235294118;0.543529411764706 0.747571764705882 0.881774117647059;0.54 0.74562 0.88086;0.536470588235294 0.743668235294118 0.879945882352941;0.532941176470588 0.741716470588235 0.879031764705882;0.529411764705882 0.739764705882353 0.878117647058823;0.525882352941176 0.737812941176471 0.877203529411765;0.522352941176471 0.735861176470588 0.876289411764706;0.518823529411765 0.733909411764706 0.875375294117647;0.515294117647059 0.731957647058823 0.874461176470588;0.511764705882353 0.730005882352941 0.873547058823529;0.508235294117647 0.728054117647059 0.872632941176471;0.504705882352941 0.726102352941177 0.871718823529412;0.501176470588235 0.724150588235294 0.870804705882353;0.497647058823529 0.722198823529412 0.869890588235294;0.494117647058824 0.720247058823529 0.868976470588235;0.490588235294118 0.718295294117647 0.868062352941176;0.487058823529412 0.716343529411765 0.867148235294118;0.483529411764706 0.714391764705882 0.866234117647059;0.48 0.71244 0.86532;0.476470588235294 0.710488235294118 0.864405882352941;0.472941176470588 0.708536470588235 0.863491764705882;0.469411764705882 0.706584705882353 0.862577647058824;0.465882352941176 0.704632941176471 0.861663529411765;0.462352941176471 0.702681176470588 0.860749411764706;0.458823529411765 0.700729411764706 0.859835294117647;0.455294117647059 0.698777647058824 0.858921176470588;0.451764705882353 0.696825882352941 0.858007058823529;0.448235294117647 0.694874117647059 0.857092941176471;0.444705882352941 0.692922352941176 0.856178823529412;0.441176470588235 0.690970588235294 0.855264705882353;0.437647058823529 0.689018823529412 0.854350588235294;0.434117647058823 0.687067058823529 0.853436470588235;0.430588235294118 0.685115294117647 0.852522352941177;0.427058823529412 0.683163529411765 0.851608235294118;0.423529411764706 0.681211764705882 0.850694117647059;0.42 0.67926 0.84978;0.416470588235294 0.677308235294118 0.848865882352941;0.412941176470588 0.675356470588235 0.847951764705882;0.409411764705882 0.673404705882353 0.847037647058823;0.405882352941176 0.671452941176471 0.846123529411765;0.402352941176471 0.669501176470588 0.845209411764706;0.398823529411765 0.667549411764706 0.844295294117647;0.395294117647059 0.665597647058824 0.843381176470588;0.391764705882353 0.663645882352941 0.842467058823529;0.388235294117647 0.661694117647059 0.841552941176471;0.384705882352941 0.659742352941177 0.840638823529412;0.381176470588235 0.657790588235294 0.839724705882353;0.377647058823529 0.655838823529412 0.838810588235294;0.374117647058824 0.653887058823529 0.837896470588235;0.370588235294118 0.651935294117647 0.836982352941176;0.367058823529412 0.649983529411765 0.836068235294118;0.363529411764706 0.648031764705882 0.835154117647059;0.36 0.64608 0.83424;0.356470588235294 0.644128235294118 0.833325882352941;0.352941176470588 0.642176470588235 0.832411764705882;0.349411764705882 0.640224705882353 0.831497647058824;0.345882352941176 0.63827294117647 0.830583529411765;0.342352941176471 0.636321176470588 0.829669411764706;0.338823529411765 0.634369411764706 0.828755294117647;0.335294117647059 0.632417647058824 0.827841176470588;0.331764705882353 0.630465882352941 0.826927058823529;0.328235294117647 0.628514117647059 0.826012941176471;0.324705882352941 0.626562352941176 0.825098823529412;0.321176470588235 0.624610588235294 0.824184705882353;0.317647058823529 0.622658823529412 0.823270588235294;0.314117647058824 0.620707058823529 0.822356470588235;0.310588235294118 0.618755294117647 0.821442352941176;0.307058823529412 0.616803529411765 0.820528235294118;0.303529411764706 0.614851764705882 0.819614117647059;0.3 0.6129 0.8187;0.296470588235294 0.610948235294118 0.817785882352941;0.292941176470588 0.608996470588235 0.816871764705882;0.289411764705882 0.607044705882353 0.815957647058823;0.285882352941176 0.605092941176471 0.815043529411765;0.282352941176471 0.603141176470588 0.814129411764706;0.278823529411765 0.601189411764706 0.813215294117647;0.275294117647059 0.599237647058824 0.812301176470588;0.271764705882353 0.597285882352941 0.811387058823529;0.268235294117647 0.595334117647059 0.810472941176471;0.264705882352941 0.593382352941177 0.809558823529412;0.261176470588235 0.591430588235294 0.808644705882353;0.257647058823529 0.589478823529412 0.807730588235294;0.254117647058823 0.587527058823529 0.806816470588235;0.250588235294118 0.585575294117647 0.805902352941177;0.247058823529412 0.583623529411765 0.804988235294118;0.243529411764706 0.581671764705882 0.804074117647059;0.24 0.57972 0.80316;0.236470588235294 0.577768235294118 0.802245882352941;0.232941176470588 0.575816470588235 0.801331764705882;0.229411764705882 0.573864705882353 0.800417647058824;0.225882352941176 0.571912941176471 0.799503529411765;0.222352941176471 0.569961176470588 0.798589411764706;0.218823529411765 0.568009411764706 0.797675294117647;0.215294117647059 0.566057647058824 0.796761176470588;0.211764705882353 0.564105882352941 0.795847058823529;0.208235294117647 0.562154117647059 0.794932941176471;0.204705882352941 0.560202352941176 0.794018823529412;0.201176470588235 0.558250588235294 0.793104705882353;0.19764705882353 0.556298823529412 0.792190588235294;0.194117647058824 0.554347058823529 0.791276470588235;0.190588235294118 0.552395294117647 0.790362352941176;0.187058823529412 0.550443529411765 0.789448235294118;0.183529411764706 0.548491764705882 0.788534117647059;0.18 0.54654 0.78762;0.176470588235294 0.544588235294118 0.786705882352941;0.172941176470588 0.542636470588235 0.785791764705882;0.169411764705882 0.540684705882353 0.784877647058824;0.165882352941176 0.538732941176471 0.783963529411765;0.162352941176471 0.536781176470588 0.783049411764706;0.158823529411765 0.534829411764706 0.782135294117647;0.155294117647059 0.532877647058823 0.781221176470588;0.151764705882353 0.530925882352941 0.780307058823529;0.148235294117647 0.528974117647059 0.779392941176471;0.144705882352941 0.527022352941176 0.778478823529412;0.141176470588235 0.525070588235294 0.777564705882353;0.137647058823529 0.523118823529412 0.776650588235294;0.134117647058823 0.521167058823529 0.775736470588235;0.130588235294118 0.519215294117647 0.774822352941176;0.127058823529412 0.517263529411765 0.773908235294118;0.123529411764706 0.515311764705882 0.772994117647059;0.12 0.51336 0.77208;0.116470588235294 0.511408235294118 0.771165882352941;0.112941176470588 0.509456470588235 0.770251764705882;0.109411764705882 0.507504705882353 0.769337647058824;0.105882352941177 0.505552941176471 0.768423529411765;0.102352941176471 0.503601176470588 0.767509411764706;0.0988235294117646 0.501649411764706 0.766595294117647;0.0952941176470588 0.499697647058824 0.765681176470588;0.091764705882353 0.497745882352941 0.764767058823529;0.0882352941176471 0.495794117647059 0.763852941176471;0.0847058823529412 0.493842352941177 0.762938823529412;0.0811764705882353 0.491890588235294 0.762024705882353;0.0776470588235294 0.489938823529412 0.761110588235294;0.0741176470588236 0.487987058823529 0.760196470588235;0.0705882352941176 0.486035294117647 0.759282352941176;0.0670588235294117 0.484083529411765 0.758368235294118;0.0635294117647058 0.482131764705882 0.757454117647059;0.0599999999999999 0.48018 0.75654;0.0564705882352942 0.478228235294118 0.755625882352941;0.0529411764705883 0.476276470588235 0.754711764705882;0.0494117647058824 0.474324705882353 0.753797647058824;0.0458823529411765 0.472372941176471 0.752883529411765;0.0423529411764706 0.470421176470588 0.751969411764706;0.0388235294117647 0.468469411764706 0.751055294117647;0.0352941176470588 0.466517647058824 0.750141176470588;0.0317647058823529 0.464565882352941 0.749227058823529;0.028235294117647 0.462614117647059 0.748312941176471;0.0247058823529411 0.460662352941176 0.747398823529412;0.0211764705882354 0.458710588235294 0.746484705882353;0.0176470588235295 0.456758823529412 0.745570588235294;0.0141176470588236 0.454807058823529 0.744656470588235;0.0105882352941177 0.452855294117647 0.743742352941177;0.00705882352941167 0.450903529411765 0.742828235294118;0.00352941176470589 0.448951764705882 0.741914117647059;0 0.447 0.741];
%% Train tree model
X1 = MFE(:,1:3);
Y1 = MFE(:,4);

% NoKFold = 10;
TreeModelCV = fitrtree(X1,Y1);
% % GaussModelCV = fitrgp(X1,Y1,'CrossVal','on','KFold',NoKFold);
%
fig = figure;
fig.Position = [185,51,1271,912];
p = plotPartialDependence(TreeModelCV,[1,3]);
xlabel('x');
ylabel('y');
view([89.9 90]);
xlim([min(MFE(:,1)),max(MFE(:,1))]);
ylim([min(MFE(:,3)),max(MFE(:,3))]);
colormap(c)
p.Children.EdgeColor = 'none';
p.Children.FaceColor = 'none';

title([],'Visible','off');

figure
p = plotPartialDependence(TreeModelCV,3,"Conditional","absolute");
p.Children(2,1).MarkerFaceColor = 'none';
p.Children(2,1).MarkerEdgeColor = 'none';
% f2 = figure;
% %
% p = partialDependence(TreeModelCV,[1,3]);
% h = heatmap(f2,p,'GridVisible','off');

%% Train Linear Regression Model
epiGDL = MFE(:,1);
epiCL = MFE(:,2);
phi = MFE(:,3);
FE = MFE(:,4);
T = table(epiGDL,epiCL,phi,FE,'VariableNames',{'epiGDL','epiCL','phi','FE'});
lm = fitlm(T,'linear');
figure
p = plotPartialDependence(lm,3,"Conditional","absolute");
p.Children(2,1).MarkerFaceColor = 'none';
p.Children(2,1).MarkerEdgeColor = 'none';
[p.Children(3:end,1).Color] = deal([0.8, 0.8, 0.8]);
[ypred,yci] = predict(lm,MFE(:,1:3));


% figure
% hold on
% plot(phi,FE,'.');
% plot(phi',ypred','r-');
% fill(phi',yci(:,1)');
% fill(phi',yci(:,2)');
fig = figure;
fig.Position = [185,51,1271,912];
t = tiledlayout(3,3);
for i = 1:3
    for j = 1:3
        if i <= j
            if i == j
                nexttile
                hold on
                PDP = plotPartialDependence(lm,i,"Conditional","absolute");
                PDP.Children(2,1).MarkerFaceColor = 'none';
                PDP.Children(2,1).MarkerEdgeColor = 'none';
                [PDP.Children(3:end,1).Color] = deal([0.8, 0.8, 0.8]);
                scatter(MFE(:,i),MFE(:,4),'o',...
                    'MarkerFaceColor',[0.3010 0.7450 0.9330],...
                    'MarkerEdgeColor',[1 1 1],...
                    'MarkerFaceAlpha',0.5,...
                    'MarkerEdgeAlpha',0.5);
            else
                nexttile
                p = plotPartialDependence(lm,[i,j]);
                %                 xlabel('x');
                %                 ylabel('y');
                view([89.9 90]);
                xlim([min(MFE(:,i)),max(MFE(:,i))]);
                ylim([min(MFE(:,j)),max(MFE(:,j))]);
                colormap(c)
                p.Children.EdgeColor = 'none';
            end
        else
            nexttile
            plot(MFE(:,j),MFE(:,i),'o','Color',[0.3010 0.7450 0.9330]);
        end
    end
end








toc
