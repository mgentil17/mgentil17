function plot_QA_QC(glid_nav,glid_sci,r_adcp,nav,sci,opt,c_curstruct,c_turbstruct,dpath)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Control Quality Plot on Glider dataset
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Glider variables
x_nav = glid_nav.data(:,1);
xi_nav = nav.data(:,1);

x_sci = glid_sci.data(:,1);
xi_sci = sci.data(:,1);
xi_opt = opt.data(:,1);

% ADCP variables
% vel_raw = find(~isnan(r_adcp.profiles.vel1));
vel_raw = find(~isnan(r_adcp.profiles.EAcntB1));
vel_clean = find(~isnan(c_curstruct.pro.vel1));
turb_raw = find(~isnan(r_adcp.profiles.EAcntB3));
turb_clean =  find(~isnan(c_turbstruct.pro.BI3));

Xraw = [length(x_nav) length(x_sci) length(x_sci) length(vel_raw) length(turb_raw)];
Xpro = [length(xi_nav) length(xi_sci) length(xi_opt) length(vel_clean) length(turb_clean)];
Var = [Xraw(1) Xpro(1); Xraw(2) Xpro(2); Xraw(3) Xpro(3); Xraw(4) Xpro(4); Xraw(5) Xpro(5)];

% Color
pink = rgb('DeepPink');
green = rgb('YellowGreen');

% Plot histogram
m = figure('units','normalized','outerposition',[0 0 1 1]);
width = 0.6;
h1 = barh(Var,width);
grid minor;

labels1 = string(h1(1).YData);
for i = 1:length(Var)
    text(Var(i,1), i-0.25,labels1(i),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
end
labels2 = string(h1(2).YData);
for i = 1:length(Var)
    text(Var(i,2), i+0.05,labels2(i),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
end

% % xtips1 = h1(1).YEndPoints + 0.3;
% % ytips1 = h1(1).XEndPoints;
% % labels1 = string(h1(1).YData);
% % text(xtips1,ytips1,labels1,'VerticalAlignment','middle');
% % 
% % xtips1 = h1(2).YEndPoints + 0.3;
% % ytips1 = h1(2).XEndPoints;
% % labels1 = string(h1(2).YData);
% % text(xtips1,ytips1,labels1,'VerticalAlignment','middle');

h1(2).FaceColor = green;
h1(2).EdgeColor = green;
h1(2).FaceAlpha = 0.5;
h1(1).FaceColor = pink;
h1(1).EdgeColor = pink;
h1(1).FaceAlpha = 0.5;

yticklabels({'Navigation data','CTD data','Bio-optical data','Current data','Turbidity data'});
% xlim([0 2.6e5]);
xlabel('Number of data','Fontweight','bold');
ytickangle(45);
a = get(gca,'YTick');
set(gca,'YTick',a,'fontsize',11,'FontWeight','bold');
legend({'Raw','Pre-processed'},'Location','Southeast');
title({'QA/QC of glider data', '(navigation, science, and acoustic bay)'});

% Save plot
cd(dpath);
saveas(m,'QA_QC_glider_data.jpeg');

end