function plot_interp(nav,sci,opt,i_glidstruct,dpath)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Plot of Data Synchronization and Interpolation
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

h = figure('units','normalized','outerposition',[0 0 1 1]);

% Colors
pink = rgb('DeepPink');
green = rgb('YellowGreen');

% Attitude variables
idok = find(~isnan(nav.data(:,11)));
x = nav.data(idok,1);
xi = i_glidstruct.Dac(:,1);
y = nav.data(idok,11);              % pitch
z = nav.data(idok,12);              % head
w = nav.data(idok,13);              % roll
yi = i_glidstruct.nav(:,5);
zi = i_glidstruct.nav(:,6);
wi = i_glidstruct.nav(:,7);

subplot 321;
plot(x,y,'o','Color',green,'Markersize',4); hold on; grid minor
plot(xi,yi,'-','Color',pink,'Linewidth',2);
xlim([xi(1) xi(2000)]);
ylabel('Pitch [°]','Fontweight','bold');
legend('pitch','pitch_i_n_t_e_r_p');
ax = gca;
ax.FontSize = 12;
set(ax,'TickDir','out');
set(ax,'xticklabel',[]);

subplot 322;
plot(x,z,'o','Color',green,'Markersize',4); hold on; grid minor
plot(xi,zi,'-','Color',pink,'Linewidth',2);
xlim([xi(1) xi(2000)]);
ylabel('Heading [°]','Fontweight','bold');
legend('head','head_i_n_t_e_r_p');
ax = gca;
ax.FontSize = 12;
set(ax,'TickDir','out');
set(ax,'xticklabel',[]);

subplot 323;
plot(x,w,'o','Color',green,'Markersize',4); hold on; grid minor
plot(xi,wi,'-','Color',pink,'Linewidth',2);
xlim([xi(1) xi(2000)]);
ylabel('Roll [°]','Fontweight','bold');
legend('roll','roll_i_n_t_e_r_p');
ax = gca;
ax.FontSize = 12;
set(ax,'TickDir','out');
set(ax,'xticklabel',[]);

clearvars x xi y yi z zi w wi

% Science variables
idok = find(~isnan(sci.data(:,2)));
x = sci.data(idok,1);
y = sci.data(idok,2);              % temperature
xi = i_glidstruct.sci(:,1);
yi = i_glidstruct.sci(:,2);

subplot 324;
plot(x,y,'o','Color',green,'Markersize',4); hold on; grid minor
plot(xi,yi,'-','Color',pink,'Linewidth',2);
xlim([xi(1) xi(2000)]);
ylabel('Temperature [°C]','Fontweight','bold');
legend('T','T_i_n_t_e_r_p');
ax = gca;
ax.FontSize = 12;
set(ax,'TickDir','out');
set(ax,'xticklabel',[]);

clearvars x xi y yi

% Bio-optical variables
idok = find(~isnan(opt.data(:,2)));
x = opt.data(idok,1);
y = opt.data(idok,3);              % BB
z = opt.data(idok,2);              % Chla
xi = i_glidstruct.sci(:,1);
yi = i_glidstruct.sci(:,6);
zi = i_glidstruct.sci(:,5);

subplot 325;
plot(x,y,'o','Color',green,'Markersize',4); hold on; grid minor
plot(xi,yi,'-','Color',pink,'Linewidth',2);
xlim([xi(1) xi(2000)]);
ylabel('Backscatter [m^-^1]','Fontweight','bold');
legend('BB','BB_i_n_t_e_r_p');
ax = gca;
ax.FontSize = 12;
set(ax,'TickDir','out');
datetick('x','HH:MM','keepticks');

subplot 326;
plot(x,z,'o','Color',green,'Markersize',4); hold on; grid minor
plot(xi,zi,'-','Color',pink,'Linewidth',2);
xlim([xi(1) xi(2000)]);
ylabel('Chlorophyll-a [\mug L^-^1]','Fontweight','bold');
legend('Chla','Chla_i_n_t_e_r_p');
ax = gca;
ax.FontSize = 12;
set(ax,'TickDir','out');
datetick('x','HH:MM','keepticks');

% Save figure
cd(dpath);
saveas(h,'interp_variables.jpeg');

end