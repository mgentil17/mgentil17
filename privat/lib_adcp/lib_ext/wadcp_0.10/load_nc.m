function load_nc
path_routines=cd;

[filein_name,filein_path]=uigetfile('*.nc','Name of the netcdf to open:');
cd(filein_path)
%filein_name='huveaune_raw_2_proc.nc';
ncin=netcdf(filein_name);
cd(path_routines);

time_adcp=ncin{'time'}(:);

t_min=min(time_adcp);
t_max=max(time_adcp);

prompt={'tmin', 'tmax'};
dlg_title='Data processing - time limits';
num_lines=1;
def={datestr(t_min),datestr(t_max)};
input1=inputdlg(prompt,dlg_title,num_lines,def);
tinf=datenum(input1(1));
tsup=datenum(input1(2));

timevec=find(time_adcp>tinf & time_adcp<tsup);
time_adcp=time_adcp(timevec);
nbens=length(time_adcp);


binpos=ncin{'binpos'}(:);
vu=ncin{'vu'}(timevec,:);
spd=ncin{'spd'}(timevec,:);
%dir=ncin{'dir'}(timevec,:);
WS=ncin{'WS'}(:);

yrange=(0:length(binpos)*2-1)*WS/2+binpos(1);

%amean=ncin{'amean'}(timevec,:);
IV=ncin{'BI'}(timevec,:);
Mgl=ncin{'Mgl'}(timevec,:);
Tb_adcp=ncin{'Tb_adcp'}(timevec);
Depth_adcp=ncin{'Depth_adcp'}(timevec);
%Batt=ncin{'Batt'}(timevec);
close(ncin)


fsize=10;
dt_tick=19;
tm_min=min(time_adcp);
tm_max=max(time_adcp);


spdmax=1;
hsmax=4;
ivmin=-70;
ivmax=-50;
mglmax=60;
upvelmax=0.1;
zmax=25;



%%
% * hydrodynamics and IV

figure('Name',strcat(filein_name(1:end-3),' -- Hydrodynamics and IV'))
subplot(2,1,1), hold on
imagesc(time_adcp,binpos,spd'), title('Current velocity (m/s)'), xlabel('Time'), ylabel('Current velocity (m/s)','FontSize',fsize);axis([tm_min tm_max 0 zmax]);
colormap(jet(20))
shading flat
set(gca,'XTick',(tm_min:(tm_max-tm_min)/10:tm_max),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tm_min:(tm_max-tm_min)/10:tm_max),dt_tick),'FontSize',fsize,'Box','on');
caxis([0 spdmax]); colorbar('East')
plot(time_adcp,Depth_adcp);axis([tm_min tm_max 0 zmax]);
hold off

% subplot(3,1,2), hold on
% imagesc(time_adcp,binpos,vu'), title('Upward velocity (m/s)'), xlabel('Time'), ylabel('Upward velocity (m/s)','FontSize',fsize);axis([tm_min tm_max 0 zmax]);
% colormap(jet(20))
% shading flat
% set(gca,'XTick',(tm_min:(tm_max-tm_min)/10:tm_max),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tm_min:(tm_max-tm_min)/10:tm_max),dt_tick),'FontSize',fsize,'Box','on');
% caxis([0 upvelmax]); colorbar('East')
% plot(time_adcp,Depth_adcp);axis([tm_min tm_max 0 zmax]);
% hold off
% 

subplot(2,1,2), hold on
imagesc(time_adcp,yrange,Mgl'), title('Concentration en MES (g/l))'), xlabel('Time'), ylabel('Concentration en MES (g/l)','FontSize',fsize);axis([tm_min tm_max 0 zmax]);
set(gca,'XTick',(tm_min:(tm_max-tm_min)/10:tm_max),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tm_min:(tm_max-tm_min)/10:tm_max),dt_tick),'FontSize',fsize,'Box','on');
caxis([0 0.1])
shading flat
colorbar('East')
plot(time_adcp,Depth_adcp);axis([tm_min tm_max 0 zmax]);
hold off

figure
imagesc(Mgl)

