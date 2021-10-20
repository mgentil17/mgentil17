function load_nc_previmer_lecroizic
warning off all

path_routines=cd;


%%
% * User-defined colors for all plots
icolor(:,1)=[162 38 109]/255;
icolor(:,2)=[0 96 255]/255;
icolor(:,3)=[80 162 84]/255;
icolor(:,4)=[192 119 38]/255;
icolor(:,5)=[217 23 222]/255;
icolor(:,6)=[167 211 247]/255;


[filein_name,filein_path]=uigetfile('*.mat','Name of the matlab file to open:');
cd(filein_path)
load(filein_name);


time=data.time.val;

t_min=min(time);
t_max=max(time);


prompt={'tmin', 'tmax'};
dlg_title='Data processing - time limits';
num_lines=1;
def={datestr(t_min),datestr(t_max)};
input1=inputdlg(prompt,dlg_title,num_lines,def);
tinf=datenum(input1(1));
tsup=datenum(input1(2));



timevec=find(time>tinf & time<tsup);
timevec=timevec(1:end-1);
time_adcp=time(timevec);

nbens=length(time);


spd=data.spd.val;
spd=spd(timevec,:);
dir=data.dir.val;
dir=dir(timevec,:);

Depth_adcp=data.Depth_adcp.val;
Depth_adcp=Depth_adcp(timevec,:);

yrange=data.yrange.val;

Mgl=data.Mgl.val;
Mgl=Mgl(timevec,:);



fsize=14;
dt_tick=19;
tm_min=min(time_adcp);
tm_max=max(time_adcp);


spdmax=1;
hsmax=4;
Mglmax=0.10;
zmax=25;


%%
% charger fichier houle

file=xlsread('PMD4_Le_Croisic_1_houle.xls');
time_houle=datenum(file(:,3),file(:,2),file(:,1))+file(:,4);
HS=file(:,6);

%%
% charger fichier seapoint
file=xlsread('PMD4_Le_Croisic_1_seapoint_fond.xls');
time_seapoint=datenum(file(:,3),file(:,2),file(:,1))+file(:,4);
seapoint=file(:,7)/1000; % passage mg/l ->g/l

%%
% charger fichier MS5
file=xlsread('PMD4_Le_Croisic_1_MS5_subsurface.xls');
time_MS5=datenum(file(:,3),file(:,2),file(:,1))+file(:,4);
MS5=file(:,9)/1000;

MS5_adcp=interp1(time_MS5,MS5,time_adcp);
seapoint_adcp=interp1(time_seapoint,seapoint,time_adcp);



for i=1:length(time_adcp);
  
    Mgl_MS5(i)=interp1(yrange,Mgl(i,:),Depth_adcp(i)-1.70);
    Mgl_MS5_3m(i)=interp1(yrange,Mgl(i,:),Depth_adcp(i)-3);
    Mgl_MS5_5m(i)=interp1(yrange,Mgl(i,:),Depth_adcp(i)-5);
    Mgl_MS5_8m(i)=interp1(yrange,Mgl(i,:),Depth_adcp(i)-8);
    
    toto=find(isnan(Mgl(i,:)));
    depth_Mgl(i)=yrange(toto(1));
    
end



%%
% * hydrodynamics and IV
% 
% figure('Name',strcat(filein_name(1:end-3),' -- Hydrodynamics and IV'))
% subplot(2,1,1), hold on
% imagesc(time_adcp,binpos,spd'), title('Current velocity (m/s)'), xlabel('Time'), ylabel('Current velocity (m/s)','FontSize',fsize);axis([tm_min tm_max 0 zmax]);
% colormap(jet(20))
% shading flat
% set(gca,'XTick',(tm_min:(tm_max-tm_min)/10:tm_max),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tm_min:(tm_max-tm_min)/10:tm_max),dt_tick),'FontSize',fsize,'Box','on');
% caxis([0 spdmax]); colorbar('East')
% plot(time_adcp,Depth_adcp);axis([tm_min tm_max 0 zmax]);
% hold off

% subplot(3,1,2), hold on
% imagesc(time_adcp,binpos,vu'), title('Upward velocity (m/s)'), xlabel('Time'), ylabel('Upward velocity (m/s)','FontSize',fsize);axis([tm_min tm_max 0 zmax]);
% colormap(jet(20))
% shading flat
% set(gca,'XTick',(tm_min:(tm_max-tm_min)/10:tm_max),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tm_min:(tm_max-tm_min)/10:tm_max),dt_tick),'FontSize',fsize,'Box','on');
% caxis([0 upvelmax]); colorbar('East')
% plot(time_adcp,Depth_adcp);axis([tm_min tm_max 0 zmax]);
% hold off
% 
figure
subplot(2,1,1), hold on
plot(time_houle,HS), title('Hauteur significative (m)','FontSize',fsize), ylabel('Hauteur significative (m)','FontSize',fsize);axis([tm_min tm_max 0 hsmax]);
set(gca,'XTick',(tm_min:(tm_max-tm_min)/10:tm_max),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tm_min:(tm_max-tm_min)/10:tm_max),dt_tick),'FontSize',fsize,'Box','on');

subplot(2,1,2), hold on
imagesc(time_adcp,yrange,Mgl'), title('Concentration en MES (g/l)','FontSize',fsize), xlabel('Time'), ylabel('Concentration en MES (g/l)','FontSize',fsize);axis([tm_min tm_max 0 zmax]);
set(gca,'XTick',(tm_min:(tm_max-tm_min)/10:tm_max),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tm_min:(tm_max-tm_min)/10:tm_max),dt_tick),'FontSize',fsize,'Box','on');
caxis([0 Mglmax])
shading flat
colorbar('East')
plot(time_adcp,Depth_adcp,'w');axis([tm_min tm_max 0 zmax]);
set(gca,'XTick',(tm_min:(tm_max-tm_min)/10:tm_max),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tm_min:(tm_max-tm_min)/10:tm_max),dt_tick),'FontSize',fsize,'Box','on');
hold off

figure
subplot(3,1,1), hold on
plot(time_houle,HS), title('Hauteur significative (m)','FontSize',fsize), ylabel('Hauteur significative (m)','FontSize',fsize);axis([tm_min tm_max 0 hsmax]);
set(gca,'XTick',(tm_min:(tm_max-tm_min)/10:tm_max),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tm_min:(tm_max-tm_min)/10:tm_max),dt_tick),'FontSize',fsize,'Box','on');


subplot(3,1,2), hold on
plot(time_MS5,MS5,'Color',icolor(:,3)), title('Concentration en MES (g/l)','FontSize',fsize), ylabel('Concentration en MES (g/l)','FontSize',fsize);axis([tm_min tm_max 0 0.1]);
plot(time_adcp,Mgl_MS5,'Color',icolor(:,2))
plot(time_adcp,Mgl_MS5_3m,'Color',icolor(:,4))
plot(time_adcp,Mgl_MS5_5m,'Color',icolor(:,5))
set(gca,'XTick',(tm_min:(tm_max-tm_min)/10:tm_max),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tm_min:(tm_max-tm_min)/10:tm_max),dt_tick),'FontSize',fsize,'Box','on');
hold off
legend('MS5 subsurface','ADCP 1.7m sous surface', 'ADCP 3m sous surface','ADCP 5m sous surface')

subplot(3,1,3), hold on
plot(time_seapoint,seapoint,'Color',icolor(:,1)), title('Concentration en MES (g/l)','FontSize',fsize), ylabel('Concentration en MES (g/l)','FontSize',fsize);axis([tm_min tm_max 0 0.1]);
plot(time_adcp,Mgl(:,1),'Color',icolor(:,2))
set(gca,'XTick',(tm_min:(tm_max-tm_min)/10:tm_max),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tm_min:(tm_max-tm_min)/10:tm_max),dt_tick),'FontSize',fsize,'Box','on');
hold off
legend('Seapoint fond','ADCP')

% 
% figure
% subplot(1,2,1)
% hold on
% plot(Mgl_MS5,MS5_adcp,'.','Color',icolor(:,1))
% plot(Mgl_MS5_3m,MS5_adcp,'.','Color',icolor(:,2))
% plot(Mgl_MS5_5m,MS5_adcp,'.','Color',icolor(:,3))
% hold off
% subplot(1,2,2)
% plot(Mgl(:,1),seapoint_adcp,'.')
% xlabel('Concentration en MES ADCP (g/l)');ylabel('Concentration en MES Seapoint (g/l)'); 

% 
% figure
% hold on
% plot(time_adcp,Depth_adcp)
% plot(time_adcp,depth_Mgl,'r')
% 
% out_data(:,1:6)=datevec(time_adcp);
% out_data(:,7)=MS5_adcp;
% out_data(:,8)=seapoint_adcp;
% out_data(:,9)=Mgl_MS5;
% out_data(:,10)=Mgl_MS5_3m;
% out_data(:,11)=Mgl_MS5_5m;
% out_data(:,12)=Mgl_MS5_8m;
% out_data(:,13)=Mgl(:,1);
% save('PMD4_le_croizic_ADCP_surface.txt','out_data','-ASCII');


cd(path_routines);
