%% WADCP PACKAGE 0.3 -- function backscatter
%%
%
% Backscatter signal processing toolbox
%
%

function plot_adcp

global nbens theta WN WS WB IMM ens_int binpos yrange xrange binpos_bt inst_mode
global time ve vn vu dir spd a1 amean RL_grid BI Mgl
global Batt Tb_adcp Depth_adcp Sb_adcp
global filein_name path_routines filein_path

tm_min=min(time);
tm_max=max(time);
dt_tick=19;
fsize=14;
if inst_mode==0
maxdepth=max(Depth_adcp);
zaff=binpos;
yrange_aff=yrange;
else
maxdepth=Inf;
zaff=binpos_bt;
yrange_aff=yrange(end)-yrange;
end

maxu=0.8; % max vitesse en m/s
maxMgl=max(Mgl(1,:));

%%
% * hydrodynamics
cd(filein_path)
figure('Name',strcat(filein_name(1:end-3),' -- Hydrodynamics'))
subplot(3,1,1)
imagesc(xrange,zaff,ve'),axis xy, title('Eastward current velocity [m/s]','FontSize',fsize)
ylabel('heigth a.b. [m]'), xlabel(num2str(max(datevec(xrange(1)))));
axis([tm_min tm_max 0 maxdepth]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');
caxis([-maxu maxu]); colorbar('East')
if inst_mode==0
hold on;plot(time,Depth_adcp,'w');axis([tm_min tm_max 0 maxdepth]);
end
hold off;
subplot(3,1,2)
imagesc(xrange,zaff,vn'),axis xy, title('Northward current velocity (m/s)','FontSize',fsize)
xlabel('Time','FontSize',fsize), ylabel('Current velocity (m/s)','FontSize',fsize);axis([tm_min tm_max 0 binpos(end)]);
ylabel('heigth a.b. [m]'), xlabel(num2str(max(datevec(xrange(1)))));
axis([tm_min tm_max 0 maxdepth]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');
caxis([-maxu maxu]); colorbar('East')
if inst_mode==0
hold on;plot(time,Depth_adcp,'w');axis([tm_min tm_max 0 maxdepth]);
end
hold off
subplot(3,1,3)
imagesc(xrange,zaff,vu'),axis xy, title('Upward current velocity (m/s)','FontSize',fsize)
ylabel('heigth a.b. [m]'), xlabel(num2str(max(datevec(xrange(1)))));
axis([tm_min tm_max 0 maxdepth]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');
caxis([-maxu/10 maxu/10]); colorbar('East')
if inst_mode==0
hold on; plot(time,Depth_adcp,'w');axis([tm_min tm_max 0 maxdepth]);
end
hold off

%%
% * hydrodynamics and temperature

figure('Name',strcat(filein_name(1:end-3),' -- Hydrodynamics and temperature'))
subplot(3,1,1)
imagesc(xrange,zaff,spd'),axis xy, title('Current velocity [m/s]','FontSize',fsize)
ylabel('heigth a.b. [m]'), xlabel(num2str(max(datevec(xrange(1)))));
axis([tm_min tm_max 0 maxdepth]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');
caxis([0 maxu]); colorbar('East')
if inst_mode==0
hold on;plot(time,Depth_adcp,'w');axis([tm_min tm_max 0 maxdepth]);
end
hold off;
subplot(3,1,2)
imagesc(xrange,zaff,dir'),axis xy, title('Current Direction [Deg]','FontSize',fsize)
ylabel('heigth a.b. [m]'), xlabel(num2str(max(datevec(xrange(1)))));
axis([tm_min tm_max 0 maxdepth]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');
caxis([0 360]); colorbar('East')
if inst_mode==0
hold on; plot(time,Depth_adcp,'w');axis([tm_min tm_max 0 maxdepth]);
end
hold off

subplot(3,1,3)
plot(time,Tb_adcp), title('Bottom temperature (DegC)','FontSize',fsize)
ylabel('T [DegC]'), xlabel(num2str(max(datevec(xrange(1)))));
axis([tm_min tm_max 5 15]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');

%%
% * hydrodynamics and BI


figure('Name',strcat(filein_name(1:end-3),' -- Hydrodynamics, Backscatter and SPM '))
subplot(3,1,1)
imagesc(xrange,zaff,spd');axis xy, title('Current velocity [m/s]','FontSize',fsize)
ylabel('heigth a.b. [m]'), xlabel(num2str(max(datevec(xrange(1)))));
axis([tm_min tm_max 0 maxdepth]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');
caxis([0 0.5]); colorbar('East')
if inst_mode==0
hold on; plot(time,Depth_adcp,'w');axis([tm_min tm_max 0 maxdepth]);
end
hold off

subplot(3,1,2)
imagesc(xrange,yrange_aff,BI), axis xy, title('Backscatter index [dB ref.1 m^3]','FontSize',fsize)
ylabel('heigth a.b. [m]'), xlabel(num2str(max(datevec(xrange(1)))));
axis([tm_min tm_max 0 maxdepth]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');
caxis([fix(min(min(BI))) fix(max(BI(1,:)))]); colorbar('East')
if inst_mode==0
hold on; plot(time,Depth_adcp,'w');axis([tm_min tm_max 0 maxdepth]);
end
hold off;
subplot(3,1,3)
imagesc(xrange,yrange_aff,Mgl*1000), axis xy, title('SPM concentration [mg/l]','FontSize',fsize)
ylabel('heigth a.b. [m]'), xlabel(num2str(max(datevec(xrange(1))))); 
axis([tm_min tm_max 0 maxdepth]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');
caxis([0 maxMgl*1000]);colorbar('East')
if inst_mode==0
hold on; plot(time,Depth_adcp,'w');axis([tm_min tm_max 0 maxdepth]);
end
hold off



%%

% * Backscatter and SPM


figure('Name',strcat(filein_name(1:end-3),' -- Amplitude, Backscatter and SPM'))
subplot(3,1,1)
imagesc(xrange,yrange_aff,RL_grid);axis xy, title('Amplitude [dB]','FontSize',fsize)
ylabel('heigth a.b. [m]'), xlabel(num2str(max(datevec(xrange(1)))));
axis([tm_min tm_max 0 maxdepth]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');
%caxis([0 0.5]); colorbar('East')
caxis([fix(min(min(RL_grid))) fix(max(RL_grid(1,:)))]); colorbar('East')
if inst_mode==0
hold on; plot(time,Depth_adcp,'w');axis([tm_min tm_max 0 maxdepth]);
end
hold off

subplot(3,1,2)
imagesc(xrange,yrange_aff,BI), axis xy, title('Backscatter index [dB ref.1 m^3]','FontSize',fsize)
ylabel('heigth a.b. [m]'), xlabel(num2str(max(datevec(xrange(1)))));
axis([tm_min tm_max 0 maxdepth]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');
caxis([fix(min(min(BI))) fix(max(BI(1,:)))]); colorbar('East')
if inst_mode==0
hold on; plot(time,Depth_adcp,'w');axis([tm_min tm_max 0 maxdepth]);
end
hold off;
subplot(3,1,3)
imagesc(xrange,yrange_aff,Mgl*1000), axis xy, title('SPM concentration [mg/l]','FontSize',fsize)
ylabel('heigth a.b. [m]'), xlabel(num2str(max(datevec(xrange(1))))); 
axis([tm_min tm_max 0 maxdepth]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)],'XMinorTick','on','TickDir','out','FontSize',fsize,'Box','on');
datetick('x',dt_tick,'keeplimits','keepticks');
caxis([0 maxMgl*1000]);colorbar('East')
if inst_mode==0
hold on; plot(time,Depth_adcp,'w');axis([tm_min tm_max 0 maxdepth]);
end
hold off


cd(path_routines)


