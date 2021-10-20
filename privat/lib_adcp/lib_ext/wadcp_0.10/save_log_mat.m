function save_log_mat

global fileout_name
global path_routines
global matsave

filein_name=fileout_name;

[fileout_name,fileout_path]=uiputfile('*.mat','Name of the netcdf file where to save formatted raw data:',filein_name(1:end-4));
cd(fileout_path);


time=datenum(matsave(:,2)+2000,matsave(:,3),matsave(:,4),matsave(:,5),matsave(:,6),matsave(:,7));
matsave(:,15)=time;
matsave=sortrows(matsave,15);

time=matsave(:,15);
Hs=matsave(:,9);
Tp=matsave(:,10);
Dp=matsave(:,11);
water_level=matsave(:,12);

H10=matsave(:,13);
T01=matsave(:,14);


adcp_wave.time.val=time;adcp_wave.time.obj='Matlab serial time';
adcp_wave.year.val=matsave(:,2)+2000;adcp_wave.year.obj='Year';
adcp_wave.month.val=matsave(:,3);adcp_wave.month.obj='Month';
adcp_wave.day.val=matsave(:,4);adcp_wave.day.obj='Day';
adcp_wave.hour.val=matsave(:,5);adcp_wave.hour.obj='Hour';
adcp_wave.minute.val=matsave(:,6);adcp_wave.minute.obj='Minute';
adcp_wave.second.val=matsave(:,7);adcp_wave.second.obj='Second';

adcp_wave.hs.val=Hs;adcp_wave.hs.obj='Significant wave height (m)';
adcp_wave.h10.val=H10;adcp_wave.h10.obj='Largest wave height (H10) (m)';
adcp_wave.tp.val=Tp; adcp_wave.tp.obj='Peak wave period (s)';
adcp_wave.t01.val=T01;adcp_wave.t01.obj='Mean wave period (s)';
adcp_wave.dp.val=Dp;adcp_wave.dp.obj='Peak wave direction (°)';
adcp_wave.depth.val=water_level;adcp_wave.depth.obj='Water depth (m)';

save(fileout_name,'adcp_wave')

cd(path_routines);


figure

subplot(3,1,1)
plot(time,water_level),xlabel('time'),ylabel('water depth h (m)')
subplot(3,1,2)
plot(time,Hs),xlabel('time'),ylabel('Significant wave height Hs (m)')
subplot(3,1,3)
plot(time,Tp),xlabel('time'),ylabel('Wave period Tp (s)')


