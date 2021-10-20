function save_log

global fileout_name
global path_routines
global matsave

filein_name=fileout_name;

[fileout_name,fileout_path]=uiputfile('*.nc','Name of the netcdf file where to save formatted raw data:',filein_name(1:end-4));
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


prompt={'Enter the netcdf header lines'};
name='netcdf header';
numlines=1;
defaultanswer={filein_name(1:end-4)};
answer=inputdlg(prompt,name,numlines,defaultanswer);

ncquiet
%filesave='adcp_frame_scope08_wave.nc';
ncadv=netcdf(fileout_name,'clobber');

ncadv.description=answer;
ncadv.author='Your Name - Ifremer';
ncadv.date=datestr(date);

ncadv('nbdata')=length(time);

ncadv{'time'}={'nbdata'}; ncadv{'time'}.units='days';ncadv{'time'}.label='Matlab serial date'; ncadv{'time'}(:)=time;
ncadv{'Year'}={'nbdata'}; ncadv{'Year'}.units='Year';ncadv{'Year'}.label='Year'; ncadv{'Year'}(:)=matsave(:,2)+2000;
ncadv{'Month'}={'nbdata'}; ncadv{'Month'}.units='Month';ncadv{'Month'}.label='Month'; ncadv{'Month'}(:)=matsave(:,3);
ncadv{'Day'}={'nbdata'}; ncadv{'Day'}.units='Day';ncadv{'Day'}.label='Day'; ncadv{'Day'}(:)=matsave(:,4);
ncadv{'Hour'}={'nbdata'}; ncadv{'Hour'}.units='Hour';ncadv{'Hour'}.label='Hour'; ncadv{'Hour'}(:)=matsave(:,5);
ncadv{'Minutes'}={'nbdata'}; ncadv{'Minutes'}.units='Minutes';ncadv{'Minutes'}.label='Minutes'; ncadv{'Minutes'}(:)=matsave(:,6);
ncadv{'Seconds'}={'nbdata'}; ncadv{'Seconds'}.units='Seconds';ncadv{'Seconds'}.label='Seconds'; ncadv{'Seconds'}(:)=matsave(:,7);


ncadv{'Hs'}={'nbdata'}; ncadv{'Hs'}.units='m';ncadv{'Hs'}.label='Significant wave height'; ncadv{'Hs'}(:)=Hs;
ncadv{'H10'}={'nbdata'}; ncadv{'H10'}.units='m';ncadv{'H10'}.label='maximum wave height'; ncadv{'H10'}(:)=H10;
ncadv{'Tp'}={'nbdata'}; ncadv{'Tp'}.units='s';ncadv{'Tp'}.label='Peak wave period'; ncadv{'Tp'}(:)=Tp;
ncadv{'T01'}={'nbdata'}; ncadv{'T01'}.units='s';ncadv{'T01'}.label='Mean wave period'; ncadv{'T01'}(:)=T01;
ncadv{'Dp'}={'nbdata'}; ncadv{'Dp'}.units='Deg';ncadv{'Dp'}.label='Peak wave direction'; ncadv{'Dp'}(:)=Dp;
ncadv{'h'}={'nbdata'}; ncadv{'h'}.units='m';ncadv{'h'}.label='Water level'; ncadv{'h'}(:)=water_level;

close(ncadv)
cd(path_routines);


figure

subplot(3,1,1)
plot(time,water_level),xlabel('time'),ylabel('water depth h (m)')
subplot(3,1,2)
plot(time,Hs),xlabel('time'),ylabel('Significant wave height Hs (m)')
subplot(3,1,3)
plot(time,Tp),xlabel('time'),ylabel('Wave period Tp (s)')


