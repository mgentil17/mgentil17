function concat_adcp_proc



veout=[];
vnout=[];
vuout=[];
spdout=[];
dirout=[];
a1out=[];
a2out=[];
a3out=[];
a4out=[];
ameanout=[];
c1out=[];
c2out=[];
c3out=[];
c4out=[];
cmeanout=[];
Tbout=[];
Sbout=[];
Battout=[];
pitchout=[];
headingout=[];
rollout=[];
dateout=[];
Dout=[];
BIout=[];
Mglout=[];

[filein_name,filein_path]=uigetfile('*.nc','Open the .nc file with the pre-processed data','MultiSelect','on');
cd(filein_path);
filetmp=filein_name;
filein_name(end)=filetmp(1);
for i=1:length(filein_name)-1
filein_name(i)=filetmp(i+1);
end

clear filetmp

for k=1:length(filein_name)

ncout=netcdf(filein_name{k});

theta=ncout{'theta'}(:);
WN=ncout{'WN'}(:);
WS=ncout{'WS'}(:);
WB=ncout{'WB'}(:);
IMM=ncout{'IMM'}(:);
ens_int=ncout{'ens_int'}(:);
f=ncout{'f'}(:);
a_t=ncout{'a_t'}(:);
ouv=ncout{'ouv'}(:);
Kc=ncout{'Kc'}(:);
EC0=ncout{'EC0'}(:);
B=ncout{'B'}(:);
SL0=ncout{'SL0'}(:);
Lt=ncout{'Lt'}(:);


binpos=ncout{'binpos'}(:);
yrange=ncout{'yrange'}(:);

veout=[veout; ncout{'ve'}(:,:)];
vnout=[vnout; ncout{'vn'}(:,:)];
vuout=[vuout; ncout{'vu'}(:,:)];


spdout=[spdout; ncout{'spd'}(:,:)];
dirout=[dirout; ncout{'dir'}(:,:)];

a1out=[a1out; ncout{'a1'}(:,:)];
a2out=[a2out; ncout{'a2'}(:,:)];
a3out=[a3out; ncout{'a3'}(:,:)];
a4out=[a4out; ncout{'a4'}(:,:)];
ameanout=[ameanout; ncout{'amean'}(:,:)];

c1out=[c1out; ncout{'c1'}(:,:)];
c2out=[c2out; ncout{'c2'}(:,:)];
c3out=[c3out; ncout{'c3'}(:,:)];
c4out=[c4out; ncout{'c4'}(:,:)];
cmeanout=[cmeanout; ncout{'cmean'}(:,:)];

BIout=[BIout; ncout{'BI'}(:,:)];
Mglout=[Mglout; ncout{'Mgl'}(:,:)];

Tbout=[Tbout; ncout{'Tb_adcp'}(:)];
Sbout=[Sbout; ncout{'Sb_adcp'}(:)];
Battout=[Battout; ncout{'Batt'}(:)];
Dout=[Dout; ncout{'Depth_adcp'}(:)];
headingout=[headingout; ncout{'heading'}(:)];
pitchout=[pitchout; ncout{'pitch'}(:)];
rollout=[rollout; ncout{'roll'}(:)];

dateout=[dateout; ncout{'time'}(:)];
close(ncout)
end


%size(vuout)
%size(vnout)
%WN
length(dateout)

[fileout_name,fileout_path]=uiputfile('*.nc','Name of the netcdf file where to save formatted raw data:');
cd(fileout_path)

ncquiet
fileout_name=strcat(fileout_name(1:end-3),'.nc');
ncout=netcdf(fileout_name,'clobber');

ncout.description='WADCP 0.3 ADP results';
ncout.author='Your Name';
%ncout.date=datestr(date);


ncout('param')=1;
ncout('nbens')=length(dateout);
ncout('WN')=WN;
ncout('yrange')=length(yrange);

ncout{'nbens'}={'param'};  ncout{'nbens'}.units='number of ensembles';ncout{'nbens'}(1)=length(dateout);
ncout{'theta'}={'param'};  ncout{'theta'}.units='beam angle (deg)';ncout{'theta'}(1)=theta;
ncout{'WN'}={'param'};  ncout{'WN'}.units='number of cells';ncout{'WN'}(1)=WN;
ncout{'WS'}={'param'};  ncout{'WS'}.units='cell size (m)';ncout{'WS'}(1)=WS;
ncout{'WB'}={'param'};  ncout{'WB'}.units='blank size (m)';ncout{'WB'}(1)=WB;
ncout{'IMM'}={'param'};  ncout{'IMM'}.units='height of the transducers from bed (m)';ncout{'IMM'}(1)=IMM;
ncout{'ens_int'}={'param'};  ncout{'ens_int'}.units='ensemble interval (s)';ncout{'ens_int'}(1)=ens_int;
ncout{'f'}={'param'};  ncout{'f'}.units='transducer frequency (kHz)';ncout{'f'}(1)=f;
ncout{'a_t'}={'param'};  ncout{'a_t'}.units='blank size (Deg)';ncout{'a_t'}(1)=a_t;
ncout{'ouv'}={'param'};  ncout{'ouv'}.units='blank size (Deg)';ncout{'ouv'}(1)=ouv;
ncout{'Kc'}={'param'};  ncout{'Kc'}.units='count/dB conversion coefficient (dB/count)';ncout{'Kc'}(1)=Kc;
ncout{'EC0'}={'param'};  ncout{'EC0'}.units='ADP internal noise';ncout{'EC0'}(1)=EC0;
ncout{'B'}={'param'};  ncout{'B'}.units='ADP noise in air (dB)';ncout{'B'}(1)=B;
ncout{'SL0'}={'param'};  ncout{'SL0'}.units='Emitted power signal (dB)';ncout{'SL0'}(1)=SL0;
ncout{'Lt'}={'param'};  ncout{'Lt'}.units='Transmit pulse length (m)';ncout{'Lt'}(1)=Lt;

ncout{'binpos'}={'WN'};  ncout{'binpos'}.units='Bin cell position from the bed (m)';ncout{'binpos'}(:)=binpos;
ncout{'yrange'}={'yrange'};  ncout{'yrange'}.units='Bin cell position from the bed (m) for BI and Mgl data';ncout{'yrange'}(:)=yrange;

ncout{'ve'}={'nbens','WN'};  ncout{'ve'}.units='Velocity (Eastward)';ncout{'ve'}(:,:)=veout;
ncout{'vn'}={'nbens','WN'};  ncout{'vn'}.units='Velocity (Northward)';ncout{'vn'}(:,:)=vnout;
ncout{'vu'}={'nbens','WN'};  ncout{'vu'}.units='Velocity (Upward)';ncout{'vu'}(:,:)=vuout;

ncout{'spd'}={'nbens','WN'};  ncout{'spd'}.units='Current Velocity';ncout{'spd'}(:,:)=spdout;
ncout{'dir'}={'nbens','WN'};  ncout{'dir'}.units='Current Direction';ncout{'dir'}(:,:)=dirout;

ncout{'a1'}={'nbens','WN'};  ncout{'a1'}.units='Backscattered signal beam 1 (count)';ncout{'a1'}(:,:)=a1out;
ncout{'a2'}={'nbens','WN'};  ncout{'a2'}.units='Backscattered signal beam 2 (count)';ncout{'a2'}(:,:)=a2out;
ncout{'a3'}={'nbens','WN'};  ncout{'a3'}.units='Backscattered signal beam 3 (count)';ncout{'a3'}(:,:)=a3out;
ncout{'a4'}={'nbens','WN'};  ncout{'a4'}.units='Backscattered signal beam 4 (count)';ncout{'a4'}(:,:)=a4out;
ncout{'amean'}={'nbens','WN'};  ncout{'amean'}.units='Backscattered signal (averaged) (count)';ncout{'amean'}(:,:)=ameanout;

ncout{'BI'}={'nbens','yrange'};  ncout{'BI'}.units='Backscatter Index (dB)';ncout{'BI'}(:,:)=BIout;
ncout{'Mgl'}={'nbens','yrange'};  ncout{'Mgl'}.units='SPM mass concentration (g/l)';ncout{'Mgl'}(:,:)=Mglout;

ncout{'c1'}={'nbens','WN'};  ncout{'c1'}.units='Correlation beam 1';ncout{'c1'}(:,:)=c1out;
ncout{'c2'}={'nbens','WN'};  ncout{'c2'}.units='Correlation beam 2';ncout{'c2'}(:,:)=c2out;
ncout{'c3'}={'nbens','WN'};  ncout{'c3'}.units='Correlation beam 3';ncout{'c3'}(:,:)=c3out;
ncout{'c4'}={'nbens','WN'};  ncout{'c4'}.units='Correlation beam 4';ncout{'c4'}(:,:)=c4out;
ncout{'cmean'}={'nbens','WN'};  ncout{'cmean'}.units='Correlation beam (averaged)';ncout{'cmean'}(:,:)=cmeanout;

ncout{'Tb_adcp'}={'nbens'};  ncout{'Tb_adcp'}.units='Bottom temperature (DegC)';ncout{'Tb_adcp'}(:)=Tbout;
ncout{'Sb_adcp'}={'nbens'};  ncout{'Sb_adcp'}.units='Bottom salinity (PSU)';ncout{'Sb_adcp'}(:)=Sbout;
ncout{'Batt'}={'nbens'};  ncout{'Batt'}.units='Batteries (V)';ncout{'Batt'}(:)=Battout;
ncout{'Depth_adcp'}={'nbens'};  ncout{'Depth_adcp'}.units='Pressure from ADCP sensor (m)';ncout{'Depth_adcp'}(:)=Dout;
ncout{'heading'}={'nbens'};  ncout{'heading'}.units='ADCP heading (Deg)';ncout{'heading'}(:)=headingout;
ncout{'pitch'}={'nbens'};  ncout{'pitch'}.units='ADCP pitch (Deg)';ncout{'pitch'}(:)=pitchout;
ncout{'roll'}={'nbens'};  ncout{'roll'}.units='ADCP roll (Deg)';ncout{'roll'}(:)=rollout;

ncout{'time'}={'nbens'};  ncout{'time'}.units='Matlab serial date and time for ensembles';ncout{'time'}(:)=dateout;
close(ncout)



