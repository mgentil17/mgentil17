function nc_concat

global time ve vn vu a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean spd dir
global path_routines filein_path ncin fileout_name fileout_path filein_name
global Batt Tb_adcp Depth_adcp Sb_adcp
global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 Lt binpos binpos_bt inst_mode

filein_name='PMD4_LCR-1_1A_AST140_1_proc.nc';
ncin=netcdf(filein_name);

time=ncin{'time'}(:);

theta=ncin{'theta'}(:);
WN=ncin{'WN'}(:);
WS=ncin{'WS'}(:);
WB=ncin{'WB'}(:);
IMM=ncin{'IMM'}(:);
ens_int=ncin{'ens_int'}(:);
f=ncin{'f'}(:);
a_t=ncin{'a_t'}(:);
ouv=ncin{'ouv'}(:);
Kc=ncin{'Kc'}(:);
EC0=ncin{'EC0'}(:);
B=ncin{'B'}(:);
SL0=ncin{'SL0'}(:);
Lt=ncin{'Lt'}(:);
inst_mode=ncin{'inst_mode'}(:);
binpos=ncin{'binpos'}(:);

yrange=ncin{'yrange'}(:);

ve=ncin{'ve'}(:,:);
vn=ncin{'vn'}(:,:);
vu=ncin{'vu'}(:,:);
spd=ncin{'spd'}(:,:);
dir=ncin{'dir'}(:,:);

a1=ncin{'a1'}(:,:);
a2=ncin{'a2'}(:,:);
a3=ncin{'a3'}(:,:);
a4=ncin{'a4'}(:,:);
amean=ncin{'amean'}(:,:);

c1=ncin{'c1'}(:,:);
c2=ncin{'c2'}(:,:);
c3=ncin{'c3'}(:,:);
c4=ncin{'c4'}(:,:);
cmean=ncin{'cmean'}(:,:);
Tb_adcp=ncin{'Tb_adcp'}(:);

BI=ncin{'BI'}(:,:);
Mgl=ncin{'Mgl'}(:,:);

Sb_adcp=ncin{'Sb_adcp'}(:);
Batt=ncin{'Batt'}(:);
Depth_adcp=ncin{'Depth_adcp'}(:);
heading=ncin{'heading'}(:);
pitch=ncin{'pitch'}(:);
roll=ncin{'roll'}(:);
close(ncin)

%%
%

filein_name='PMD4_LCR-1_1A_AST140_2_proc.nc';
ncin=netcdf(filein_name);

time=[time; ncin{'time'}(:)];

ve=[ve; ncin{'ve'}(:,:)];
vn=[vn; ncin{'vn'}(:,:)];
vu=[vu; ncin{'vu'}(:,:)];
spd=[spd; ncin{'spd'}(:,:)];
dir=[dir; ncin{'dir'}(:,:)];

a1=[a1; ncin{'a1'}(:,:)];
a2=[a2; ncin{'a2'}(:,:)];
a3=[a3; ncin{'a3'}(:,:)];
a4=[a4; ncin{'a4'}(:,:)];
amean=[amean; ncin{'amean'}(:,:)];

c1=[c1; ncin{'c1'}(:,:)];
c2=[c2; ncin{'c2'}(:,:)];
c3=[c3; ncin{'c3'}(:,:)];
c4=[c4; ncin{'c4'}(:,:)];
cmean=[cmean; ncin{'cmean'}(:,:)];
Tb_adcp=[Tb_adcp; ncin{'Tb_adcp'}(:)];

size(BI)

BI=[BI; ncin{'BI'}(:,:)];
Mgl=[Mgl; ncin{'Mgl'}(:,:)];

size(BI)

Sb_adcp=[Sb_adcp; ncin{'Sb_adcp'}(:)];
Batt=[Batt; ncin{'Batt'}(:)];
Depth_adcp=[Depth_adcp; ncin{'Depth_adcp'}(:)];
heading=[heading; ncin{'heading'}(:)];
pitch=[pitch; ncin{'pitch'}(:)];
roll=[roll; ncin{'roll'}(:)];
close(ncin)

filein_name='PMD4_LCR-1_1B_AST140_proc.nc';
ncin=netcdf(filein_name);

time=[time; ncin{'time'}(:)];

ve=[ve; ncin{'ve'}(:,:)];
vn=[vn; ncin{'vn'}(:,:)];
vu=[vu; ncin{'vu'}(:,:)];
spd=[spd; ncin{'spd'}(:,:)];
dir=[dir; ncin{'dir'}(:,:)];

a1=[a1; ncin{'a1'}(:,:)];
a2=[a2; ncin{'a2'}(:,:)];
a3=[a3; ncin{'a3'}(:,:)];
a4=[a4; ncin{'a4'}(:,:)];
amean=[amean; ncin{'amean'}(:,:)];

c1=[c1; ncin{'c1'}(:,:)];
c2=[c2; ncin{'c2'}(:,:)];
c3=[c3; ncin{'c3'}(:,:)];
c4=[c4; ncin{'c4'}(:,:)];
cmean=[cmean; ncin{'cmean'}(:,:)];
Tb_adcp=[Tb_adcp; ncin{'Tb_adcp'}(:)];

BI=[BI; ncin{'BI'}(:,:)];
Mgl=[Mgl; ncin{'Mgl'}(:,:)];

size(BI)

Sb_adcp=[Sb_adcp; ncin{'Sb_adcp'}(:)];
Batt=[Batt; ncin{'Batt'}(:)];
Depth_adcp=[Depth_adcp; ncin{'Depth_adcp'}(:)];
heading=[heading; ncin{'heading'}(:)];
pitch=[pitch; ncin{'pitch'}(:)];
roll=[roll; ncin{'roll'}(:)];
close(ncin)

nbens=length(time);

ncquiet
fileout_name='PMD4_le_croizic.nc';
ncout=netcdf(fileout_name,'clobber');

ncout.description='WADCP 0.3 ADP results';
ncout.author='Your Name';
ncout.date=datestr(date);


ncout('param')=1;
ncout('nbens')=nbens;
ncout('WN')=WN;
ncout('yrange')=length(yrange);

ncout{'nbens'}={'param'};  ncout{'nbens'}.units='number of ensembles';ncout{'nbens'}(1)=nbens;
ncout{'theta'}={'param'};  ncout{'theta'}.units='beam angle (deg)';ncout{'theta'}(1)=theta;
ncout{'WN'}={'param'};  ncout{'WN'}.units='number of cells';ncout{'WN'}(1)=WN;
ncout{'WS'}={'param'};  ncout{'WS'}.units='cell size (m)';ncout{'WS'}(1)=WS;
ncout{'WB'}={'param'};  ncout{'WB'}.units='blank size (m)';ncout{'WB'}(1)=WB;
ncout{'IMM'}={'param'};  ncout{'IMM'}.units='height of the transducers from bed (m)';ncout{'IMM'}(1)=IMM;
ncout{'ens_int'}={'param'};  ncout{'ens_int'}.units='ensemble interval (s)';ncout{'ens_int'}(1)=ens_int;
ncout{'f'}={'param'};  ncout{'f'}.units='transducer frequency (kHz)';ncout{'f'}(1)=f;
ncout{'a_t'}={'param'};  ncout{'a_t'}.units='Diameter of transducer (m)';ncout{'a_t'}(1)=a_t;
ncout{'ouv'}={'param'};  ncout{'ouv'}.units='beam spreading (Deg)';ncout{'ouv'}(1)=ouv;
ncout{'Kc'}={'param'};  ncout{'Kc'}.units='count/dB conversion coefficient (dB/count)';ncout{'Kc'}(1)=Kc;
ncout{'EC0'}={'param'};  ncout{'EC0'}.units='transducer internal noise (count)';ncout{'EC0'}(1)=EC0;
ncout{'B'}={'param'};  ncout{'B'}.units='transducer internal noise (dB)';ncout{'B'}(1)=B;
ncout{'SL0'}={'param'};  ncout{'SL0'}.units='Emitted power signal level (dB)';ncout{'SL0'}(1)=SL0;
ncout{'Lt'}={'param'};  ncout{'Lt'}.units='Transmit pulse length (m)';ncout{'Lt'}(1)=Lt;

ncout{'binpos'}={'WN'};  ncout{'binpos'}.units='Bin cell position from the bed (m)';ncout{'binpos'}(:)=binpos;
ncout{'yrange'}={'yrange'};  ncout{'yrange'}.units='Bin cell position from the bed for BI and Mgl data (m)';ncout{'yrange'}(:)=yrange;




ncout{'ve'}={'nbens','WN'};  ncout{'ve'}.units='Velocity (Eastward)';ncout{'ve'}(:,:)=ve;
ncout{'vn'}={'nbens','WN'};  ncout{'vn'}.units='Velocity (Northward)';ncout{'vn'}(:,:)=vn;
ncout{'vu'}={'nbens','WN'};  ncout{'vu'}.units='Velocity (Upward)';ncout{'vu'}(:,:)=vu;

ncout{'spd'}={'nbens','WN'};  ncout{'spd'}.units='Current Velocity';ncout{'spd'}(:,:)=spd;
ncout{'dir'}={'nbens','WN'};  ncout{'dir'}.units='Current Direction';ncout{'dir'}(:,:)=dir;

ncout{'a1'}={'nbens','WN'};  ncout{'a1'}.units='Backscattered signal beam 1 (count)';ncout{'a1'}(:,:)=a1;
ncout{'a2'}={'nbens','WN'};  ncout{'a2'}.units='Backscattered signal beam 2 (count)';ncout{'a2'}(:,:)=a2;
ncout{'a3'}={'nbens','WN'};  ncout{'a3'}.units='Backscattered signal beam 3 (count)';ncout{'a3'}(:,:)=a3;
ncout{'a4'}={'nbens','WN'};  ncout{'a4'}.units='Backscattered signal beam 4 (count)';ncout{'a4'}(:,:)=a4;
ncout{'amean'}={'nbens','WN'};  ncout{'amean'}.units='Backscattered signal (averaged) (count)';ncout{'amean'}(:,:)=amean;

ncout{'c1'}={'nbens','WN'};  ncout{'c1'}.units='Correlation beam 1';ncout{'c1'}(:,:)=c1;
ncout{'c2'}={'nbens','WN'};  ncout{'c2'}.units='Correlation beam 2';ncout{'c2'}(:,:)=c2;
ncout{'c3'}={'nbens','WN'};  ncout{'c3'}.units='Correlation beam 3';ncout{'c3'}(:,:)=c3;
ncout{'c4'}={'nbens','WN'};  ncout{'c4'}.units='Correlation beam 4';ncout{'c4'}(:,:)=c4;
ncout{'cmean'}={'nbens','WN'};  ncout{'cmean'}.units='Correlation beam (averaged)';ncout{'cmean'}(:,:)=cmean;


ncout{'BI'}={'nbens','yrange'};  ncout{'BI'}.units='Backscatter Index (dB)';ncout{'BI'}(:,:)=BI;
ncout{'Mgl'}={'nbens','yrange'};  ncout{'Mgl'}.units='SPM mass concentration (g/l)';ncout{'Mgl'}(:,:)=Mgl;

ncout{'Tb_adcp'}={'nbens'};  ncout{'Tb_adcp'}.units='Bottom temperature (DegC)';ncout{'Tb_adcp'}(:)=Tb_adcp;
ncout{'Sb_adcp'}={'nbens'};  ncout{'Sb_adcp'}.units='Bottom salinity (PSU)';ncout{'Sb_adcp'}(:)=Sb_adcp;
ncout{'Batt'}={'nbens'};  ncout{'Batt'}.units='Batteries (Count)';ncout{'Batt'}(:)=Batt;
ncout{'Depth_adcp'}={'nbens'};  ncout{'Depth_adcp'}.units='Pressure from ADCP sensor (m)';ncout{'Depth_adcp'}(:)=Depth_adcp;
ncout{'heading'}={'nbens'};  ncout{'heading'}.units='ADCP heading (Deg)';ncout{'heading'}(:)=heading;
ncout{'pitch'}={'nbens'};  ncout{'pitch'}.units='ADCP pitch (Deg)';ncout{'pitch'}(:)=pitch;
ncout{'roll'}={'nbens'};  ncout{'roll'}.units='ADCP roll (Deg)';ncout{'roll'}(:)=roll;

ncout{'time'}={'nbens'};  ncout{'time'}.units='Matlab serial time for ensembles';ncout{'time'}(:)=time;
close(ncout)