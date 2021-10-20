function concat_adcp_mat



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
timeout=[];
Dout=[];

path_routines=cd;

[filein_name,filein_path]=uigetfile('*.mat','Open the .mat file with the pre-processed data','MultiSelect','on');
cd(filein_path);

filetmp=filein_name;
filein_name(end)=filetmp(1);
for i=1:length(filein_name)-1
filein_name(i)=filetmp(i+1);
end



clear filetmp

for k=1:length(filein_name)

    load(filein_name{k});
    
    timeout=[timeout; data.time.val];
    
       
theta=data.theta.val;
WN=data.WN.val;
WS=data.WS.val;
WB=data.WB.val;
IMM=data.IMM.val;
ens_int=data.ens_int.val;
f=data.f.val;
a_t=data.a_t.val;
ouv=data.ouv.val;
Kc=data.Kc.val;
EC0=data.EC0.val;
B=data.B.val;
SL0=data.SL0.val;
Lt=data.Lt.val;
inst_mode=data.inst_mode.val;
if inst_mode==1
binpos_bt=data.binpos_bt.val;
end
binpos=data.binpos.val;

veout=[veout; data.ve.val];
vnout=[vnout; data.vn.val];
vuout=[vuout; data.vu.val];
spdout=[spdout; data.spd.val];
dirout=[dirout; data.dir.val];

a1out=[a1out; data.a1.val];
a2out=[a2out; data.a2.val];
a3out=[a3out; data.a3.val];
a4out=[a4out; data.a4.val];
ameanout=[ameanout; data.amean.val];

c1out=[c1out; data.c1.val];
c2out=[c2out; data.c2.val];
c3out=[c3out; data.c3.val];
c4out=[c4out; data.c4.val];
cmeanout=[cmeanout; data.cmean.val];

Tbout=[Tbout; data.Tb_adcp.val];
Sbout=[Sbout; data.Sb_adcp.val];
Battout=[Battout; data.Batt.val];
Dout=[Dout; data.Depth_adcp.val];%+0.625; % cf document IXsurvey, capteur de pression à 32.5cm du fond
headingout=[headingout; data.heading.val];
pitchout=[pitchout; data.pitch.val];
rollout=[rollout; data.roll.val];
    
    
end


[fileout_name,fileout_path]=uiputfile('*.mat','Name of the netcdf file where to save formatted raw data:');
cd(fileout_path)

data.description='WADCP 0.5 ADP results'; data.date=datestr(date);

data.nbens.obj='Number of ensembles'; data.nbens.unit='NA'; data.nbens.val=length(timeout);
data.WN.obj='Number of cells'; data.WN.unit='NA'; data.WN.val=WN;
data.theta.obj='Beam angle';data.theta.unit='Degrees / N';data.theta.val=theta;
data.WS.obj='Cell size'; data.WS.unit='m'; data.WS.val=WS;
data.WB.obj='Blank size'; data.WB.unit='m'; data.WB.val=WB;
data.IMM.obj='Height of the transducers from bed'; data.IMM.unit='m'; data.IMM.val=IMM;
data.ens_int.obj='Ensemble interval'; data.ens_int.unit='s'; data.ens_int.val=ens_int;
data.f.obj='Transducer frequency'; data.f.unit='kHz'; data.f.val=f;
data.a_t.obj='Diameter of transducer '; data.a_t.unit='m'; data.a_t.val=a_t;
data.ouv.obj='Beam spreading'; data.ouv.unit='Degrees'; data.ouv.val=ouv;
data.Kc.obj='Count/dB conversion coefficient'; data.Kc.unit='dB/count'; data.Kc.val=Kc;
data.EC0.obj='Transducer internal noise'; data.EC0.unit='count'; data.EC0.val=EC0;
data.B.obj='transducer internal noise'; data.B.unit='dB'; data.B.val=B;
data.SL0.obj='Emitted power signal level'; data.SL0.obj='dB'; data.SL0.val=SL0;
data.Lt.obj='Transmit pulse length'; data.Lt.unit='m'; data.Lt.val=Lt;
data.inst_mode.obj='Instrument mode. 0: moored, 1: survey'; data.inst_mode.unit='NA'; data.inst_mode.val=inst_mode;


if inst_mode==1
    data.binpos.obj='Bin cell position from the sensor'; data.binpos.unit='m'; data.binpos.val=binpos;
    data.binpos_bt.obj='Bin cell position from the bed for BT mode'; data.binpos_bt.unit='m'; data.binpos_bt.val=binpos_bt;

else
    data.binpos.obj='Bin cell position from the sensor'; data.binpos.unit='m'; data.binpos.val=binpos;
end

data.ve.obj='Velocity (Eastward)'; data.ve.unit='m/s'; data.ve.val=veout;
data.vn.obj='Velocity (Northward)'; data.vn.unit='m/s'; data.vn.val=vnout;
data.vu.obj='Velocity (Upward)'; data.vu.unit='m/s'; data.vu.val=vuout;

data.spd.obj='Current Velocity'; data.spd.unit='m/s'; data.spd.val=spdout;
data.dir.obj='Current Direction'; data.dir.unit='degrees 0 East anticlockwise'; data.dir.val=dirout;

data.a1.obj='Backscattered signal beam 1'; data.a1.unit='count'; data.a1.val=a1out;
data.a2.obj='Backscattered signal beam 2'; data.a2.unit='count'; data.a2.val=a2out;
data.a3.obj='Backscattered signal beam 3'; data.a3.unit='count'; data.a3.val=a3out;
data.a4.obj='Backscattered signal beam 4'; data.a4.unit='count'; data.a4.val=a4out;
data.amean.obj='Backscattered signal (averaged)'; data.amean.unit='count'; data.amean.val=ameanout;

data.c1.obj='Correlation beam 1'; data.c1.unit='NA';data.c1.val=c1out;
data.c2.obj='Correlation beam 2'; data.c2.unit='NA';data.c2.val=c2out;
data.c3.obj='Correlation beam 3'; data.c3.unit='NA';data.c3.val=c3out;
data.c4.obj='Correlation beam 4'; data.c4.unit='NA';data.c4.val=c4out;
data.cmean.obj='Correlation beam (averaged)'; data.cmean.unit='NA';data.cmean.val=cmeanout;

data.Tb_adcp.obj='Bottom temperature'; data.Tb_adcp.unit='Celcius Degrees'; data.Tb_adcp.val=Tbout;
data.Sb_adcp.obj='Bottom salinity'; data.Sb_adcp.unit='PSU'; data.Sb_adcp.val=Sbout;
data.Batt.obj='Batteries power'; data.Batt.unit='Count'; data.Batt.val=Battout;
data.Depth_adcp.obj='Pressure from ADCP sensor'; data.Depth_adcp.unit='m'; data.Depth_adcp.val=Dout;
data.heading.obj= 'ADCP heading'; data.heading.unit='Degrees'; data.heading.val=headingout;
data.pitch.obj='ADCP pitch'; data.pitch.unit='Degrees'; data.pitch.val=pitchout;
data.roll.obj='ADCP roll'; data.roll.unit='Degrees'; data.roll.val=rollout;

data.time.obj='Matlab serial time for ensembles'; data.time.unit='Days since 1900/1/1'; data.time.val=timeout;

save(fileout_name,'data')

cd(path_routines)

