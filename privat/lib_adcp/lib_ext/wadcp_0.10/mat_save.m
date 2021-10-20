function mat_save

global time ve vn vu dir spd a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean
global path_routines filein_path fileout_name fileout_path filein_name filein
global Batt Tb_adcp Depth_adcp Sb_adcp
global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 Lt binpos inst_mode binpos_bt inst_type
global mode_exp

cd(filein_path);

display('enter save routine')


if mode_exp==0
[fileout_name,fileout_path]=uiputfile('*.mat','Name of the netcdf file where to save formatted raw data:',filein_name(1:end-4));
cd(fileout_path)
end
fileout_name=strcat(filein,'_data.mat');


data.description='WADCP 0.5 ADP results'; data.date=datestr(date);


data.nbens.obj='Number of ensembles'; data.nbens.unit='NA'; data.nbens.val=nbens;
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
data.inst_type.obj='Instrument type: Nortek, Sontek, RDI defined as string'; data.inst_type.unit='NA'; data.inst_type.val=inst_type;

if inst_mode==1
    data.binpos.obj='Bin cell position from the sensor'; data.binpos.unit='m'; data.binpos.val=binpos;
    data.binpos_bt.obj='Bin cell position from the bed for BT mode'; data.binpos_bt.unit='m'; data.binpos_bt.val=binpos_bt;

else
    data.binpos.obj='Bin cell position from the sensor'; data.binpos.unit='m'; data.binpos.val=binpos;
end

data.ve.obj='Velocity (Eastward)'; data.ve.unit='m/s'; data.ve.val=ve;
data.vn.obj='Velocity (Northward)'; data.vn.unit='m/s'; data.vn.val=vn;
data.vu.obj='Velocity (Upward)'; data.vu.unit='m/s'; data.vu.val=vu;

data.spd.obj='Current Velocity'; data.spd.unit='m/s'; data.spd.val=spd;
data.dir.obj='Current Direction'; data.dir.unit='degrees 0 East anticlockwise'; data.dir.val=dir;

data.a1.obj='Backscattered signal beam 1'; data.a1.unit='count'; data.a1.val=a1;
data.a2.obj='Backscattered signal beam 2'; data.a2.unit='count'; data.a2.val=a2;
data.a3.obj='Backscattered signal beam 3'; data.a3.unit='count'; data.a3.val=a3;
data.a4.obj='Backscattered signal beam 4'; data.a4.unit='count'; data.a4.val=a4;
data.amean.obj='Backscattered signal (averaged)'; data.amean.unit='count'; data.amean.val=amean;

data.c1.obj='Correlation beam 1'; data.c1.unit='NA';data.c1.val=c1;
data.c2.obj='Correlation beam 2'; data.c2.unit='NA';data.c2.val=c2;
data.c3.obj='Correlation beam 3'; data.c3.unit='NA';data.c3.val=c3;
data.c4.obj='Correlation beam 4'; data.c4.unit='NA';data.c4.val=c4;
data.cmean.obj='Correlation beam (averaged)'; data.cmean.unit='NA';data.cmean.val=cmean;

data.Tb_adcp.obj='Bottom temperature'; data.Tb_adcp.unit='Celcius Degrees'; data.Tb_adcp.val=Tb_adcp;
data.Sb_adcp.obj='Bottom salinity'; data.Sb_adcp.unit='PSU'; data.Sb_adcp.val=Sb_adcp;
data.Batt.obj='Batteries power'; data.Batt.unit='Count'; data.Batt.val=Batt;
data.Depth_adcp.obj='Pressure from ADCP sensor'; data.Depth_adcp.unit='m'; data.Depth_adcp.val=Depth_adcp;
data.heading.obj= 'ADCP heading'; data.heading.unit='Degrees'; data.heading.val=heading;
data.pitch.obj='ADCP pitch'; data.pitch.unit='Degrees'; data.pitch.val=pitch;
data.roll.obj='ADCP roll'; data.roll.unit='Degrees'; data.roll.val=roll;

data.time.obj='Matlab serial time for ensembles'; data.time.unit='Days since 1900/1/1'; data.time.val=time;


save(fileout_name,'data')
data=[];
clear data
cd(path_routines)
