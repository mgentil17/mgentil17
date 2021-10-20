function nc_save

global time ve vn vu dir spd a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean
global path_routines filein_path fileout_name fileout_path filein_name filein
global Batt Tb_adcp Depth_adcp Sb_adcp
global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 Lt binpos inst_mode binpos_bt mode_exp

cd(filein_path);
if mode_exp==0
[fileout_name,fileout_path]=uiputfile('*.nc','Name of the netcdf file where to save formatted raw data:',filein_name(1:end-4));
cd(fileout_path)
end

ncquiet
fileout_name=strcat(filein,'.nc');
ncout=netcdf(fileout_name,'clobber');

ncout.description='WADCP 0.3 ADP results';
ncout.author='Your Name';
ncout.date=datestr(date);


ncout('param')=1;
ncout('nbens')=nbens;
ncout('WN')=WN;

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
ncout{'inst_mode'}={'param'};  ncout{'inst_mode'}.units='Instrument mode 0: moored, 1: survey';ncout{'inst_mode'}(1)=inst_mode;

if inst_mode==1
ncout{'binpos'}={'WN'};  ncout{'binpos'}.units='Bin cell position from the sensor (m)';ncout{'binpos'}(:)=binpos;    
ncout{'binpos_bt'}={'WN'};  ncout{'binpos_bt'}.units='Bin cell position from the bed for BT (m)';ncout{'binpos_bt'}(:)=binpos_bt;
else
ncout{'binpos'}={'WN'};  ncout{'binpos'}.units='Bin cell position from the bed (m)';ncout{'binpos'}(:)=binpos;
end

ncout{'ve'}={'nbens','WN'};  ncout{'ve'}.units='Velocity (Eastward) m/s';ncout{'ve'}(:,:)=ve;
ncout{'vn'}={'nbens','WN'};  ncout{'vn'}.units='Velocity (Northward) m/s';ncout{'vn'}(:,:)=vn;
ncout{'vu'}={'nbens','WN'};  ncout{'vu'}.units='Velocity (Upward) m/s';ncout{'vu'}(:,:)=vu;

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

ncout{'Tb_adcp'}={'nbens'};  ncout{'Tb_adcp'}.units='Bottom temperature (DegC)';ncout{'Tb_adcp'}(:)=Tb_adcp;
ncout{'Sb_adcp'}={'nbens'};  ncout{'Sb_adcp'}.units='Bottom salinity (PSU)';ncout{'Sb_adcp'}(:)=Sb_adcp;
ncout{'Batt'}={'nbens'};  ncout{'Batt'}.units='Batteries (Count)';ncout{'Batt'}(:)=Batt;
ncout{'Depth_adcp'}={'nbens'};  ncout{'Depth_adcp'}.units='Pressure from ADCP sensor (m)';ncout{'Depth_adcp'}(:)=Depth_adcp;
ncout{'heading'}={'nbens'};  ncout{'heading'}.units='ADCP heading (Deg)';ncout{'heading'}(:)=heading;
ncout{'pitch'}={'nbens'};  ncout{'pitch'}.units='ADCP pitch (Deg)';ncout{'pitch'}(:)=pitch;
ncout{'roll'}={'nbens'};  ncout{'roll'}.units='ADCP roll (Deg)';ncout{'roll'}(:)=roll;

ncout{'time'}={'nbens'};  ncout{'time'}.units='Matlab serial time for ensembles';ncout{'time'}(:)=time;
close(ncout)
cd(path_routines)
