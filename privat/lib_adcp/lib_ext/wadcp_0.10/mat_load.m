function mat_load

global time ve vn vu a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean spd dir
global path_routines filein_path ncin fileout_name fileout_path filein_name
global Batt Tb_adcp Depth_adcp Sb_adcp 
global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 Lt binpos binpos_bt inst_mode inst_type mode_exp postproc_choice tminproc tmaxproc

if isnan(fileout_name)
 cd(filein_path)
    load(filein_name);


else
    cd(fileout_path)
    load(fileout_name);
    
end
cd(path_routines)


time=data.time.val;

inst_type=data.inst_type.val;


t_min=min(time);
t_max=max(time);

if mode_exp==0
prompt={'tmin', 'tmax'};
dlg_title='Data processing - time limits';
num_lines=1;
def={datestr(t_min),datestr(t_max)};
input1=inputdlg(prompt,dlg_title,num_lines,def);
tinf=datenum(input1(1));
tsup=datenum(input1(2));

else
   
    if postproc_choice==1
        tsup=t_max;
        tinf=t_min;
        
    else
        
        tinf=datenum(tminproc);
        tsup=datenum(tmaxproc);
    end
        
    
end


timevec=find(time>tinf & time<tsup);
timevec=timevec(1:end-1);

time2=time(timevec);

% timediff=diff(time2);
% findtimediff=find(abs(timediff-mean(timediff))>0.01);
% 
% if ~isempty(findtimediff)
% 
% timevec2=[];
% i1=1;
% i2=1;
% 
% for i=1:length(findtimediff);
%     i2=findtimediff(i)-2;
%  
%     timevec2=[timevec2; timevec(i1:i2)];
%     i1=findtimediff(i)+1;
% end
%     
% timevec=timevec2;
% end
time=time(timevec);


nbens=length(time);

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

ve=data.ve.val;
ve=ve(timevec,:);
vn=data.vn.val;
vn=vn(timevec,:);
vu=data.vu.val;
vu=vu(timevec,:);
spd=data.spd.val;
spd=spd(timevec,:);
dir=data.dir.val;
dir=dir(timevec,:);

a1=data.a1.val;
a1=a1(timevec,:);
a2=data.a2.val;
a2=a2(timevec,:);
a3=data.a3.val;
a3=a3(timevec,:);
a4=data.a4.val;
a4=a4(timevec,:);

amean=data.amean.val;
amean=amean(timevec,:);


c1=data.c1.val;
c1=c1(timevec,:);
c2=data.c2.val;
c2=c2(timevec,:);
c3=data.c3.val;
c3=c3(timevec,:);
c4=data.c4.val;
c4=c4(timevec,:);
cmean=data.cmean.val;
cmean=cmean(timevec,:);

Tb_adcp=data.Tb_adcp.val;
Tb_adcp=Tb_adcp(timevec);

Sb_adcp=data.Sb_adcp.val;
Sb_adcp=Sb_adcp(timevec);
% Batt=data.Batt.val;
% Batt=Batt(timevec);
Batt=ones(size(Tb_adcp));
Depth_adcp=data.Depth_adcp.val;%+0.625; % cf document IXsurvey, capteur de pression à 62.5cm du fond
Depth_adcp=Depth_adcp(timevec);
heading=data.heading.val;
heading=heading(timevec);
pitch=data.pitch.val;
pitch=pitch(timevec);
roll=data.roll.val;
roll=roll(timevec);
data=[];
clear data
cd(path_routines)

