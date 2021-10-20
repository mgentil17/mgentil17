function nc_load

global time ve vn vu a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean spd dir
global path_routines filein_path ncin fileout_name fileout_path filein_name
global Batt Tb_adcp Depth_adcp Sb_adcp
global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 Lt binpos binpos_bt inst_mode mode_exp

if isnan(fileout_name)
    

cd(filein_path)
ncin=netcdf(filein_name);
else
    cd(fileout_path)
    ncin=netcdf(fileout_name);
    
end
cd(path_routines)


time=ncin{'time'}(:);

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
time=time(timevec);

nbens=length(time);

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
binpos_bt=ncin{'binpos_bt'}(:);
binpos=ncin{'binpos'}(:);



ve=ncin{'ve'}(:,:);
ve=ve(timevec,:);
vn=ncin{'vn'}(:,:);
vn=vn(timevec,:);
vu=ncin{'vu'}(:,:);
vu=vu(timevec,:);
spd=ncin{'spd'}(:,:);
spd=spd(timevec,:);
dir=ncin{'dir'}(:,:);
dir=dir(timevec,:);

a1=ncin{'a1'}(:,:);
a1=a1(timevec,:);
a2=ncin{'a2'}(:,:);
a2=a2(timevec,:);
a3=ncin{'a3'}(:,:);
a3=a3(timevec,:);
a4=ncin{'a4'}(:,:);
a4=a4(timevec,:);
amean=ncin{'amean'}(:,:);
amean=amean(timevec,:);

c1=ncin{'c1'}(:,:);
c1=c1(timevec,:);
c2=ncin{'c2'}(:,:);
c2=c2(timevec,:);
c3=ncin{'c3'}(:,:);
c3=c3(timevec,:);
c4=ncin{'c4'}(:,:);
c4=c4(timevec,:);
cmean=ncin{'cmean'}(:,:);
cmean=cmean(timevec,:);

Tb_adcp=ncin{'Tb_adcp'}(:);
Tb_adcp=Tb_adcp(timevec,:);

Sb_adcp=ncin{'Sb_adcp'}(:);
Sb_adcp=Sb_adcp(timevec,:);
Batt=ncin{'Batt'}(:);
Batt=Batt(timevec,:);
Depth_adcp=ncin{'Depth_adcp'}(:)+0.625; % cf document IXsurvey, capteur de pression à 32.5cm du fond
Depth_adcp=Depth_adcp(timevec,:);
heading=ncin{'heading'}(:);
heading=heading(timevec,:);
pitch=ncin{'pitch'}(:);
pitch=pitch(timevec,:);
roll=ncin{'roll'}(:);
roll=roll(timevec,:);
close(ncin)
cd(path_routines)

