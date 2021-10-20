function open_rdi_binary

global time ve vn vu dir spd a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean
global path_routines filein_path fileout_name fileout_path filein_name
global Batt Tb_adcp Depth_adcp Sb_adcp
global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 Lt binpos inst_mode binpos_bt inst_type mode_exp
if mode_exp==0
filetemp=[filein_path,filein_name];
else
    filetemp=[filein_path,filesep,filein_name];
end

VAcc=LitFichierWH(filetemp,0,0);
%VAcc=LitFichierWH('C:\Users\rverney\Desktop\TP_matlab_2012\WAVES043_000_000_CUR.PD0',0,0);
fid=fopen(filetemp,'rb');fseek(fid,0,'eof');taille=ftell(fid);fclose(fid);
%fid=fopen('C:\Users\rverney\Desktop\TP_matlab_2012\WAVES043_000_000_CUR.PD0','rb');fseek(fid,0,'eof');taille=ftell(fid);fclose(fid);
NbPings=taille/double(VAcc(1).VariableData.NBytes+2);

V=LitFichierWH(filetemp,0,NbPings);
%V=LitFichierWH('C:\Users\rverney\Desktop\TP_matlab_2012\WAVES043_000_000_CUR.PD0',0,NbPings);
time=datenum([V(3).VariableData.Year' V(3).VariableData.Month' V(3).VariableData.Day' V(3).VariableData.Hour' V(3).VariableData.Min' V(3).VariableData.Sec']);
nbens=length(time);
WN=V(2).VariableData.NCells(1);
theta=V(2).VariableData.BeamAngle(1);
WS=V(2).VariableData.CellDepth(1);
WB=V(2).VariableData.BlankDepth(1);

ens_int=diff(time);
ens_int=ens_int(1)*86400;


heading=(V(3).VariableData.Heading)';
pitch=(V(3).VariableData.Pitch)';
roll=(V(3).VariableData.Roll)';
Sb_adcp=(V(3).VariableData.Salinity)';
Tb_adcp=(V(3).VariableData.Temp)';


Depth_adcp=(V(3).VariableData.TransducerDepth)';

ve=V(4).VariableData(1:4:end);
vn=V(4).VariableData(2:4:end);
vu=V(4).VariableData(3:4:end);
[dir,spd]=cart2pol(ve,vn);
dir=mod(90-dir*180/pi,360);

ve=(reshape(ve,[],nbens))';
vu=(reshape(vu,[],nbens))';
vn=(reshape(vn,[],nbens))';
spd=(reshape(spd,[],nbens))';
dir=(reshape(dir,[],nbens))';

c1=V(5).VariableData(1:4:end);
c2=V(5).VariableData(2:4:end);
c3=V(5).VariableData(3:4:end);
c4=V(5).VariableData(4:4:end);

%%
clength=length(c1);
%%


cmean=(c1+c2+c3+c4)/4;

c1=(reshape(c1,[],nbens))';
c2=(reshape(c2,[],nbens))';
c3=(reshape(c3,[],nbens))';
c4=(reshape(c4,[],nbens))';

cmean=(reshape(cmean,[],nbens))';
a1=V(6).VariableData(1:4:end)/0.45;
a2=V(6).VariableData(2:4:end)/0.45;
a3=V(6).VariableData(3:4:end)/0.45;
a4=V(6).VariableData(4:4:end)/0.45;
%%
%a1=a1(end-clength+1:end);
% a2=a2(end-clength+1:end);
% a3=a3(end-clength+1:end);
% a4=a4(end-clength+1:end);
% %%
amean=(a1+a2+a3+a4)/4;

a1=(reshape(a1,[],nbens))';

% a1b=(reshape(a1,25,[]));
% size(c1)
% size(a1b)
% pause

a2=(reshape(a2,[],nbens))';
a3=(reshape(a3,[],nbens))';
a4=(reshape(a4,[],nbens))';
amean=(reshape(amean,[],nbens))';

inst_freq=num2str(V(2).VariableData.Frequency(1));
if mode_exp==0
%% RDI ADCP
% * deployment characteristics

prompt = {'Enter the distance from the bed - IMM (m)'};
dlg_title = 'Deployment characteristics';
num_lines = 1;
def = {'0.5'};
tmp = inputdlg(prompt,dlg_title,num_lines,def);
IMM=eval(char(tmp{1}));
end

R1=V(2).VariableData.Bin1Dstnc(1);  % position milieu premiere cellule courant [m]
Lt=V(2).VariableData.XmtLength(1);                   % transmit pulse [m]



binpos=R1+WS*(0:WN-1)+IMM;
if inst_mode==1
binpos_bt=-binpos+binpos(end);
end

if mode_exp==0

if strcmp(inst_freq,'1200')

f=1228.8;           % frequence de l'ADCP [kHz]
a_t=0.051;         % diametre du transducteur [m] (exemple de valeurs par défaut : RDI 79cm pour 300kHz 51cm pr 1200kHz)
theta=20;       % orientation du faisceau [deg]
ouv=0.99;         % ouverture angulaire du faisceau en Emission+Reception [deg] (exemple de valeurs par défaut : RDI 2.87deg pr 300kHz; 0.99deg pr 1200kHz)

Kc=0.42;    % relation dB/counts (valeur par defaut autour de 0.43)
EC0=45;     % niveau de bruit interne [counts] (valeur par defaut de 40 a 65 counts, determinable dans air)
B=70;       % niveau de bruit [dB/1microPa] (valeur par defaut de 70 a 96)
SL0=217;    %niveau emis [dB/1microPa] (valeur par defaut RDI : 216 pr 300kHz; 217 pr 1200kHz)

elseif strcmp(inst_freq,'600')
    
f=614.4;           % frequence de l'ADCP [kHz]
a_t=0.0762;         % diametre du transducteur [m] (exemple de valeurs par défaut : RDI 79cm pour 300kHz 51cm pr 1200kHz)
theta=20;       % orientation du faisceau [deg]
ouv=3.0;         % ouverture angulaire du faisceau en Emission+Reception [deg] (exemple de valeurs par défaut : RDI 2.87deg pr 300kHz; 0.99deg pr 1200kHz)

Kc=0.42;    % relation dB/counts (valeur par defaut autour de 0.43)
EC0=45;     % niveau de bruit interne [counts] (valeur par defaut de 40 a 65 counts, determinable dans air)
B=70;       % niveau de bruit [dB/1microPa] (valeur par defaut de 70 a 96)
SL0=217;    %niveau emis [dB/1microPa] (valeur par defaut RDI : 216 pr 300kHz; 217 pr 1200kHz)

elseif strcmp(inst_freq,'300')
    
f=300;           % frequence de l'ADCP [kHz]
a_t=0.079;         % diametre du transducteur [m]
theta=20;       % orientation du faisceau [deg]
ouv=2.87;         % ouverture angulaire du faisceau en Emission+Reception [deg]

Kc=0.42;    % relation dB/counts (valeur par defaut autour de 0.43)
EC0=45;     % niveau de bruit interne [counts] (valeur par defaut de 40 a 65 counts, determinable dans air)
B=70;       % niveau de bruit [dB/1microPa] (valeur par defaut de 70 a 96)
SL0=217;    %niveau emis [dB/1microPa] (valeur par defaut RDI : 216 pr 300kHz; 217 pr 1200kHz)

elseif strcmp(inst_freq,'150')
    
f=150;           % frequence de l'ADCP [kHz]
a_t=0.079;         % diametre du transducteur [m]
theta=20;       % orientation du faisceau [deg]
ouv=2.87;         % ouverture angulaire du faisceau en Emission+Reception [deg]

Kc=0.42;    % relation dB/counts (valeur par defaut autour de 0.43)
EC0=45;     % niveau de bruit interne [counts] (valeur par defaut de 40 a 65 counts, determinable dans air)
B=70;       % niveau de bruit [dB/1microPa] (valeur par defaut de 70 a 96)
SL0=217;    %niveau emis [dB/1microPa] (valeur par defaut RDI : 216 pr 300kHz; 217 pr 1200kHz)
end


prompt = {'Kc : dB/counts conversion coefficient (ranging from 0.4 and 0.45)';...
'EC0 : Internal noise [count] i.e. Echo value in air (ranging from 40 to 65)';...
'B : noise [dB/1microPa] (ranging from 70 to 96)';...
'SL0: Source level [dB/1µPa/1m] (ranging from 214 to 218)'};
dlg_title = 'Summary of ADCP characteristics - can be changed by users';
num_lines = 1;
def = {num2str(Kc);num2str(EC0);num2str(B);num2str(SL0)};

tmp = inputdlg(prompt,dlg_title,num_lines,def);
Kc=eval(char(tmp{1}));
EC0=eval(char(tmp{2}));
B=eval(char(tmp{3}));
SL0=eval(char(tmp{4}));

end
