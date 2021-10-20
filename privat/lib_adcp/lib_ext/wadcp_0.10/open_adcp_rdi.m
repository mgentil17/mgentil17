%% WADCP PACKAGE 0.3 -- function open_adcp_rdi
%%
%
% Backscatter signal processing toolbox
%
%

function open_adcp_rdi

global time ve vn vu dir spd a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean
global filein path_routines filein_path
global Batt Tb_adcp Depth_adcp Sb_adcp
global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 R1 Lt binpos binpos_bt inst_mode
global inst_freq mode_exp


cd(filein_path);

load(strcat(filein,'.mat'));
cd(path_routines);

%% decription des variables

time=datenum(SerYear+2000,SerMon,SerDay,SerHour,SerMin,SerSec);
% 
% 
% toto=ones(length(time),1)/86400/2;
% toto(1:2:end)=0;


for i=2:length(time)
    if time(i)==time(i-1)
       time(i)=time(i)+1/86400/2;
    end
end


WS=RDIBinSize;
nbens=length(SerEnsembles);
WN=length(SerBins);


%%
% *Données sur l'orientation de l'adcp dans l'espace
pitch=AnP100thDeg/100;
roll=AnR100thDeg/100;
heading=AnH100thDeg/100;

%%
% Parametres auxiliaires
Tb_adcp=AnT100thDeg/100;
Sb_adcp=35.*ones(size(Tb_adcp));
Batt=AnBatt;
Depth_adcp=AnDepthmm/1000; %en m

%% chargement des variables
% * vitesse
ve=SerEmmpersec/1000; % en m/s
vn=SerNmmpersec/1000; % en m/s
vu=SerVmmpersec/1000; % en m/s

dir=SerDir10thDeg/10; % en degrés
spd=SerMagmmpersec/1000; % en m


%%
% * Amplitude signal rétrodiffusé
a1=SerEA1cnt;
a2=SerEA2cnt;
a3=SerEA3cnt;
a4=SerEA4cnt;

amean=0.25*(a1+a2+a3+a4);
%%
% * Corrélation

c1=SerC1cnt;
c2=SerC2cnt;
c3=SerC3cnt;
c4=SerC4cnt;

cmean=0.25*(c1+c2+c3+c4);

clear SerYear SerMon SerDay SerHour SerMin SerSec RDIBinSize SerBins SerEnsembles SerBins AnP100thDeg
clear AnR100thDeg AnH100thDeg AnT100thDeg AnBatt AnDepthmm SerEmmpersec SerNmmpersec SerVmmpersec SerDir10thDeg
clear SerMagmmpersec SerEA1cnt SerEA2cnt SerEA3cnt SerEA4cnt SerC1cnt SerC2cnt SerC3cnt SerC4cnt

%%
% * additional parameters

if mode_exp==0  % fi mode_exp==1 read from wadcp_config.txt

%% RDI ADCP
% * deployment characteristics

prompt = {'Enter the distance from the bed - IMM (m)';'Enter the blanking distance - WB (m)'};
dlg_title = 'Deployment characteristics';
num_lines = 1;
def = {'0.5';'0.4'};
tmp = inputdlg(prompt,dlg_title,num_lines,def);
IMM=eval(char(tmp{1}));
WB=eval(char(tmp{2}));



%%
% * instrument features

inst_freq_txt=char('1200','600','300');
[inst_freq_choice,presult]=listdlg('PromptString','ADCP Frequency (kHz)','ListString',inst_freq_txt);
if inst_freq_choice==1
inst_freq='1200';
elseif inst_freq_choice==2
inst_freq='600';
elseif inst_freq_choice==3
inst_freq='300';
end

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


R1=RDIBin1Mid;                  % position milieu premiere cellule courant [m]
Lt=R1-WB;                   % transmit pulse [m]
ens_int=RDIEnsInterval;         % ensemble interval [s]

binpos=R1+WS*(0:WN-1)+IMM;
if inst_mode==1
binpos_bt=-binpos+binpos(end);
end
clear RDIBin1Mid RDIEnsInterval
