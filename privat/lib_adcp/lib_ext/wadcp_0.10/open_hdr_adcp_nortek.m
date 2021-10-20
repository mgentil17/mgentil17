%% WADV PACKAGE 7.1 -- function open_hdr_adcp_nortek
%%
% Used to load the header file corresponding to the data file loaded
function open_hdr_adcp_nortek

%%
% *Declaration of global variables* 

global filein path_routines filein_path
global theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 R1 Lt binpos
global inst_freq mode_exp

%% Load the file
% 
cd(filein_path)
fid=fopen(strcat(filein,'hdr'));

%cd(path_routines);
line=1;

%% Read information
% 
while line<50
tline=fgetl(fid);


%%
% *Read number of cells*
if line==11
    tline(32:end)
WN=str2double(tline(32:end));
end
%%
% *Read size of cells*
if line==12
    tline(32:end-2)
WS=str2double(tline(32:end-2))/100;
end
%%
% *Read blank size*
if line==16
    tline(32:end-1)
WB=str2double(tline(32:end-1));
end

%%
% *profile interval*
if line==10
    tline(38:end-3)
ens_int=str2double(tline(38:end-3));
end

line=line+1;
end
fclose(fid)


if mode_exp==0 % if mode_exp==1 values read from wadcp_config.txt
%%
% * additional parameters

%% NORTEK ADCP
% * deployment characteristics

prompt = {'Enter the distance from the bed - IMM (m)'};
dlg_title = 'Deployment characteristics';
num_lines = 1;
def = {'0.2'};
tmp = inputdlg(prompt,dlg_title,num_lines,def);
IMM=eval(char(tmp{1}));



%%
% * instrument features

inst_freq_txt=char('600','1000');
[inst_freq_choice,presult]=listdlg('PromptString','ADCP Frequency (kHz)','ListString',inst_freq_txt);
inst_freq=inst_freq_txt(inst_freq_choice,:);


if strcmp(inst_freq,'1000')

f=1000;           % frequence de l'ADCP [kHz]
a_t=0;         % diametre du transducteur [m] (exemple de valeurs par défaut)
theta=25;       % orientation du faisceau [deg]
ouv=1.7;         % ouverture angulaire du faisceau en Emission+Reception [deg] (exemple de valeurs par défaut : RDI 2.87deg pr 300kHz; 0.99deg pr 1200kHz)

Kc=0.42;    % relation dB/counts (valeur par defaut autour de 0.43)
EC0=45;     % niveau de bruit interne [counts] (valeur par defaut de 40 a 65 counts, determinable dans air)
B=70;       % niveau de bruit [dB/1microPa] (valeur par defaut de ?)
SL0=196;    %niveau emis [dB/1microPa] (valeur par defaut donnée par Nortek)

%elseif strcmp(inst_freq,'600')
    elseif eval(inst_freq)==600
f=600;           % frequence de l'ADCP [kHz]
a_t=0;         % diametre du transducteur [m] (exemple de valeurs par défaut)
theta=25;       % orientation du faisceau [deg]
ouv=3.4;         % ouverture angulaire du faisceau en Emission+Reception [deg] (exemple de valeurs par défaut : RDI 2.87deg pr 300kHz; 0.99deg pr 1200kHz)

Kc=0.42;    % relation dB/counts (valeur par defaut autour de 0.43)
EC0=45;     % niveau de bruit interne [counts] (valeur par defaut de 40 a 65 counts, determinable dans air)
B=70;       % niveau de bruit [dB/1microPa] (valeur par defaut ?)
SL0=195;    %niveau emis [dB/1microPa] (valeur par defaut donnée par Nortek)
   
end


 

prompt = {'Kc : dB/counts conversion coefficient (ranging from 0.4 and 0.45)';'EC0 : Internal noise i.e. count value in air (ranging from 40 to 65)';'B : noise [dB/1microPa] (ranging from 70 to 96)'};
dlg_title = 'Summary of ADCP characteristics - can be changed by users';
num_lines = 1;

def = {num2str(Kc);num2str(EC0);num2str(B)};
tmp = inputdlg(prompt,dlg_title,num_lines,def);
Kc=eval(char(tmp{1}));
EC0=eval(char(tmp{2}));
B=eval(char(tmp{3}));
end

Lt=WS;                   % transmit pulse [m] from AWAC Documentation
binpos=WB+WS*(1:WN)+IMM;

cd(path_routines)
 