%% WADV PACKAGE 7.1 -- function open_ctl_adcp_sontek
%%
% Used to load the header file corresponding to the data file loaded
function open_ctl_adcp_sontek

%%
% *Declaration of global variables* 

% _deployment parameters

global filein path_routines filein_path
global theta WN WS WB IMM ens_int f a_t ouv Kc NC0 B NE0 R1 Lt nbens binpos
global inst_freq mode_exp

%% Load the file
% 
cd(filein_path)
fid=fopen(strcat(filein,'ctl'));

%cd(path_routines);
line=1;

%% Read information
% 
while line<50
tline=fgetl(fid);

if line==4
%%
% *Read number of profiles*    
    nbens=str2double(tline(32:end));
end    
if line==17
%%
% *Read angle*    
    theta=str2double(tline(32:36));
end    

%%
% *Read number of cells*
if line==40
WN=str2double(tline(32:end));
end
%%
% *Read size of cells*
if line==41
WS=str2double(tline(32:36));
end
%%
% *Read blank size*
if line==42
WB=str2double(tline(32:36));
end
%%
% *Read height of the transducers from bed*
if line==43
IMM=str2double(tline(32:36));
end
%%
% *profile interval*
if line==46
ens_int=str2double(tline(32:35));
end

line=line+1;
end
fclose(fid)

cd(path_routines)

%%
% * additional parameters

%% SONTEK ADCP
% * deployment characteristics

R1=WB+WS;                  % position milieu premiere cellule courant [m]
Lt=R1-WB;                   % transmit pulse [m]
binpos=R1+WS*(0:WN-1)+IMM;


if mode_exp==0 % if mode_exp==1 values taken from wacp_config.txt
%%
% * instrument features

inst_freq_txt=char('500','750','1000','1500','3000');
[inst_freq_choice,presult]=listdlg('PromptString','ADCP Frequency (kHz)','ListString',inst_freq_txt);
inst_freq=inst_freq_txt(inst_freq_choice,:);
    


if strcmp(inst_freq,'3000')
%% !!! ATTENTION, DEFAULT VALUES BUT NOT CORRECT !!!
f=3000;           % frequence de l'ADCP [kHz]
a_t=0.051;         % diametre du transducteur [m] (exemple de valeurs par défaut)
%theta=25;       % orientation du faisceau [deg] déjà lu
ouv=99;         % ouverture angulaire du faisceau en Emission+Reception [deg] (exemple de valeurs par défaut donnée par Sontek)

Kc=0.42;    % relation dB/counts (valeur par defaut autour de 0.43)
NC0=45;     % niveau de bruit interne [counts] (valeur par defaut de 40 a 65 counts, determinable dans air)
B=70;       % niveau de bruit [dB/1microPa] (valeur par defaut de ?)
NE0=217;    %niveau emis [dB/1microPa] (valeur par defaut donnée par Sontek)

elseif strcmp(inst_freq,'1500')
    
f=1500;           % frequence de l'ADCP [kHz]
a_t=0.051;         % diametre du transducteur [m] (exemple de valeurs par défaut)
%theta=25;       % orientation du faisceau [deg] déjà lu
ouv=0.99;         % ouverture angulaire du faisceau en Emission+Reception [deg] (exemple de valeurs par défaut donnée par Sontek)

Kc=0.42;    % relation dB/counts (valeur par defaut autour de 0.43)
NC0=45;     % niveau de bruit interne [counts] (valeur par defaut de 40 a 65 counts, determinable dans air)
B=70;       % niveau de bruit [dB/1microPa] (valeur par defaut ?)
NE0=217;    %niveau emis [dB/1microPa] (valeur par defaut donnée par Sontek)
   
elseif strcmp(inst_freq,'1000')
    
f=1000;           % frequence de l'ADCP [kHz]
a_t=0.051;         % diametre du transducteur [m] (exemple de valeurs par défaut)
%theta=25;       % orientation du faisceau [deg] déjà lu
ouv=0.99;         % ouverture angulaire du faisceau en Emission+Reception [deg] (exemple de valeurs par défaut donnée par Sontek)

Kc=0.42;    % relation dB/counts (valeur par defaut autour de 0.43)
NC0=45;     % niveau de bruit interne [counts] (valeur par defaut de 40 a 65 counts, determinable dans air)
B=70;       % niveau de bruit [dB/1microPa] (valeur par defaut ?)
NE0=217;    %niveau emis [dB/1microPa] (valeur par defaut donnée par Sontek)

elseif strcmp(inst_freq,'750')
    
f=750;           % frequence de l'ADCP [kHz]
a_t=0.051;         % diametre du transducteur [m] (exemple de valeurs par défaut)
%theta=25;       % orientation du faisceau [deg] déjà lu
ouv=0.99;         % ouverture angulaire du faisceau en Emission+Reception [deg] (exemple de valeurs par défaut donnée par Sontek)

Kc=0.42;    % relation dB/counts (valeur par defaut autour de 0.43)
NC0=45;     % niveau de bruit interne [counts] (valeur par defaut de 40 a 65 counts, determinable dans air)
B=70;       % niveau de bruit [dB/1microPa] (valeur par defaut ?)
NE0=217;    %niveau emis [dB/1microPa] (valeur par defaut donnée par Sontek)

elseif strcmp(inst_freq,'500')
    
f=500;           % frequence de l'ADCP [kHz]
a_t=0.051;         % diametre du transducteur [m] (exemple de valeurs par défaut)
%theta=25;       % orientation du faisceau [deg] déjà lu
ouv=0.99;         % ouverture angulaire du faisceau en Emission+Reception [deg] (exemple de valeurs par défaut donnée par Sontek)

Kc=0.42;    % relation dB/counts (valeur par defaut autour de 0.43)
NC0=45;     % niveau de bruit interne [counts] (valeur par defaut de 40 a 65 counts, determinable dans air)
B=70;       % niveau de bruit [dB/1microPa] (valeur par defaut ?)
NE0=217;    %niveau emis [dB/1microPa] (valeur par defaut donnée par Sontek)

end


prompt = {'Kc : dB/counts conversion coefficient (ranging from 0.4 and 0.45)';'NC0 : Internal noise i.e. count value in air (ranging from 40 to 65)';'B : noise [dB/1microPa] (ranging from 70 to 96)'};
dlg_title = 'Summary of ADCP characteristics - can be changed by users';
num_lines = 1;
def = {num2str(Kc);num2str(NC0);num2str(B)};
tmp = inputdlg(prompt,dlg_title,num_lines,def);
Kc=eval(char(tmp{1}));
NC0=eval(char(tmp{2}));
B=eval(char(tmp{3}));
end


 