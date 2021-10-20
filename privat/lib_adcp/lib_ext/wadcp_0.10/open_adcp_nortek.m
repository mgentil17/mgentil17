function open_adcp_nortek

global time ve vn vu a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean na1 na2 na3 na4 spd dir
global filein path_routines filein_path
global Batt Tb_adcp Depth_adcp Sb_adcp
global nbens ens_int

cd(filein_path);

%% lecture du fichier .csv, sans l'en tête et avec 6 colonnes pour le temps
%% (JJ MM AA HH MM SS)

filetmp=load(strcat(filein(1:end-1),'-P-PCSV.mat.csv'));

time=datenum(filetmp(:,3),filetmp(:,2),filetmp(:,1),filetmp(:,4),filetmp(:,5),filetmp(:,6));
heading=filetmp(:,8);
pitch=filetmp(:,9);
roll=filetmp(:,10);
Depth_adcp=filetmp(:,11);
Tb_adcp=filetmp(:,12);
Sb_adcp=35.*ones(size(Tb_adcp));
Batt=filetmp(:,7);

figure
plot(time)

%lecture des capteurs analogiques : turbidimetre wetlabs et altimetre
%tritech

analog1=filetmp(:,13);
analog4=filetmp(:,14);
nbens=length(time);
% filetmp=load(strcat(filein,'whd'));
% 
% 
% %% decription des variables
% 
% time=datenum(filetmp(:,3),filetmp(:,1),filetmp(:,2),filetmp(:,4),filetmp(:,5),filetmp(:,6));
% nbens=length(time);
% 
% %%
% % *Données sur l'orientation de l'adcp dans l'espace
% pitch=filetmp(:,13);
% roll=filetmp(:,14);
% heading=filetmp(:,12);
% 
% %%
% % Parametres auxiliaires
% Tb_adcp=filetmp(:,17);
% Sb_adcp=35.*ones(size(Tb_adcp));
% Batt=filetmp(:,10);
% Depth_adcp=filetmp(:,15); %en m
% 
% %% chargement des variables
% % * noise amplitude
% 
% na1=filetmp(:,19);
% na2=filetmp(:,20);
% na3=filetmp(:,21);
% na4=filetmp(:,22);
% 
% 
% clear filetmp



%% chargement des variables
% * vitesse
ve=load(strcat(filein,'v1')); % en m
vn=load(strcat(filein,'v2')); % en m
vu=load(strcat(filein,'v3')); % en m

size(time)
size(ve)

% on met toutes les variables de la dimension des profils courants

% 
% nbens=length(ve(:,1));
% time2=time(1)+(0:nbens-1)*ens_int/(3600*24)+1/86400;
% time2=time2';
% 
% toto=[1; diff(time)];
% Tb_adcp(toto==0)=[];
% Sb_adcp(toto==0)=[];
% Batt(toto==0)=[];
% Depth_adcp(toto==0)=[];
% heading(toto==0)=[];
% pitch(toto==0)=[];
% roll(toto==0)=[];
% time(toto==0)=[];


[dir,spd]=cart2pol(ve,vn);
dir=dir*180/pi;


% Tb_adcp=interp1(time,Tb_adcp,time2);
% Sb_adcp=interp1(time,Sb_adcp,time2);
% Batt=interp1(time,Batt,time2);
% Depth_adcp=interp1(time,Depth_adcp,time2);
% heading=interp1(time,heading,time2);
% pitch=interp1(time,pitch,time2);
% roll=interp1(time,roll,time2);
% time=time2;
% clear time2

%%
% * Amplitude signal rétrodiffusé
a1=load(strcat(filein,'a1'));
a2=load(strcat(filein,'a2'));
a3=load(strcat(filein,'a3'));
a4=zeros(size(a1));
a4(:,:)=NaN;



amean=1/3*(a1+a2+a3);

%%
% * Corrélation

c1=zeros(size(a1));
c1(:,:)=NaN;
c2=zeros(size(a1));
c2(:,:)=NaN;
c3=zeros(size(a1));
c3(:,:)=NaN;
c4=zeros(size(a1));
c4(:,:)=NaN;


cmean=1/3*(c1+c2+c3);

cd(path_routines);


