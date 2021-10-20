function open_adcp_sontek

global time ve vn vu dir spd a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean nbens
global path_routines filein_path filein 
global Batt Tb_adcp Depth_adcp Sb_adcp
global theta WN WS WB IMM ens_int f a_t ouv Kc NC0 B NE0 Lt binpos

open_ctl_adcp_sontek

cd(filein_path);
%% decription des variables
filetmp=load(strcat(filein,'hdr'));
%filetmp=filetmp;
time=datenum(filetmp(:,2),filetmp(:,3),filetmp(:,4),filetmp(:,5),filetmp(:,6),filetmp(:,7));

%%
% *Données sur l'orientation de l'adcp dans l'espace
pitch=filetmp(:,10);
roll=filetmp(:,11);

heading=filetmp(:,12);

%%
% Parametres auxiliaires
Tb_adcp=filetmp(:,13);
Sb_adcp=35.*ones(size(Tb_adcp));
Batt=filetmp(:,20);
Depth_adcp=filetmp(:,14); %en m
display('auxiliary parameters ok')

%% chargement des variables
% * vitesse
ve=load(strcat(filein,'ve')); % en m
vn=load(strcat(filein,'vn')); % en m
vu=load(strcat(filein,'vu')); % en m
ve=ve(:,2:end)/100;
vn=vn(:,2:end)/100;
vu=vu(:,2:end)/100;

display('velocity ok')

%%
% * signal noise ratio
sn1=load(strcat(filein,'sn1')); % en count
sn2=load(strcat(filein,'sn2')); % en count
sn3=load(strcat(filein,'sn3')); % en count
sn1=sn1(:,2:end);
sn2=sn2(:,2:end);
sn3=sn3(:,2:end);

display('SNR ok')

dir=load(strcat(filein,'dir')); % en degrés
spd=load(strcat(filein,'spd')); % en cm/s
dir=dir(:,2:end);
spd=spd(:,2:end)/100;

%%
% * Amplitude signal rétrodiffusé
a1=load(strcat(filein,'a1'));
a2=load(strcat(filein,'a2'));
a3=load(strcat(filein,'a3'));
a4=zeros(size(a1));
a4(:,:)=NaN;

a1=a1(:,2:end);
a2=a2(:,2:end);
a3=a3(:,2:end);
a4=a4(:,2:end);

amean=1/3*(a1+a2+a3);

display('amplitude ok')

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

display('correlation ok')

clear filetmp

cd(path_routines)
