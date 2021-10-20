%% WADCP PACKAGE 0.3 -- function backscatter
%%
%
% Backscatter signal processing toolbox
%
% _Written by C. Tessier, IFREMER |contact : caroline.tessier@ifremer.fr 

function backscatter(bsc)

%warning off all

%%
% _deployment parameters_
global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 R1 Lt 
global time xrange yrange nt nk dz binpos binpos_bt inst_mode inst_type mode_exp varia_batt qst_nf backs_calib_choice 
global Batt Tb_adcp Depth_adcp Sb_adcp Mgl R thetaR  ve spd vn
global RL_grid SL_grid TL_geo Cst_geo TL_w BI fidlog
keyboard
%global adjust_sdt
% 
% datevec(time(734:736))
% datevec(time(1492:1497))
% datevec(time(1642:1647))
% datevec(time(2595:2599))
% datevec(time(2641:2645))
% 
% figure
% plot((1:length(time)-1),diff(time))

%%
% *auxiliary variables*

thetaR=theta*pi/180;            % beam angle in radians 0.3491 rad       
lambda=1500/f/1000;             % mean wave length 
%omega=2*pi*f*1000;              % angular velocity
%k=2*pi/lambda;                  % wave number
%r0=k*(a_t/2)^2/2;               % limite theorique champs proche =1.67 m
r0=a_t^2/lambda/2;              % en pratique: limite champs proche =1.06m
r0=1;
%PSI=8*pi*(lambda/pi/a_t)^2;     % ouverture solide equivalente theorique
phi=ouv*pi/180;                 %ouverture equivalente mesuree
PSI=pi*(phi/2)^2;               %ouverture solide equivalente mesuree


%%
% *Cells location from ADCP*
%(intensite dans derniere moitie cellule si cellules petites, 
% ou dernier quart si cellules grandes)

for ik=1:WN
   R(ik)=WB+Lt+(ik-1)*WS+1/4*WS;      % backscatter (centre derniere moitie cellule)
   Ruv(ik)=WB+Lt+(ik-1)*WS;           % current velocity (centre cellule)
end
Rv=R/cos(thetaR);                     % wave travel distance
Vol=PSI*Rv.^2*(1/2*WS);               % cell volume for backscatter signal
Cst_geo=10*log10(PSI*1/2*WS);            % Constant Cst_geo=-42.32 a f=1250 kHz et 2theta=0.99deg
                                
%% Backscatter processing
% *Echo received signal [Counts]
EC=bsc(:,:);
bsc=[];
clear bsc

% *Received Level [dB/1microPa]
RL=B+(EC-EC0).*Kc;

if mode_exp==0
%%
% *Correction of the influence of batteries power* 

qst_proc=questdlg('Do you have variability of batteries power ? (if Yes, please CHECK the relation of Source Level SL line 60 in backscatter.m)', ...
 'Control of Source Level','Yes','No','No');
 switch qst_proc
     case 'Yes'
         varia_batt=1;
     case 'No'
         varia_batt=0;
 end

    
end
 
if varia_batt==1 
    
    if strcmp(inst_type,'Nortek')
        display('Battery depletion Nortek formula')
      SL=SL0+20*log10(Batt/13.2); %specifique PrevmierD4 Nortek 
    else
    display('cf C Tessier lab observations - PhD Thesis')
% Profiles where batteries power was acceptable
I=find(Batt>=100); %Batt en count ! RDI files
d_bat=time(I);
tensionADC=Batt(I);
%calib in tank water ADCP RDI 1200 kHZ n�4285 : 
%SL=-2.59*10^-4*tensionADC.^2+0.1278.*tensionADC+203.62;
%calib in tank water ADCP RDI 1200 kHZ mooring BZHS 2005 :
%SL=-1.149*10^-4*tensionADC.^2+0.0707.*tensionADC+207.47;  %beam1
SL=-1.212*10^-4*tensionADC.^2+0.0729.*tensionADC+206.78;  %beam3

    end
elseif varia_batt==0
% prompt = {'SL0: Source level [dB/1�Pa/1m] (ranging from 214 to 218)'};
% dlg_title = 'ADCP characteristics - can be changed by users';
% num_lines = 1;
% def = {num2str(SL0)};
% tmp = inputdlg(prompt,dlg_title,num_lines,def);
% SL0=eval(char(tmp{1})); 
display('Battery depletion not accounted for')
SL=SL0*ones(length(time),1);
end



tensionADC=[];
clear tensionADC


%%
% *Sound speed calculation*
Cf=celeriteChen(Depth_adcp,Tb_adcp,Sb_adcp); 


%%
% *Data interpolation on a user-defined regular grid 
% (time of ADCP, vertical position of backscatter data adjusted)

% -vertical axis:
 y=IMM+R; 
 dz=1/2*WS;          % vertical step [m] because value R is in the center of the last middle of the cell and is available for this last middle of the cell
 yrange=(y(1):dz:max(y)+dz); % to be just on the vertical level refered at the bottom 

 %modifs RV temporaires
 
 yrange=y;
 %dz=WS;
 %yrange=y;
% -time axis:
 xrange=time;                        %raw time of ADCP data en TU

 % dt=mean(time(2:end)-time(1:end-1)); %not used, just to check that dt is constant and time is monotone
%xrange=time(1):dt:time(end);
% dt=0.003472;   % 5mn;
% dt=0.006944;   %10mn;
%-------------------------------------
 nt=length(xrange); 
 nk=length(yrange);
 [X,Y]=meshgrid(xrange,yrange); %used for the 2D plot

 %vertical interpolation to place backscatter data every middle cell :
 %RL_grid=interp2(time,y,RL',X,Y);            %Interpolated Received Level 
 
 RL_grid=RL';
 
 %time interpolation of energy
% if varia_batt==1 
    SL_grid(1:nt)=interp1(time,SL,xrange);        %Interpolated Emitted Level
%  elseif varia_batt==0
%     SL_grid(1:nt)=SL0;  %if constant Emitted level (to be ask to user)
%  end

X=[];
Y=[];
NR=[];
clear X Y NR

%%
% * Transmission loss (along beam path) * 
% TL_mean = TL_geo +TL_w = Total loss -20log(R)
% because simplification of 10log10(Vol) to compute Cst_geo
%
% water absorption alpha_w [Francois & Garrison, 1985]

% RV : attention pas depth mais profondeur de la cellule
% for i=1:length(yrange)
alpha_wdb=equ_att_son_garrison(f,Depth_adcp,Tb_adcp,Sb_adcp);

if mode_exp==0
% near field correction [Downing et al. 1995]
%!!!! if ADCP is a Zed-Hed, so no ringing, blank very small and nothing to correct !!!
qst_nf=questdlg('Do you want to use Near Field correction for the estimation of signal attenuation ? (NO needed if the ADCP is a Zed-Hed, ask to constructor if yours is one or no)', ...
    'Correction of Near Field propagation','Yes','No','Yes');
end

switch qst_nf
     case 'Yes'
        Z_chp=((yrange-IMM)/cos(thetaR))/r0;
        PSI_chp=(1+1.35*Z_chp+(2.5*Z_chp).^3.2)./(1.35*Z_chp+(2.5*Z_chp).^3.2);
     case 'No'
        PSI_chp(1:nk)=1.0;       
 end
 
%geometrical attenuation for spherical spreading (-20log(R)):
TL_geo=20*log10((yrange-IMM)/cos(thetaR).*PSI_chp); 
TL_geo=TL_geo(1:nk)'*ones(1,nt); % to have a 2D matrix now

% water attenuation profiles 
    %if assumed constant profile : 

    TL_w=2*alpha_wdb(1:nt)*(yrange-IMM)/cos(thetaR);
    
    TL_w=TL_w';

    %if mean constant value : 
    %TL_w=2*mean(alpha_wdb)*(yrange-IMM)/cos(thetaR);     
    %TL_w=TL_w(1:nk)'*ones(1,nt); % to have a 2D matrix now

%% si attenation par eau    
% *Backscatter index from ADCP measurements*
BI=zeros(nk,nt);

for i=1:nk,  
    BI(i,1:nt)=RL_grid(i,1:nt)-SL_grid(1:nt)+TL_geo(i,1:nt)+TL_w(i,1:nt)-Cst_geo;
end

    
    
    %% attenuation par le s�diment et eau
%     backscatter_iter
    


%%




display('Scattering Index processed')
display('Choice of calibration methodology now')
if mode_exp==0
backs_calib_txt=char('Empirical [TBD series]','Empricial [TBD profiles]',...
    'Therorical Iterative method ','Mixte [TBD series]','Mixte [TBD profiles]','Already calibrated','None');
[backs_calib_choice,presult]=listdlg('PromptString','choice of methodology for backscatter signal calibration','ListString',backs_calib_txt);
end



    fprintf(fidlog,'Signal parameters (Kc, EC0, B, SL0) : %2.2f, %2.2f, %2.2f, %2.2f\n', [Kc EC0 B SL0]);
    fprintf(fidlog,'Battery effect : %1.0f\n', varia_batt);    
    fprintf(fidlog,'Near field correction : %s\n', qst_nf);       
    fprintf(fidlog,'Turbidity Calibration option : %1.0f\n', backs_calib_choice);    


switch backs_calib_choice
    case 1
% *empirical method to evaluate sediment concentration from TBD time series*
     adjust_sdt=0;
     backscatter_emp(1);
    case 2 
% *empirical method to evaluate sediment concentration from TBD vertical profiles*
     adjust_sdt=0;
    backscatter_emp(2);
    case 3 
% *Iterative method to evaluate sediment concentration AND attenuation*
    backscatter_iter;
    case 4 
% *Empirical (with TBD series) + Iterative methods to evaluate sediment concentration AND attenuation*
    adjust_sdt=1;
    backscatter_iter_empir_tbd_bottom;
    case 5
% *Empirical (with TBD profiles) + Iterative methods to evaluate sediment concentration AND attenuation*
    %backscatter_iter_empir_tbd_profile;
    case 6
     adjust_sdt=0;
     backscatter_emp(3);
    case 7
        Mgl=zeros(size(BI));
end
% 
% figure
% imagesc(xrange,yrange,BI)
% ylim([0 max(Depth_adcp)])
% set(gca,'YDir','normal')
% title('BI backscatter.m')
% datetick('x',19)
% caxis([-70 -20])
% colorbar
% shading flat

figure
subplot(211)
imagesc(xrange,binpos,ve')
ylim([0 30])
set(gca,'YDir','normal')
title('Current speed')
datetickzoom
caxis([-1 1])
colorbar
shading flat

subplot(212)
imagesc(xrange,binpos,vn')
ylim([0 30])
set(gca,'YDir','normal')
title('Current speed')
datetickzoom
caxis([-1 1])
colorbar
shading flat
% save('BI.mat','BI','time','yrange')







