function open_adcp

global time tmin tmax moy_ampl VE VN VV Dir Norme ampl_b1 ampl_b2 ampl_b3 ampl_b4 moy_cor nbin nbens hbin
global filein_name path_routines filein_path
global RDIBinSize RDIEnsInterval RDISecPerPing RDIPingsPerEns RDIBin1Mid
global AnBatt AnT100thDeg AnDepthmm


[filein_name,filein_path]=uigetfile({'*.mat';'*.*'},'Open data file....');
cd(filein_path);
load(filein_name);
cd(path_routines);

%% decription des variables

time=datenum(SerYear,SerMon,SerDay,SerHour,SerMin,SerSec);
hbin=RDIBin1Mid+RDIBinSize*(SerBins-1);
nbens=length(SerEnsembles);
nbin=length(SerBins);
bin=(1:nbin);
tmin=min(time);
tmax=max(time);
%%
% *Données sur l'orientation de l'adcp dans l'espace
pitch=AnP100thDeg/100;
roll=AnR100thDeg/100;
heading=AnH100thDeg/100;

%% chargement des variables
% * vitesse
VE=SerEmmpersec/1000; % en m
VN=SerNmmpersec/1000; % en m
VV=SerVmmpersec/1000; % en m

Dir=SerDir10thDeg/10; % en degrés
Norme=SerMagmmpersec/1000; % en m
%%
% * Pression
heau=AnDepthmm/1000; %en m
%%
% * Amplitude signal rétrodiffusé
ampl_b1=SerEA1cnt;
ampl_b2=SerEA2cnt;
ampl_b3=SerEA3cnt;
ampl_b4=SerEA4cnt;

moy_ampl=0.25*(ampl_b1+ampl_b2+ampl_b3+ampl_b4);
%%
% * Corrélation

cor_b1=SerC1cnt;
cor_b2=SerC2cnt;
cor_b3=SerC3cnt;
cor_b4=SerC4cnt;

moy_cor=0.25*(cor_b1+cor_b2+cor_b3+cor_b4);
%%
% * percent good

pg_b1=SerPG1;
pg_b2=SerPG2;
pg_b3=SerPG3;
pg_b4=SerPG4; % utilisé pour valider mesures si PG4>25%

%cor2_b1=cor_b1;
%cor2_b2=cor_b2;
%cor2_b3=cor_b3;
%cor2_b4=cor_b4;
%moy_cor2=moy_cor;

%% FILTRE DE SURFACE -- recherche de la dernière cellule avant la surface
%
% * plusieurs possibilités : 
%   (1) à partir du capteur de pression
%   (2) à partir du max de signal rétrodiffusé
%   (3) à partir d'un travail sur la corrélation des faisceaux


bv=char('Pression','Backscatter','Corrélation');
[vchoice,vresult]=listdlg('PromptString','Choisir le filtre de surface','ListString',bv);
nbv=length(vchoice);
if nbv>1
    diplay('attention choisir un seul filtre'); 
    bv=char('Pression','Backscatter','Corrélation');
    [vchoice,vresult]=listdlg('PromptString','Choisir le filtre de surface','ListString',bv);
    nbv=length(vchoice);
end


    pos=zeros(nbens,1);

if vchoice==1
    
    for i=1:nbens
        tmp=find(hbin<=heau(i));
        if isempty(tmp)
           pos(i)=NaN;
           continue
        else
        pos(i)=tmp(length(tmp)); % ou -1 pou s'assurer d'être réellement dans l'eau 
            if pos(i)<=0
            pos(i)=NaN;
            end
        end
    end

%% 
% * rv
    
elseif vchoice==2
        for i=1:nbens
           tmp=find(moy_ampl(i,:)==max(moy_ampl(i,:)));
           pos(i)=tmp(1);
           if pos(i)<=0
               pos(i)=NaN;
           end
        end
%%
% * Nicolas Chini        
%elseif vchoice==2
% il existe un echo même si pas d'eau au dessus de l'adcp
% dans ce cas l'écho est de l'ordre de 90dB
% on met à nan les échos inférieurs à 90dB
%moy_ampl(find(moy_ampl<90)) =nan; % 90 valeur bidon


%for i=1:nbens
%    dev_moy_ampl = sign(diff(moy_ampl(i,:))); % dev_moy_ampl vaut soit 1 soit 0 soit -1
%    n_dev = length(dev_moy_ampl);
    % on recherche le signe de dev_moy_ampl
    %i
%    change_sign = find(abs(diff(dev_moy_ampl))>=1);
%    figure
%    subplot(2,1,1)
%
%    plot(bin, moy_ampl(i,:))
%    title(['ens ' num2str(i)]);
%    subplot(2,1,2)
%    hold on
%    plot(bin(1:end-1), dev_moy_ampl)
%    plot(bin(change_sign), dev_moy_ampl(change_sign),'ro')
%    if isempty(change_sign)
        % le signal d'echo est monotone
 %       disp('  ')
 %       disp([ ' ens : ' num2str(i)])
 %       disp('signal d echo monotone  ')
 %       pos(i) = nan;
 %   elseif dev_moy_ampl(1:change_sign(1))<1 & length(change_sign)>1
 %           disp('  ')
 %           disp([ 'ens : ' num2str(i)])
 %           disp('youri  ')
 %           pos(i) = change_sign(2)+1;
 %   elseif dev_moy_ampl(1:change_sign(1))==1 & length(change_sign)>1
            % le signal est croissant
 %           disp('  ')
 %           disp([ 'ens : ' num2str(i)])
 %           disp('signal d echo croissant  ')
 %           pos(i) = change_sign(1)+1;
 %   else
 %           disp('  ')
 %           disp([ 'changement de signe a bin: ' num2str(change_sign(1))])
 %           disp('  ')
 %       pos(i) = change_sign(1)+1;

    
%%       

elseif vchoice==3
    pos=zeros(nbens,1);
tmp_cor=zeros(nbens,nbin-1);

for i=1:nbens
   
   tmp_cor(i,:)=moy_cor(i,2:nbin)-moy_cor(i,1:nbin-1);
   val=4; %% paramètre variable
   for j=2:nbin-2
       
       if (tmp_cor(i,j)<-val && tmp_cor(i,j-1)>val) || (tmp_cor(i,j-1)>val && tmp_cor(i,j+1)>val)
            pos(i)=j; % ou +1...à voir
            break
       end
       
       if j==nbin-2
            pos(i)=nbin-1;
       end
   end  
end
pos2=pos;
pos(pos==nbin-1)=NaN;
delta_pos=pos(2:nbens)-pos(1:nbens-1);

for i=1:nbens-1
    if abs(delta_pos(i))>=2
    %k=1;
    %while abs(pos(i+1+k)-pos(i))>2
    %k=k+1;
        if abs(pos(i+2)-pos(i))<=2
        pos(i+1)=floor(0.5*(pos(i)+pos(i+2)));
        elseif abs(pos(i+3)-pos(i))<=2
        pos(i+1)=floor(0.5*(pos(i)+pos(i+3)));
        else
        pos(i+1)=floor(0.5*(pos(i)+pos(i+4)));
        end
    %pos(i+1)=floor(0.5*(pos(i)+pos(i+k+1)));
    delta_pos=pos(2:nbens)-pos(1:nbens-1);
    
    end
end


nanpos=isnan(pos);

for i=2:nbens-1
   if isnan(pos(i))
       pos(i)=(pos(i+1)+pos(i-1))/2;
   end

   %if isnan(pos(i)) & pos(i+1)>1
   %pos(i+1)=NaN;
   %end
   
   %if isnan(pos(i)) & pos(i-1)>1
   %pos(i-1)=NaN;
   %end
       
end

end

%% application du filtre choisi
%

for i=1:nbens
moy_ampl(i,bin>pos(i))=NaN;
VE(i,bin>pos(i))=NaN;
VN(i,bin>pos(i))=NaN;
VV(i,bin>pos(i))=NaN;
Dir(i,bin>pos(i))=NaN;
Norme(i,bin>pos(i))=NaN;
ampl_b1(i,bin>pos(i))=NaN;
ampl_b2(i,bin>pos(i))=NaN;
ampl_b3(i,bin>pos(i))=NaN;
ampl_b4(i,bin>pos(i))=NaN;
moy_cor(i,bin>pos(i))=NaN;
if isnan(pos(i))
moy_ampl(i,:)=NaN;
VE(i,:)=NaN;
VN(i,:)=NaN;
VV(i,:)=NaN;
Dir(i,:)=NaN;
Norme(i,:)=NaN;
moy_cor(i,:)=NaN;
ampl_b1(i,:)=NaN;
ampl_b2(i,:)=NaN;
ampl_b3(i,:)=NaN;
ampl_b4(i,:)=NaN;
end
end