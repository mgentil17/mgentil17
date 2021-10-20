%% WADCP PACKAGE 0.3 -- function backscatter_iter_empir_tbd_bottom
%%
%
% Backscatter signal processing toolbox
%
% _Written by C. Tessier, IFREMER |contact : caroline.tessier@ifremer.fr 

function backscatter_iter_empir_tbd_bottom
% * iterative method and calibration con TBD * 

global nt nk RL_grid SL_grid TL_geo Cst_geo TL_w BI Mgl
global IMM yrange dz xrange

%%
% * sediment absorption constantes * 
%   alpha_s=2*(zeta_vdB+zeta_ddB)*M 
%% PICLOD corrected
global f nf
diameter=55; % diametre du sediment en suspension en micron
nf=2; % choix de la dimension fractale : entre 1 et 3, en milieu naturel varie entre 1.6 et 2.5
    [zeta_vdB,zeta_ddB,Cst_sdt]=sdt_absorption(f,diameter);

%%
%%[zeta_vdB,zeta_ddB,Cst_sdt]=sdt_absorption;
%% End of correction

%%

%FIRST STEP : adjust empirical calibration for the first cell time series
%             with TBD data (after computing of sediment attenuation )

% - estimate ivcal from TBD/BI, choose of calibration period
  [ivcal,pval,turb_alt,t_calib0,turb_int0,nens_calib,t_calib,turb_int]=cal_turbidity ;

   M0=turb_int';       %g/l        %tbd data time serie 
   BI0=BI(turb_alt,t_calib);       %vector BI (first cell, tbd calib time period)
   M0_adcp=10.^((ivcal(1)*BI0+ivcal(2))./10); %g/l
   %a=ivcal(1);b=ivcal(2);
   %M0=(10^((a*BI0+b)/10))*10^-3;
   d0=min(yrange-IMM); %first cell-transducer distance 
% - calculate qs and new BI1 
    qs0=2*(zeta_vdB+zeta_ddB)*M0*d0; %sediment attenuation
    BI1=BI0+qs0 ; 
    %calculate NEW coefficient calibration 10log(Mtbd)=a*BI+b
    ivcal=polyfit(BI1,10*log10(M0),1);
    %Obtained ADCP values of concentration
    pval=polyval(ivcal,BI1); %10log10(g/l)
    M1=10.^(pval./10); %g/l
 %%%   figure; plot(t_calib,turb_int,'r',t_calib,M0_adcp','g',t_calib,M1,'b');legend('TBD','M0adcp','M1adcp',0);
    %correlation matrix
    coef=corrcoef(pval,10*log10(turb_int'));
    % relative error 
    er=abs(10.^(pval/10)-turb_int');
    er1=sum(er./turb_int')/nens_calib;
    er2=sqrt(sum(er.^2)/sum(turb_int'.^2));  %=17.8%  
    % figure plot 
    %cal_turbidity_plot
    
%- on garde les valeurs calculees au niveau du TBD 
Mgl(turb_alt,t_calib)=M1;
BI1(turb_alt,t_calib)=BI1;
%- on garde les paramètres de calibration BI/TBD tbd en g/l
a=ivcal(1);b=ivcal(2);
    
% %OLD TEST of iterations of TBD empirical calibration----------------    
% ite=0;
% epsM=100;
% while epsM>=0.001
%       ite=ite+1;
% % - calculate qs and new BI1 
%     qs0=2*(zeta_vdB+zeta_ddB)*M0*d0; %sediment attenuation
%     BI1=BI0+qs0' ; 
% % - new estimation of ivcal from TBD/BI1 
% %   [ivcal,pval]=cal_turbidity_ite(BI1,t_calib0,turb_int0) ;
%     %calculate coefficient calibration 10log(Mtbd)=a*BI+b
%     %ivcal=polyfit(BI(turb_alt,t_calib0),10*log10(turb_int0)',1);
%     ivcal=polyfit(BI1,10*log10(M0)',1);
%     %ADCP values (10log(Madcp))
%     pval=polyval(ivcal,BI1)'; %10log10(g/l)
%     epsM=max(abs(M0-10.^(pval./10))); 
%     M0=10.^(pval./10); 
%     BI0=BI1;
% end
%     cst_calib=BI1-BI(turb_alt,t_calib0);
%     figure;plot(BI0,cst_calib)
% %-------------------------------------------------------------
    

%%
%SECOND STEP : iterative calculation of profiles BI/qs and then use 10log10(M)=a BI +b


hh = waitbar(0,'Iteration of profiles in progress...');

clear qs MM M M1 %(M et M1 en [g/l]; MM=10log10(g/l))
%M(1:nk,1:nt)=NaN;
%MM(1:nk,1:nt)=NaN;

for t=1:nt,
waitbar(t/nt,hh,strcat('Iteration of profiles in progress...  ',int2str(t/nt*100),'%'));

%premiere cellule ADCP : on repart de BI de base avec parametres calib TBD precedent
M0=10^((a*BI(turb_alt,t)+b)/10);  %g/l  %concentration initiale comme si qs pris en compte mais avec valeur calib finale
qs0=2*(zeta_vdB+zeta_ddB)*M0*d0;  %recalcul qs sur toute la serie ADCP avec calibration approximative 
BI0=BI(turb_alt,t)+qs0;            %recalcul BI pour corriger de qs 

M1=M0; %g/l
%qs(turb_alt,t)=qs0;
%BI1(turb_alt,t)=BI0;
%for i=turb_alt+1:nk, %iteration profils pour estimer qs/BI 
for i=turb_alt:nk, %iteration profils pour estimer qs/BI 
   epsM=100;
   ite=0;
   while epsM>=0.000001 %= 0.001mg/l
      ite=ite+1;
      if i==turb_alt 
         qs(turb_alt,t)=qs0;
         BI1(turb_alt,t)=BI0;
      else
      qs(i,t)=qs(i-1,t)+2*(zeta_vdB+zeta_ddB)*M1*dz;
      %BI0=RL_grid(i,t)-SL_grid(t)+TL_geo(i,t)+TL_w(i,t)-Cst_geo;%-Cst_sdt;
      BI1(i,t)=BI(i,t)+qs(i,t);
      end 
      MM(i,t)=(a*(BI1(i,t))+b); %10log10[g/l]
      M(i,t)=10^(MM(i,t)/10);   %g/l
      epsM=abs((M(i,t)-M1));    %g/l
      M1=M(i,t);
   end
   itef(i,t)=ite; %final number of iteration 
end
end
Mgl=M; %g/l
BI_init=BI; % initial value of BI before to correct from sediment attenuation
BI=BI1;     %!!!!!!!!reactualisation backscatter index BI !!!!!!!!

%FIGURES TO CONTROL RESULTS
%----------------------------------------------------
%figure; imagesc(Mm1*1000);axis xy; caxis([0 100]);colorbar
figure;imagesc(itef); axis xy; caxis([0 4]);colorbar;title('number of iterations');
figure;imagesc(BI_init-BI1); axis xy; colorbar;

figure; plot(xrange(t_calib),turb_int,'r',xrange(t_calib),M0_adcp,'b',xrange,Mgl(turb_alt,:),'g');
legend('TBD','M0adcp calib tbd','Madcp ite',0);Ylim([0 max(turb_int)]);
set(gca,'Xlim',[fix(xrange(t_calib(1))) fix(xrange(t_calib(end))+1)],'Xtick',[fix(xrange(t_calib(1))):1:fix(xrange(t_calib(end))+1)]);
datetick('x',19,'keeplimits','keepticks');
xlabel(num2str(max(datevec(xrange(1)))));

%-----------------------------------------------------



close(hh)
display('End of Mixte 1 Method')

