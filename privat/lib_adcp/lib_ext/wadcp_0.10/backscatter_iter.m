
%% WADCP PACKAGE 0.1 -- function backscatter_iter
%%
%
% Backscatter signal processing toolbox
%
% _Written by C. Tessier, IFREMER |contact : caroline.tessier@ifremer.fr 

function backscatter_iter

global nt nk RL_grid SL_grid TL_geo Cst_geo TL_w Mgl BI f
global IMM yrange dz xrange R WS thetaR nf time_Mgl

%%
% * sediment absorption constants * 
%   alpha_s=2*(zeta_vdB+zeta_ddB)*M 

% [zeta_vdB,zeta_ddB,Cst_sdt]=sdt_absorption;

deltaR=ones(size(R))*WS;
deltaR(1)=R(1);
%% parametre du sediment
diameter=55; % diametre du sediment en suspension en micron
nf=2; % choix de la dimension fractale : entre 1 et 3, en milieu naturel varie entre 1.6 et 2.5

%%

    [zeta_vdB,zeta_ddB,Cst_sdt]=sdt_absorption(f,diameter);
  %  [zeta_vdB(id),zeta_ddB(id),Cst_sdt(id)]=sdt_absorption_primary(f,diameter(id));



%%

% * iterative method * 
hh = waitbar(0,'Iterative Method in progress...');

%clear qs MM M M1 %(M et M1 en [kg/m3==g/l]; MM=10log10(M))
d0=min(yrange-IMM); %first cell-transducer distance (ex 2005 : transd-Y(7)=0.1+4*0.25=1.10m)
% M0=10*10^-3; % initial concentration in this area

    
%%
subpoints=1;


Mgl=ones(nk,floor(nt/subpoints))*NaN;
BI=ones(nk,floor(nt/subpoints))*NaN;
time_Mgl=zeros(floor(nt/subpoints),1);
%%


%qs0=2*(zeta_vdB+zeta_ddB)*M0*d0;
%qs(1,1:nt)=qs0;
% M(1:nk,1:nt)=NaN;
% MM(1:nk,1:nt)=NaN;
for k=1:floor(nt/subpoints),
    t=k*subpoints;
    time_Mgl(k)=xrange(t);
%MM1=10*log10(1*10^-3); %1 mg/l
waitbar(t/nt,hh,strcat('Iterative Method in progress...--',int2str(t/nt*100),'%'));
%M1=1*10^-3; %1 mg/l
    for i=1:nk
    
        spmloop=0.001;
        deltaspm=100000;
%         iter=0;
    
%% PICLOD correction!
% Need to cum sum the Ctmps along the path
Ctmp = 0;

        while deltaspm>0.001
        
            Ctmp=Ctmp + 2*(zeta_vdB+zeta_ddB)*spmloop*deltaR(i)/cos(thetaR)';  %% PICLOD corrected
%            Ctmp=2*(zeta_vdB+zeta_ddB)*spmloop*deltaR(i)/cos(thetaR)';
            IVtmp=RL_grid(i,t)-SL_grid(t)+TL_geo(i,t)-Cst_geo+TL_w(i,t)+Ctmp;
            % CMEStmp(id)=10.^((aiv*IVtmp(id)+biv)/10);
            CMEStmp=10.^(IVtmp/10)/Cst_sdt;
            deltaspm=abs(CMEStmp-spmloop);
            spmloop=CMEStmp;
            % iter=iter+1;
        end

	ctmp_check(i, k)=Ctmp;
        Mgl(i,k)=CMEStmp;
        BI(i,k)=IVtmp;

    end
end
%Mgl=M;
%Mgl=(10.^(MM./10));

% %%
% 
% for t=1:1200:nt,
%     
%     timetmp=xrange(t);
% %MM1=10*log10(1*10^-3); %1 mg/l
% waitbar(t/nt,hh,strcat('Iterative Method in progress...--',int2str(t/nt*100),'%'));
% M1=1*10^-3; %1 mg/l
% for i=1:nk
%     
%     spmloop=0.001;
%     deltaspm=100000;
%     iter=0;
% 
%     
%     
%     %IV_water=RL_grid(i,t)-SL_grid(t)+TL_geo(i,t)-Cst_geo+TL_w(i,t);
%     %CMES_water_tmp= 10.^(IV_water/10)/Cst_sdt;  
%     
%     while deltaspm>0.001
%         
%     Ctmp=2*(zeta_vdB+zeta_ddB)*spmloop*deltaR(i)/cos(thetaR)';
%     IVtmp=RL_grid(i,t)-SL_grid(t)+TL_geo(i,t)-Cst_geo+TL_w(i,t)+Ctmp;
%    % CMEStmp(id)=10.^((aiv*IVtmp(id)+biv)/10);
%     CMEStmp=10.^(IVtmp/10)/Cst_sdt;
%     deltaspm=abs(CMEStmp-spmloop);
%     spmloop=CMEStmp;
%     iter=iter+1
%     end
% 
%     
%         CMEStmp_fixe(i,t)=CMEStmp;
%         IVtmp_fixe(i,t)=IVtmp;
%         
%     
%     
% end
% end
% %%



close(hh)
display('End of Iterative Method')

%%

% %%%%% nombre d'iterations : 
% figure; colormap(jet(5));imagesc(itef); axis xy; caxis([0 5]);colorbar; %---> 1-3 iterations
% %%%%% attenuation liee aux particules dans l eau (qs=cumulee ; delta_qs=diff):
% %figure; imagesc(qs); axis xy; caxis([0 2*10^-3]);colorbar;
% for i=2:nk,
% delta_qs(i,1:nt)=qs(i,1:nt)-qs(i-1,1:nt);
% end
% figure; 
% imagesc(xrange,yrange,delta_qs); 
% axis xy; colorbar;caxis([0 0.05]);colorbar; %delta_qs<0.01dB
% %%%%% concentration massique obtenue :
% figure; imagesc(Mgl*1000);axis xy;colorbar;caxis([0 200]);colorbar
% 


figure
imagesc(time_Mgl,yrange,BI)
colorbar
caxis([-70 -20])
set(gca,'YDir','normal')


figure
imagesc(time_Mgl,yrange,Mgl)
colorbar
set(gca,'YDir','normal')
caxis([0 0.1])


figure
hold on
plot(time_Mgl,Mgl(6,:),'r')
plot(time_Mgl,Mgl(8,:),'k')
plot(time_Mgl,Mgl(10,:),'b')
datetickzoom

