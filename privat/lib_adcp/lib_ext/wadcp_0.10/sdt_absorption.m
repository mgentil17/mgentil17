
%% WADCP PACKAGE 0.3 -- function [zeta_vdB,zeta_ddB,Cst_sdt]=sdt_absorption
%%
%
% Backscatter signal processing toolbox
%
% _Written by C. Tessier, IFREMER |contact : caroline.tessier@ifremer.fr 

function [zeta_vdB,zeta_ddB,Cst_sdt]=sdt_absorption(f,rm)
global f  nf
%%

% * sediment and water properties to evaluate sigbar and sigbartot * 

% a_s=40;% default particules diameter [µm]
rm=rm/2*10^-6;      %radius en meters needed !! 
rho_s=2650;         % default individual particles density [kg/m^3]
c_s=4500;           % default sound celerity in particles [m/s]
rho_w=1024;          % default water density [kg/m^3]
c0=1505;            % default sound celerity in water [m/s]

% fractal : 
%nf=2.0; % 2.4 caudebec
dp=0.000004;

rho_s=rho_w+(rho_s-rho_w)*(dp/(rm*2))^(3-nf);

% rho_s
% pause

%modèle de porosité particules élémentaires
% rb=rm-3*10^-6;
% rc=rm+3*10^-6;
% u_g=2*rb;
% rho_w=1024;
% rho_g=2700;
% N=1-0.63*(rm^3/rc^3);
% rho_s=N*rho_w+(1-N)*rho_g;
% Kw=2.25*10^9;
% Kg=1.47*10^10;
% c_s=sqrt(Kw*Kg/(rho_s*(N*Kg+(1-N)*Kw)));
% mu0=2*10^9;
% u0=1000*10^-6;
% c_s=sqrt(c0^2+(mu0/rho_s)*(u_g/u0)^(1/3));

%théorie fractale
% rb=rm-3*10^-6;
% % rc=rm+3*10^-6;
% u_g=2*rb;
% rho_w=1024;
% rho_g=2700;
% rho_s=(rho_g-rho_w)*(4*10^-6/(2*rm));
% N=(rho_s-rho_g)/(rho_w-rho_g);
% Kw=2.25*10^9;
% Kg=1.47*10^10;
% c0=sqrt(Kw*Kg/(rho_s*(N*Kg+(1-N)*Kw)));
% mu0=2*10^9;
% u0=1000*10^-6;
% c_s=sqrt(c0^2+(mu0/rho_s)*(u_g/u0)^(1/3));
% c_s=4500;

% prompt = {'a_s : particules mean/equivalent diameter [µm]';'rho_s : particles density [kg/m^3]';'c_s : sound celerity in particles';'rho0 : water density [kg/m^3]';'c0 : sound celerity in water [m/s]'};
% dlg_title = 'Summary of particles and water characteristics - can be changed by users';
% num_lines = 1;
% def = {num2str(a_s);num2str(rho_s);num2str(c_s);num2str(rho0);num2str(c0) };
% tmp = inputdlg(prompt,dlg_title,num_lines,def);
% a_s=eval(char(tmp{1}));
% rho_s=eval(char(tmp{2}));
% c_s=eval(char(tmp(3)));
% rho0=eval(char(tmp{4}));
% c0=eval(char(tmp{5}));

%-----------------------------------
% rm=rm/2*10^-6;      %radius en meters needed!
v_s=4/3*pi*(rm).^3; %individual volume en m^3
%-----------------------------------

% ask user to choose nature of particles  (2 default case and one user defined) 
% qst_proc=questdlg('Which kind of particles ?','Particles Caracteristics','Mineral','Zooplankton','User defined','Mineral');
% switch qst_proc
%     case 'Mineral'
%         id_sed=1;
%     case 'Zooplankton'
%         id_sed=2;
%     case 'User defined'
%         id_sed=3;
% end
% 
% % ask user to choose theorical model (1) Thorne or (2) Tessier (from Stanton)
% qst_proc=questdlg('Which Acoustic Scattering Model ?','Acoustic Scattering Model',...
%     'Tessier 2006 (all particles)','Thorne 2002 (sand)','Tessier 2006 (all particles');
% switch qst_proc
%     case 'Thorne 2002 (sand)'
%         id_mod=1;
%     case 'Tessier 2006 (all particles)'
%         id_mod=2;
% end
% 
% if ((id_mod==1) & (id_sed~=1)) 
%     display('!!! WARNING !!! ')
%     display('Acoustic Scattering Model of Thorne valid for MINERAL particles only')
% end

% evaluation of scattering cross-sections 
[sigbar,sigbartot]=sigmabar(f,rm,rho_s,c_s,2);%2 stanton
% 
rho_sm=rho_s-rho_w;

% sediment caracterisation in the Scatter Index IV : IV=10log10(M*sigbar/rho_s/vs) 
%Cst_sdt=10*log10(sigbar/rho_s/v_s); % version ctessier
Cst_sdt=(sigbar/(rho_sm*v_s));
%%
% * sediment absorption constantes * 

%-absorption visqueuse particules zeta_v [Urick,1948]
lambda=1500/f/1000;  % mean wave length 
omega=2*pi*f*1000;  %vitesse angulaire (rad/s)              
nu=1.3*10^-6;%viscosite cinematique eau nu (m^2/s)
beta=(omega/2/nu)^0.5;
teta=0.5*(1+9./(2*beta.*rm));
s=(9./(4*beta.*rm)).*(1+1./(beta.*rm));
gg=rho_s/rho_w;
zeta_v=(2*pi/lambda*(gg-1)^2)/(2*rho_s)*(s/(s.^2+(gg+teta).^2));
zeta_vdB=20*log10(exp(1))*zeta_v;

%-absorption due a diffusion par particules zeta_d
zeta_ddB=20/2*log10(exp(1))*(sigbartot/rho_sm/v_s);

%-amortissement total du au sediment alpha_s
%qs=2*(zeta_vdB+zeta_ddB)*M*dR;


