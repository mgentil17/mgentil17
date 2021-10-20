%function sigmabar.m
%Compute the individual effective mean Backscattering Cross Section 
% and the global Scattering Cross Section (for the scattering on a particle integrated in all directions)
% !!! here, independent of the vertical axis and written for ONE mean size of particles a_s !!!
%choicePart= 1 : mineral
%            2 : organic (biological o detritical)
%            3 : user defined
%choiceMod=  1 : model of Thorne et al. 2002 (calibrated for sand particles)
%            2 : model Tessier 2006, constructed from Stanton 1998 (fluid sphere, non rigid particles)
% f [kHz]      : nominal frecuency of ADCP 
% a_s [m]      : radius of particle
% rho_s [kg/m3]: density of particle
% c_s [m/s]    : sound velocity in particle
%------------------------------------------------------------------------------

function [sigbar,sigbartot]=sigmabar(f,a_s,rho_s,c_s,choiceMod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0=1505;rho0=1024; %(T=15deg,S=34psu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=2*pi/(c0/f/1000);
ka=k.*a_s;
ka2=ka.^2;
ka4=ka.^4;
%%
% *
% switch choicePart
%     case 1 %mineral---------------------------
%     %c1=3800;%e=16.1551;
%     c1=4500;
%     rho1=2650;
% %     c1=1680;
% %     rho1=1130;
%     g=rho1./rho0;
%     h=c1/c0;      
%     e=g*h^2; 
%     case 2 %zooplankton, Stanton 1998)---------
%    % c1= ;
%     rho1=rho0*1.04;
%     e=1.12; 
%     g=rho1/rho0;
%     h=sqrt(e/g); %h=1.0377;
  %user defined ----------------------
    c1=c_s;
    rho1=rho_s;
    g=rho1./rho0;
    h=c1/c0;      
    e=g*h^2;  
% end

A=(e-1)/(3*e)+(g-1)./(2*g+1); %A=0.57
R=(g*h-1)/(g*h+1);            %R=0.73
B=((e-1)/(3*e))^2+1/3*((g-1)./(2*g+1))^2; %B=0.12



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch (choiceMod)
      
  case 1  %Thorne et al., 1993, 2002  
          %model calibrated for sand particles

Kalpha43=2*B; %(Kalpha=0.18 sables; Kalpha=0.002 zoo)
Kf=2*A;       %(Kf=1.14 sables ; Kf=0.0779 zoo)
C0=1.1*( 1-0.25*exp(-( (ka-1.4)./0.5 ).^2) ).*( 1+0.37*exp(-( (ka-2.8)./2.2 ).^2) ); %for sand mixture
%C0=1.1; %for zoo?

%to compute the SCATTER INDEX :
    %Fm=FORM function
    Fm=C0.*(Kf*ka2)./(1+Kf*ka2); 
    %L=Scatter LENGHT
    L=Fm.*a_s/2;
    %Backscattering cross section
    sigbar=L^2;
%to compute sediment attenuation : 
    %XHI=NORMALISED TOTAL scattering cross section
    XHI=(Kalpha43.*ka4)./(1+1.3*ka2+Kalpha43.*ka4);  
    %sig=TOTAL scattering cross section
    sigbartot= (2*pi*a_s.^2).*XHI;    
    
    case 2   %Tessier 2006, constructed from Stanton, 1998 
            %for non rigid (fluid) spheres (mineral o organic) 
            
%xi=2*sqrt(2)*A/R;
xi=sqrt(2)*A/R;
%sigbar=(A^2).*(a_s.^2).*(ka4./(1+xi.*ka2+(2*(A^2)/R^2)*ka4)); 
sigbar=(A^2).*(a_s.^2).*(ka4./(1+xi.*ka2).^2); 
%asymptotes:
sigR=A^2*(a_s.^2).*ka4; %Rayleigh Backscattering
sigG=R^2*(a_s.^2)/2;    %Geometric Scattering

%xitot=2*sqrt(2*B)/R;
xitot=sqrt(2*B)/R;
sigbartot=4*pi*B*a_s.^2.*(ka4./(1+xitot.*ka2).^2);

 end


