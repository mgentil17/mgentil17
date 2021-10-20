function [i_glidstruct] = thermal_lag(glid,prof,dpath)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Glider Thermal Lag Correction
% According to the work of Garrau et al. (2011)
%
% As this process is very time consuming, it is only run on a subset of 50
% dives when the toolbox goes through an autonomous run. An average is
% performed on alpha and tau parameters, used to correct the thermal lag.
% 
%  References:
%    Garau, B.; Ruiz, S.; G. Zhang, W.; Pascual, A.; Heslop, E.;
%    Kerfoot, J.; and Tintoré, J.; 2011:
%    Thermal Lag Correction on Slocum CTD Glider Data.
%    Journal of Atmospheric and Oceanic Technology, vol. 28, pages 1065-1071.
%
% Functions used are provided by the SOCIB Toolbox and write by Joan Paul
% Beltran (https://github.com/socib/glider_toolbox)
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%==========================================================================
% Define Variables
%==========================================================================
time = glid.sci(:,1);
temp = glid.sci(:,2);                       % temperature
pres = glid.sci(:,3);                       % pressure
cond = glid.sci(:,4);                       % conductivity

%==========================================================================
% Identification of pair casts
%==========================================================================
% Cast variable
cast = glid.cast;

% Remove surfacings
idbad = find(cast==3);
cast(idbad) = [];
time(idbad) = [];
temp(idbad) = [];
pres(idbad) = [];
cond(idbad) = [];

% Difference of casts
dc = diff(cast);

% Index glider inflection between downcast/upcasts
id = find(dc==1 | dc==-1);
idd = id(1:prof);                       % pairs of 50 dives
idd = [1;idd];

%==========================================================================
% Thermal Lag Parameters
%==========================================================================
params = [];
for ii = 1:prof-2                    % 50 dives
    
    % First half dive
    time1 = time(idd(ii):idd(ii+1));
    temp1 = temp(idd(ii):idd(ii+1));
    pres1 = pres(idd(ii):idd(ii+1));
    cond1 = cond(idd(ii):idd(ii+1));
    
    % Second half dive
    time2 = time(idd(ii+1)+1:idd(ii+2));
    temp2 = temp(idd(ii+1)+1:idd(ii+2));
    pres2 = pres(idd(ii+1)+1:idd(ii+2));
    cond2 = cond(idd(ii+1)+1:idd(ii+2));
    
    % Find thermal lag parameters 
    tmp = findThermalLagParams(time1, cond1, temp1, pres1, ...
    time2, cond2, temp2, pres2);
    params = [params;tmp];
end

%==========================================================================
% Thermal Lag Correction
%==========================================================================
% Average parameters for correction
paramsA = mean(params);                 

% Variables to correct
time = glid.sci(:,1);
temp = glid.sci(:,2);                       % temperature
cond = glid.sci(:,4);                       % conductivity
pres = glid.sci(:,3);                       % pressure

% Thermal lag correction
[temp_inside, cond_outside] = correctThermalLag(time, cond, temp, paramsA);

%==========================================================================
% Plot a profile
%==========================================================================
% Colors
pink = rgb('DeepPink');
green = rgb('YellowGreen');

% Plot a profile
h = figure('units','normalized','outerposition',[0 0 1 1]);

subplot 121; plot(temp(idd(5):(idd(6))),-pres(idd(5):(idd(6))),'Color',green,'Linewidth',2);
hold on; grid minor;
plot(temp_inside(idd(5):(idd(6))),-pres(idd(5):(idd(6))),'--','Color',pink,'Linewidth',2);
legend('raw','Corrected');
ylabel('Pressure [db]','Fontweight','bold');
xlabel('Temperature [°C]','Fontweight','bold');
title('Thermal Lag Correction');
ax = gca;
ax.FontSize = 12;

subplot 122; plot(cond(idd(5):(idd(6))),-pres(idd(5):(idd(6))),'Color',green,'Linewidth',2);
hold on; grid minor;
plot(cond_outside(idd(5):(idd(6))),-pres(idd(5):(idd(6))),'--','Color',pink,'Linewidth',2);
legend('raw','Corrected');
xlabel('Conductivity','Fontweight','bold');
title('Thermal Lag Correction');
ax = gca;
ax.FontSize = 12;

cd(dpath);
saveas(h,'thermal_lag_corr.jpeg');

%==========================================================================
% Allocate Outputs
%==========================================================================
i_glidstruct = glid;
i_glidstruct.sci(:,2) = temp_inside;
i_glidstruct.sci(:,4) = cond_outside;

end

