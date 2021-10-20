function [i_param_derived] = derived_param_v3(i_glidstruct,OU,calib_SPM,PR)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Derived Physaical and Bio-Optical Parameters
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%==========================================================================
% Define Variables
%==========================================================================
glid = i_glidstruct;
time = glid.sci(:,1);
T = glid.sci(:,2);                          % temperature
cndr = glid.sci(:,4);                       % conductivity ratio
P = glid.sci(:,3);                          % pressure
lon = glid.nav(:,3);                        % longitude
lat = glid.nav(:,2);                        % latitude
BB = glid.sci(:,6);                         % optical backscatter
m = calib_SPM.slope;                        % SPM calibration slope
b = calib_SPM.y_int;                        % SPM calibration y-intercept

%==========================================================================
% Derived Params
%==========================================================================
SP = sw_salt(cndr/sw_c3515,T,P);            % practical salinity
dens = sw_dens(SP,T,P);                     % density
TP = sw_ptmp(SP,T,P,PR);                    % potential temperature
[SA, ~] = gsw_SA_from_SP(SP,P,lon,lat);     % absolute salinity
CT = gsw_CT_from_t(SA,T,P);                 % conservative temperature
sigma0 = gsw_sigma0(SA,CT);                 % potential density anomaly  

% Derivation of optical backscattering
if OU==1
    bbp = BB;
else
    S = nanmean(SP);                                  % mean salinity of the dataset
    lambda = 700;                                     % scatter wavelength of sensor
    beta = BB;                                        % Beta is the scatter at lambda
    [~, bbp] = Beta2Bb_v2(S, lambda, beta);
end

% Derived SPM concentration
SSC = m.*bbp + b;

%==========================================================================
% Allocate Outputs
%==========================================================================
i_param_derived.time = time;
i_param_derived.P = P;
i_param_derived.SP = SP;
i_param_derived.dens = dens;
i_param_derived.TP = TP;
i_param_derived.SA = SA;
i_param_derived.CT = CT;
i_param_derived.sigma0 = sigma0;
i_param_derived.bbp = bbp;
i_param_derived.SSC = SSC;
end