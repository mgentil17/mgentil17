%
%/////////////////////////////////////////////////////////////////////////%
% USER DEFINE PARAMETERS
% To update Glider_Adcp_main_program.m
%
% Note : Parameters to be adjusted according to the glider deployment and
% the scientific processes studied.
% For the different steps please refer to the documentation.
%
% version 1.0
% 2021
% Mathieu Gentil / CEFREM / mathieu.gentil@univ-perp.fr
% from Gaël Many and Pierre Cauchy work.
%/////////////////////////////////////////////////////////////////////////%

%%
%==========================================================================
% DEPLOYMENT FEATURES
%==========================================================================
glid_l = {'theque'};                                                       % glider name
deploy = {'matugli_2017'};                                                 % name pf deployment
deploy_year = {'2017'};                                                    % year of deployment


%%
%==========================================================================
% STEP 0 and 7 (INPUTS & OUTPUTS)
%==========================================================================
% Define the path of the toolbox
disp('Select path of the toolbox');
selpath = uigetdir;

% Access to the Toolbox and to the data
addpath(genpath([selpath '/Glider_ADCP_Toolbox_v0']));

% Folder adress to LOAD Glider raw data (CTD, Navigation)
% Be careful, the files must be matrices (.mat)
adressdbd = [selpath '/Glider_ADCP_Toolbox_v0/glider_data/Raw/dbd'];     % raw glider navigation
adressebd = [selpath '/Glider_ADCP_Toolbox_v0/glider_data/Raw/ebd'];     % raw glider ctd

% Folder adress to LOAD ADCP raw data 
% Be careful, the files must be .PD0
adressadcp = [selpath '/Glider_ADCP_Toolbox_v0/glider_data/Raw/PD0'];

% Folder adress to LOAD Suspended Particulates Matter (SPM) data calibration
adressSPM_calib = [selpath '/Glider_ADCP_Toolbox_v0/glider_data/Raw/SPM'];

% Path to save processed data
dpath = [selpath '/Glider_ADCP_Toolbox_v0/glider_data/Processed'];

% Configuration of ADCP used
% If more than one configuration of the adcp is used during the deployment
% indicate the reference
WN_ref = 40;                % number of bin
BinSize_ref = 1;            % size of each bin [m]
sal_val = 38;               % set the salinity (PSU) values for the velocity sound 


%%
%==========================================================================
% STEP 1 (CLEANING)
%==========================================================================
%..........................................................................
% Conversion factor for glider data
%..........................................................................  
% Time conversion
sinceEden = datenum(1970,1,2) - datenum(0000,1,1);
SEC = 24*60*60;
% Radian to degree
r2d = 57.2958;  
% Science unit
CX = 10;                    % factor 10 for conductivity ratio and pressure

%..........................................................................
% Cleaning filters for glider data
%..........................................................................
% Threshold gap or constant pressure
threshold_gap_PP = 1e-3;                       % [db]

% Magnetic declination variation
% Estimate a mean value of magnetic declination  : 1.50 deg E
% (uncertainties ~ 0.33 deg)
% Source : https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#declination
mag_declination = 1.50;  

% Coordinate outliers
lat_O = [-90 90];                              % [deg]
lon_O = [-180 180];                            % [deg]

%..........................................................................
% Cleaning ADCP data
%..........................................................................
time_lag_trick = 0;             % Deal with bad TE / TP config. 0/1
ptch_offset    = 11;			% Hardware tilt. +11° pitch.
Kc             = 0.61;			% Conversion factor from counts do dB.
cnt_tresh      = 64;			% Empirical correlation treshold
vel_tresh      = 0.75;          % Velocity maximal
% db_tresh       = 75;            % Decibel treshold (Todd et al 2016; Merckelbach et al 2018)
ofs_BT         = 0.15;          % Empirical bottom track offset limit 


%%
%==========================================================================
% STEP 2 (SYNCHRONIZATION & INTERPOLATION)
%==========================================================================
% Define the reference time (start/end)
time_sync_start = 7.367356855447916e+05;            % section 6
time_sync_end = 7.367429545262731e+05;              % section 8                

% Define the time step for interpolation
int_timestep = 10/86400;                            % 10 [s];

% Number of points used to smooth data
nb_smooth = 1;


%%
%==========================================================================
% STEP 4/5 (DERIVED PARAMETERS/ABSOLUTE VELOCITIES COMPARISONS)
%==========================================================================
%..........................................................................
% Derived Parameters : Physical and Bio-Optical
%..........................................................................
% Optical units
% OU = 1;                                             % [NTU]
OU = 0;                                             % [m^-^1]

% Reference pressure
PR = 0;                                              % [db]

% Treshold on descent/ascent glider speed
% to index the glider casts
% Empirical value
thresh_cast = 0.15;                                                          % if navigation sensor is choosed to depth alignment                              
% thresh_cast = 0.1;                                                         % if science sensor is choosed to depth alignment  

% filter on depth vaiable to index the glider casts
% Empirical value
filt_order = 10;                                     % order of median

% Thermal lag correction between pairs of profiles
% Very time consumming
prof = 100;                                           % pairs of 50 downcasts/upcasts

%..........................................................................
% Derived Parameters: Glider Flight Model
%..........................................................................
% Temporal period for the glider parameter optimization
idx_time = 1*86400;             % optimization every day

% Pressure level above which data is replaced by NaN for the vertical
% glider velocities (Merckelbach et al., 2009)
threshold_P_W = 10; 

% Glider parameters that will be optimized
% First parameter of param0 : glider volume (Vg, cm3)
% Second parameter of param0 : hull compressibility (eps, Pa-1)
% Third paramater of param0 : parasite drag (Cd0, rad-2)
Vg  = 0.0649;
eps = 6.1e-10;
Cd0 = 0.1;
mg  = 66.66;
param0 = [Vg*1e2,eps*1e10,Cd0*1e2,mg*1e-1];

% Define the option for the optimization to compute flight model
%**************************************************************************
% % Help:
% % First try with default values [1];
%     optimset('PlotFcns',@optimplotfval);
% % if the results are not satisfactory, play on Tolerance parameters [2];
%     options = optimset('TolFun',1e-2,'TolX',1e-2,'PlotFcns',@optimplotfval);
% % Finally if needed, you can play on the iteration number [3];
%     options = optimset('TolFun',1e-2,'TolX',1e-2,'PlotFcns',@optimplotfval,'MaxIter',100);
% % For futher tuning, see fminsearch help. 
%**************************************************************************
% Config choosed for the deployment [2]:
options = optimset('TolFun',1e-2,'TolX',1e-2,'PlotFcns',@optimplotfval);

% If you want check the result of flight model with the initial
% conditions (choose 1, otherwise 0)
initial_comp = 1;

%..........................................................................
% Derived Parameters: Absolute Velocities
%..........................................................................
% Size of stack and moving window 'nwin'
% Be careful it depends on the ADCP sampling strategy
SizeStack = 2*BinSize_ref;                           % Size of stacking layer
if SizeStack<2*BinSize_ref; f=msgbox('Warning : BinStack should be >= 2*binsize'); end
nwin = 2*BinSize_ref;                                % Size of half window

% Stacking of profile data per user-defined layers 
BinStack = 2*BinSize_ref;

% Altimeter value
alt_lim = 3;                                         % the glider goes up 3 m from the bottom (altimeter)

% The number of function evaluations
stddev = 0.066;                                      % Standart deviation (in m s-1) of single ping measurements for 1 m cell size at 600 KHz 

% Monte Carlo simulations
disp('Warning: Monte-Carlo iterations must be > 1');
disp('Warning: Monte-Carlo simulation is a long process');
nMC = 3;

% Smooth velocity data with smooth2a function                                    
nwinc = 2;                                          % number of points to smooth columns (profiles on horizontal)
nwinr = 4;                                          % number of points used to smooth rows (bins on vertical)

% Temporal interpolation of absolute velocity data
dt = 9000/86400;                                                            % ~ 2,5 h

% Maximum depth of the deployment
max_dpth = 200;                                      % Max depth of the Gulf of Lion's shelf


%%
%==========================================================================
% STEP 6 (MAPPING AND GRIDDING)
%==========================================================================
% Reference positions
startlat = 43.333226; startlon = 4.844482;                             % mouth of the Rhône river
endlat = 42.903569; endlon = 4.848985;                                 % limit Rhône continental shelf/canyon (isobath 100m)
proj_angle = 90;                                                       % to determine the projection angle

% Glider distance step interpolation
% Be careful the vertical resolution must be egal or superior of the data
% resolution
diststep = 0.3;                 % 300m
dpthstep = 1;                   % 1m

%%
%==========================================================================
% STEP 7 (PLOT)
%==========================================================================
plot_out = 1;                              % display (1) or not (0) plot

Ymax = 150;                                % depth limit [m]


