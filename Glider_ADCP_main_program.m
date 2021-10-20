%
%//////////////////////////////////////////////////////////////////////////
% Estimation of currents and turbidity in the Coastal area from Glider
% mounted ADCP
%                                                  
% Main Script
% Version 1.0
% 2021
% Mathieu Gentil / CEFREM / mathieu.gentil@univ-perp.fr
% from Gaël Many and Pierre Cauchy work.
%//////////////////////////////////////////////////////////////////////////


%%
%==========================================================================
% STEP 0:  LOAD AND CONCATENATE DATA
%==========================================================================

%--------------------------------------------------------------------------
% USER PARAMETERS LOADING
Glider_ADCP_define_param
%--------------------------------------------------------------------------

% ebd (ctd sensor) and dbd (navigation) Glider data
[r_glider_nav, r_glider_sci] = dbd_ebd_concat_v3(adressdbd,adressebd);

% Adcp data
[r_adcp_concat] = adcp_concat_v4(adressadcp, WN_ref, BinSize_ref, sal_val);
r_adcp = r_adcp_concat;

% Output : L0 raw and concatenate data
L0_raw_glider_data.nav = r_glider_nav;
L0_raw_glider_data.sci = r_glider_sci;
L0_raw_glider_data.adcp = r_adcp_concat;
cd(dpath);
save ('L0_raw_glider_data.mat','L0_raw_glider_data');


%%
%==========================================================================
% STEP 1: CLEANING GLIDER DATA
%==========================================================================
% Variables
glid_nav = r_glider_nav;
glid_sci = r_glider_sci;
thresh = threshold_gap_PP;

% Units conversion
[nav, sci] = convert_units(glid_nav,glid_sci,sinceEden,SEC,r2d,CX);

% Cleaning glider data
[nav,out_nav] = filt_nav(nav,mag_declination,lat_O,lon_O);
[sci,opt,out_sci,out_opt] = filt_sci(sci,thresh);

% Cleaning ADCP data
[c_curstruct, c_turbstruct] = QA_QC_adcp_v2(r_adcp,time_lag_trick,ptch_offset,...
    Kc,cnt_tresh,vel_tresh,ofs_BT);

if plot_out == 1
    % Plot QA/QC Glider data
    plot_QA_QC(glid_nav,glid_sci,r_adcp,nav,sci,opt,c_curstruct,c_turbstruct,dpath);    
end

% Output : L1 clean glider data
L1_QA_QC_glider_data.nav = nav;
L1_QA_QC_glider_data.sci = sci;
L1_QA_QC_glider_data.opt = opt;
L1_QA_QC_glider_data.adcp_curstruct = c_curstruct;
L1_QA_QC_glider_data.adcp_turbstruct = c_turbstruct;
save ('L1_QA_QC_glider_data.mat','L1_QA_QC_glider_data');


%%
%==========================================================================
% STEP 2: SYNCHRONIZATION AND INTERPOLATION
%==========================================================================
f = msgbox({'Processing:'; 'Only on the small part of the deployment (example)'},'Warn','warn');
    
% Reference time vector 
time_sync = [time_sync_start: int_timestep : time_sync_end];

% Define variables
nav = L1_QA_QC_glider_data.nav;
sci = L1_QA_QC_glider_data.sci;
opt = L1_QA_QC_glider_data.opt;
curstruct = L1_QA_QC_glider_data.adcp_curstruct;
turbstruct = L1_QA_QC_glider_data.adcp_turbstruct;

% Data interpolation and synchronization
[i_glidstruct,i_curstruct,i_turbstruct] = interp_time_series_v3(nav,sci,opt,curstruct,turbstruct,time_sync,nb_smooth);

if plot_out == 1
    % Plot the results of interpolations
    plot_interp(nav,sci,opt,i_glidstruct,dpath);
end

L2_interp_glider_time_series.glid = i_glidstruct;
L2_interp_glider_time_series.adcp_current = i_curstruct;
L2_interp_glider_time_series.adcp_turbidity = i_turbstruct;
cd(dpath);
save ('L2_interp_glider_time_series.mat','L2_interp_glider_time_series');


%%
%==========================================================================
% STEP 3: USER DEFINE SECTIONS
%==========================================================================
% Be careful, this is a manual part to be performed by the user
% depending on the selected deployment
f = msgbox({'User define sections:'; 'Manual part that must be adapted according to the deployment'; 'Check in delim_glid_section function (from line 117)'},'Warn','warn');
pause(5);

promptMessage = sprintf('Do you want to Continue processing,\or Cancel to avort processing?');
button = questdlg(promptMessage, 'Continue', 'Continue', 'Cancel', 'Continue');
if strcmpi(button, 'Cancel')
    return; % Or break or continue
end

% Delimit Glider sections
glid = L2_interp_glider_time_series.glid;
[ind_start_sec, ind_end_sec] = delim_glid_section_v2(glid);


%%
%==========================================================================
% STEP 4 : DERIVED PARAMETERS
%==========================================================================
% Define variables
i_glidstruct = L2_interp_glider_time_series.glid;
i_curstruct = L2_interp_glider_time_series.adcp_current;
i_turbstruct = L2_interp_glider_time_series.adcp_turbidity;

%..........................................................................
% Derived Parameters - 1: Glider Casts & Depth Alignment & Thermal Lag
% Correction
%..........................................................................
[i_glidstruct,i_curstruct,i_turbstruct] = glider_casts(i_glidstruct,i_curstruct,i_turbstruct,thresh_cast,filt_order);

% Thermal Lag
glid = i_glidstruct;
[i_glidstruct] = thermal_lag(glid,prof,dpath);

%..........................................................................
% Derived Parameters - 2: Physical and Bio-Optical
%..........................................................................
cd(adressSPM_calib);
load('calib_mat1617.mat');
calib_SPM = calib_mat1617;
[i_param_derived] = derived_param_v3(i_glidstruct,OU,calib_SPM,PR);

%..........................................................................
% Derived Parameters - 3: Glider Flight Model
%..........................................................................
thresh = threshold_P_W;
in_comp = initial_comp;
[i_GuwM] = glid_underwater_motion_v2(i_glidstruct,i_param_derived,param0,idx_time,thresh,options,in_comp,plot_out,dpath);

%..........................................................................
% Derived Parameters - 4: Absolute Velocities
%..........................................................................
% Flag upcasts and surfacings by NaN
% on each glider sections
ind_start = ind_start_sec;
ind_end = ind_end_sec;
i_param = i_param_derived;
glid = i_glidstruct;
[i_downcasts_sections,i_sections] = flag_downcast_v4(glid, i_curstruct, i_turbstruct, i_GuwM, i_param, ind_start, ind_end, dpath);

dataset = i_downcasts_sections;

% Choose a constrain
fn = {'Bottom Track', 'Depth Average Current'};
[constrain, ~] = listdlg('PromptString',{'Select a constrain.',...
    'Only one constain can be selected at a time.',''},...
    'SelectionMode','single','ListString',fn);

% % % Method to derive absolute velocities from glider flight model: Direct
% % % method
% % c_glider_data = L1_QA_QC_glider_data.nav;                               % clean navigation glider data
% % [dm_vel] = cur_direct_method_MC_v0(dataset,constrain,alt_lim,c_glider_data,SizeStack,BinStack,nwin,nMC,stddev);

% Method to derive absolute velocities from shear stress: Shear Stress
% method
c_glider_data = L1_QA_QC_glider_data.nav;
[sm_vel] = cur_shear_method_MC_v0(dataset,constrain,alt_lim,c_glider_data,SizeStack,BinStack,nwin,nMC,stddev,plot_out);

% Smoothing velocity data
[sm_vel_filt, DAC_cur] = smooth_vel(sm_vel,constrain,max_dpth,nwinc,nwinr,plot_out,dpath);


%%
%==========================================================================
% STEP 5 : COMPARISON OF ABSOLUTE VELOCITIES
%==========================================================================
% % % Plot Absolute Velocities
% % plot_panel_dm_vs_sm_v0(dm_vel_filt,sm_vel_filt,DAC_cur,dpath);

% Comparison of both methods with a semi-external parameter : depth average
% current 
[i_vel, stats_vel] = plot_panel_vel_vs_DAC_v0(sm_vel_filt, DAC_cur, dt, plot_out, dpath);
i_curstruct_vf = sm_vel_filt;

% % % Absolute current method selection
% % disp('#### if results between two methods are similar, prefer shear method: better surface resolution ####');
% % fn = {'Shear Method', 'Direct Method'};
% % [method, ~] = listdlg('PromptString',{'Select a method for absolute velocities.',...
% %     'Only one method can be selected at a time.',''},...
% %     'SelectionMode','single','ListString',fn);
% % 
% % if method == 1
% %     i_curstruct_vf = sm_vel_filt;
% % else
% %     i_curstruct_vf = dm_vel_filt;
% % end

% Output : L3 Derived Parameters
L3_derived_params = i_downcasts_sections;

names = fieldnames(i_curstruct_vf);
for ii = 1:length(names)
    L3_derived_params.(names{ii}).adcp_processed = i_curstruct_vf.(names{ii}); 
end
cd(dpath);
save ('L3_derived_params.mat','L3_derived_params');


%%
%==========================================================================
% STEP 6: MAPPING AND GRIDDING
%==========================================================================
% Glider data projection
glidstruct = i_sections;
[proj_glid] = distglider_v0(glidstruct,startlon,startlat,endlon,endlat,proj_angle,plot_out,dpath);

% Glider data gridding
proj = proj_glid;
[glidmat] = interpdatadist_v0(proj,i_curstruct_vf,i_sections,diststep,dpthstep,alt_lim,deploy);

% Derivation of spatial parameters
[glidmat_pro] = der_spatial_param_v0(glidmat,deploy);

% Output: L4 Mapping
L4_glider_mapping = glidmat_pro;
cd(dpath);
save('L4_glider_mapping.mat','L4_glider_mapping');


%%
if plot_out == 1
%==========================================================================
% STEP 7: PLOT
%==========================================================================
glidmat_pro = L4_glider_mapping;

% Select the variable to plot
fn = {'Temperature','Absolute Salinity','Density Anomaly','Brunt-Vaisala Frequency',...
    'Horizontal Current','E-W Current','N-S Current','Richardson Number',...
    'Optical Backscater','Acoustic Backscatter','SPM concentration','chlorophyll-a'};
[Variable, ~] = listdlg('PromptString',{'Select a variable to plot.',...
    'Only one variable can be selected at a time.',''},...
    'SelectionMode','single','ListString',fn);

% Configuration for plot
[X,Xbis,Y,Z,lat,Col,noS,n] = config_plot(glidmat_pro,Variable);

% Plot glider data
Var = Variable;
plot_glider_time_series(X,Xbis,Y,Ymax,Z,lat,Col,n,Var,glid_l,deploy,dpath);
end

%..........................................................................