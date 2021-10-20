function [r_glider_nav, r_glider_sci] = dbd_ebd_concat_v3(adressdbd,adressebd)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% Select and concatenate ebd/dbd data in a structure
% to obtain a glider and science matrix for a whole deployment
%
% version 3.0
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%


%==========================================================================
%% GLIDER NAVIGATION (DBD DATA)
%==========================================================================
dir_racine_dbd = adressdbd;                                                % dbd file path
cd(dir_racine_dbd);
file_dbd = dir('*.mat');                                                   % number of .mat files

% LOAD DATA
for ii = 1:length(file_dbd)
    names_dbd = file_dbd(ii).name;                                         % name of files
    F = fullfile(dir_racine_dbd,names_dbd);                                % path of file
    Struct = load(F);                                                      % load file in a structure
    names_Struct = fieldnames(Struct);
    
    % Create a structure
    eval(['dbd_sensors' '=Struct.' names_Struct{end} '.sensors;']);        % sensors
    eval(['dbd_data' '=Struct.' names_Struct{end} '.data;']);              % data
    
    % EXTRACTION OF THE DESIRED VARIABLES
    % Index identification
    ind_nav_time = dbd_sensors.m_present_time;
    ind_nav_ltgps = dbd_sensors.m_gps_lat;                                 % gps latitude
    ind_nav_lggps = dbd_sensors.m_gps_lon;                                 % gps longitude
    ind_nav_ltwpt = dbd_sensors.c_wpt_lat;                                 % latitude waypoint (glider command)
    ind_nav_lgwpt = dbd_sensors.c_wpt_lon;                                 % longitude waypoint (glider command)
    ind_nav_lt = dbd_sensors.m_lat;                                        % latitude (computed by the glider) dead reckoning
    ind_nav_lg = dbd_sensors.m_lon;                                        % longitude (computed by the glider)
    ind_nav_vx = dbd_sensors.m_water_vx;                                   % speed (east-west) between two surfacings
    ind_nav_vy = dbd_sensors.m_water_vy;                                   % speed (north-south) between two surfacings
    ind_nav_balla = dbd_sensors.m_ballast_pumped;                          % glider ballasting
    ind_nav_pitch = dbd_sensors.m_pitch;                                   % pitch of glider
    ind_nav_head = dbd_sensors.m_heading;                                  % heading of glider
    ind_nav_roll = dbd_sensors.m_roll;                                     % roll of glider
    ind_nav_pp = dbd_sensors.m_depth;                                      % depth of glider
    
    % Index of final water
    % to compare with water velocity
    ind_nav_vx_fin = dbd_sensors.m_final_water_vx;
    ind_nav_vy_fin = dbd_sensors.m_final_water_vy;
        
    % Select data corresponding to sensors
    ind_nav_sensors = [ind_nav_time ind_nav_ltgps ind_nav_lggps ind_nav_ltwpt ind_nav_lgwpt ind_nav_lt ...
        ind_nav_lg ind_nav_vx ind_nav_vy ind_nav_balla ind_nav_pitch ind_nav_head ind_nav_roll ind_nav_pp ind_nav_vx_fin ind_nav_vy_fin];
    var_nav = {'time','ltgps','lggps','ltwpt','lgwpt','lt','lg','vx','vy','balla','pitch','head','roll','depth','vx_final','vy_final'};
    
    r_nav_data = dbd_data(:,ind_nav_sensors);
    
    % Stores navigation variables and data in a structure
    eval(['glider_nav.sensors' '=var_nav;']);
    eval(['glider_nav.data.' names_Struct{end} '=r_nav_data;']);
    close all
end

% CONCATENATE DATA
cell = struct2cell(glider_nav);                                            % convert struct to cell
data_cell = struct2cell(cell{end});                                        % select cell of data
dbd_concat = vertcat(data_cell{:});                                        % concatenates all cells
[~,col_nav] = size(glider_nav.sensors);


%==========================================================================
%% GLIDER SCIENCE (EBD DATA)
%==========================================================================
dir_racine_ebd = adressebd;                                                % ebd file path
cd(dir_racine_ebd);
file_ebd = dir('*.mat'); 

% LOAD DATA
for ii = 1:length(file_ebd)
    names_ebd = file_ebd(ii).name;                                         % name of files
    F = fullfile(dir_racine_ebd,names_ebd);                                % path of file
    Struct = load(F);                                                      % load file in a structure
    names_Struct = fieldnames(Struct);
    
    % Load the ebd sensors into a structure
    eval(['ebd_sensors' '=Struct.' names_Struct{end} '.sensors;']);
    eval(['ebd_data' '=Struct.' names_Struct{end} '.data;']);
    
    % Check the number of science sensors
    [~,sensors_number] = size(ebd_data);
    
    if sensors_number > 27                                                 % below 27, there are no variables we are looking for
        % Find the index of the ebd variables studied
        ind_sci_time = ebd_sensors.sci_m_present_time;                     % time
        ind_sci_TT = ebd_sensors.sci_water_temp;                           % temperature
        ind_sci_PP = ebd_sensors.sci_water_pressure;                       % pressure of ctd
        ind_sci_RT = ebd_sensors.sci_water_cond;                           % conductivity ratio
        
        % Identification of optical sensor used
        try
            ind_sci_CHL = ebd_sensors.sci_flbbcd_chlor_units;              % chlorophylle-a
            ind_sci_BB = ebd_sensors.sci_flbbcd_bb_units;                  % backscatter
            ind_sci_CDOM = ebd_sensors.sci_flbbcd_cdom_units;              % colored dissolved organic matter
            disp('FLBBCD SENSOR')
        catch
            disp(['### No flbbcd sensor ###'])
            ind_sci_flbbcd = NaN(size(ind_sci_time));
        end
        
        try
            ind_sci_CHL = ebd_sensors.sci_flntu_chlor_units;               % chlorophylle-a
            ind_sci_BB = ebd_sensors.sci_flntu_turb_units;                 % backscatter
            ind_sci_CDOM = ebd_sensors.sci_flbbcd_cdom_units;              % colored dissolved organic matter
            disp('FLNTU SENSOR')
        catch
            disp(['### No CDOM on flntu sensor ###'])
        end
        
        % Select data corresponding to sensors
        if exist('ind_sci_flbbcd')
            ind_sci_sensors = [ind_sci_time ind_sci_TT ind_sci_PP ind_sci_RT ind_sci_CHL ind_sci_BB];
            var_sci = {'time','TT','PP','RT','CHL','BB'};
            opt_unit = 'NTU';
            
        else
            ind_sci_sensors = [ind_sci_time ind_sci_TT ind_sci_PP ind_sci_RT ind_sci_CHL ind_sci_BB ind_sci_CDOM];
            var_sci = {'time','TT','PP','RT','CHL','BB','CDOM'};
            opt_unit = 'm-1 sr-1';
        end
        
        r_sci_data = ebd_data(:,ind_sci_sensors);
        
        % Stores science variables and data in a structure
        eval(['glider_sci.sensors' '=var_sci;']);
        eval(['glider_sci.data.' names_Struct{end} '=r_sci_data;']);
        close all
        
    else
        disp(['### The science file is corrupted, it is not retained for this study ###'])
    end 
end

% CONCATENATE DATA
cell = struct2cell(glider_sci);                                            % convert struct to cell
data_cell = struct2cell(cell{end});                                        % select cell of data
ebd_concat = vertcat(data_cell{:});                                        % concatenates all cells
[~,col_sci] = size(glider_sci.sensors);


%==========================================================================
%% ALLOCATE OUTPUTS
%==========================================================================
r_glider_nav.sensors = glider_nav.sensors;                               % glider navigation sensors
r_glider_nav.data = dbd_concat;                                          % glider navigation data on the whole deployment

r_glider_sci.sensors = glider_sci.sensors;                               % glider navigation sensors
r_glider_sci.data = ebd_concat;                                          % glider navigation data on the whole deployment

end