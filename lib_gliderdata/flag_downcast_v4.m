function [i_downcasts_sections, i_sections] = flag_downcast_v4(i_glidstruct,i_curstruct,i_turbstruct,i_GuwM,i_param,ind_start,ind_end,dpath)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% Create a matrix on the deployment with all datasets
% Flag by NaN upcasts and surfacings of glider
% The on-board ADCP being tilted 11° forward of the glider 
% (to compensate for the pitch of the glider during its downcast) 
% renders the instrument unusable during the glider's upcast.
%
% version 2
% 10/2019
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%


%==========================================================================
% DEFINE VARIABLES
%==========================================================================
glid = i_glidstruct;
cur = i_curstruct;
turb = i_turbstruct;
GuwM = i_GuwM;
param = i_param;

%==========================================================================
% DATA PREPARATION
%==========================================================================
% Convert struct2mat
[int_curser, sensor_curser] = struct2mat(cur.ser);
[int_turbser, sensor_turbser] = struct2mat(turb.ser);
[int_GuwM, sensor_GuwM] = struct2mat(GuwM);
[int_param, sensor_param] = struct2mat(param);

% Suppress the optimized parameters (Vg, Cd0, eps and mg) in GuwM
% matrix
int_GuwM(:,end-3:end) = [];
sensor_GuwM(end-3:end) = [];

% Alignment of the matrices on GuwM matrix
% (due to the difference operation required to calculate GuwM)
glid.nav = glid.nav(2:end,:);
glid.sci = glid.sci(2:end,:);
glid.Dac = glid.Dac(2:end,:);
glid.cast = glid.cast(2:end,:);

int_param = int_param(2:end,:);

int_curser = int_curser(2:end,:);
int_turbser = int_turbser(2:end,:);

cur.pro.vel1 = cur.pro.vel1(:,2:end);
cur.pro.vel2 = cur.pro.vel2(:,2:end);
cur.pro.vel3 = cur.pro.vel3(:,2:end);
cur.pro.vel4 = cur.pro.vel4(:,2:end);
cur.pro.Dpth = cur.pro.Dpth(:,2:end);
cur.pro.time = cur.pro.time(:,2:end);

turb.pro.BI1 = turb.pro.BI1(:,2:end);
turb.pro.BI2 = turb.pro.BI2(:,2:end);
turb.pro.BI3 = turb.pro.BI3(:,2:end);
turb.pro.BI4 = turb.pro.BI4(:,2:end);
turb.pro.Dpth1 = turb.pro.Dpth1(:,2:end);
turb.pro.Dpth2 = turb.pro.Dpth2(:,2:end);
turb.pro.Dpth3 = turb.pro.Dpth3(:,2:end);
turb.pro.Dpth4 = turb.pro.Dpth4(:,2:end);

%==========================================================================
% DIVISION OF THE DEPLOYMENT MATRIX INTO SECTION MATRICES
%==========================================================================
for i = 1:length(ind_start)
    num = num2str(i);
    
    ind_s = ind_start(i);                       % section start index
    ind_e = ind_end(i);                         % section end index

    % Division of matrices according to indices
    dataset.nav = glid.nav(ind_s:ind_e,:);
    dataset.Dac = glid.Dac(ind_s:ind_e,:);
    dataset.sci = glid.sci(ind_s:ind_e,:);
    dataset.cast = glid.cast(ind_s:ind_e);
    dataset.param = int_param(ind_s:ind_e,:);
    dataset.curser = int_curser(ind_s:ind_e,:);
    dataset.turbser = int_turbser(ind_s:ind_e,:);
    dataset.Guwm = int_GuwM(ind_s:ind_e,:);
    dataset.curpro.vel1 = cur.pro.vel1(:,(ind_s:ind_e)');
    dataset.curpro.vel2 = cur.pro.vel2(:,(ind_s:ind_e)');
    dataset.curpro.vel3 = cur.pro.vel3(:,(ind_s:ind_e)');
    dataset.curpro.vel4 = cur.pro.vel4(:,(ind_s:ind_e)');
    dataset.curpro.Dpth = cur.pro.Dpth(:,(ind_s:ind_e)');
    dataset.curpro.time = cur.pro.time(:,(ind_s:ind_e)');
    dataset.turbpro.BI1 = turb.pro.BI1(:,(ind_s:ind_e)');
    dataset.turbpro.BI2 = turb.pro.BI2(:,(ind_s:ind_e)');    
    dataset.turbpro.BI3 = turb.pro.BI3(:,(ind_s:ind_e)');    
    dataset.turbpro.BI4 = turb.pro.BI4(:,(ind_s:ind_e)');    
    dataset.turbpro.Dpth1 = turb.pro.Dpth1(:,(ind_s:ind_e)');    
    dataset.turbpro.Dpth2 = turb.pro.Dpth2(:,(ind_s:ind_e)');    
    dataset.turbpro.Dpth3 = turb.pro.Dpth3(:,(ind_s:ind_e)'); 
    dataset.turbpro.Dpth4 = turb.pro.Dpth4(:,(ind_s:ind_e)'); 
    
    % Store sections in the common variables
    eval(['allcasts.section' num '=dataset;']);
    
    % Flag surfacings and upcasts by NaN
    ind_up = find(dataset.cast==2);
    ind_surf = find(dataset.cast==3);

    dataset.nav(ind_up,:) = NaN;       dataset.nav(ind_surf,:) = NaN;
    dataset.Dac(ind_up,:) = NaN;       dataset.Dac(ind_surf,:) = NaN;
    dataset.sci(ind_up,:) = NaN;       dataset.sci(ind_surf,:) = NaN;
    dataset.cast(ind_up,:) = NaN;      dataset.cast(ind_surf,:) = NaN;
    dataset.curser(ind_up,:) = NaN;    dataset.curser(ind_surf,:) = NaN;
    dataset.Guwm(ind_up,:) = NaN;      dataset.Guwm(ind_surf,:) = NaN;

    dataset.curpro.vel1(:,ind_up') = NaN;     dataset.curpro.vel1(:,ind_surf') = NaN;
    dataset.curpro.vel2(:,ind_up') = NaN;     dataset.curpro.vel2(:,ind_surf') = NaN;
    dataset.curpro.vel3(:,ind_up') = NaN;     dataset.curpro.vel3(:,ind_surf') = NaN;
    dataset.curpro.vel4(:,ind_up') = NaN;     dataset.curpro.vel4(:,ind_surf') = NaN;
    dataset.curpro.Dpth(:,ind_up') = NaN;     dataset.curpro.Dpth(:,ind_surf') = NaN;
    dataset.curpro.time(:,ind_up') = NaN;     dataset.curpro.time(:,ind_surf') = NaN;

    dataset.cfg_adcp = cur.cfg;

    dataset.sensors.curser = sensor_curser;
    dataset.sensors.Guwm = sensor_GuwM;
    dataset.sensors.param = sensor_param;
    
    % Store sections in the common variables
    eval(['downcasts.section' num '=dataset;']);   
end

%==========================================================================
% Allocate Outputs
%==========================================================================
i_downcasts_sections = downcasts;
i_sections = allcasts;

end