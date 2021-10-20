function [i_vel, stats_vel] = plot_panel_vel_vs_DAC_v0(sm_vel_filt, DAC_cur, dt, plot_out, dpath)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% Comparison of direct and shear method with independant measurements (DAC) in
% order to select an absolute velocity method estimation
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sm = sm_vel_filt;
dac = DAC_cur;
names = fieldnames(sm);                       % number of sections

for ii = 1:length(names)
    
    %======================================================================
    % DEFINE VARIABLES ON USER-DEFINED SECTION
    %======================================================================
    time = sm.(names{ii}).time;
    dpth = sm.(names{ii}).dpth;
  
    % Depth vecor
    max_dpth = (max(dpth));
    [~,I] = max(max_dpth);
    z = dpth(:,I);
    
    % DAC
    lat = dac.(names{ii}).lat;
    lon = dac.(names{ii}).lon;
    vx = dac.(names{ii}).vx;
    vy = dac.(names{ii}).vy;
    
    % Shear method
    u_sm = sm.(names{ii}).u;
    v_sm = sm.(names{ii}).v;
    
    %======================================================================
    % METHODS VS DEPTH AVERAGE CURRENT
    %======================================================================
    %......................................................................
    % Part 1 : Preparation of adcp and glider current data
    %......................................................................
    umean_sm = nanmedian(u_sm,1);                         % shear method
    vmean_sm = nanmedian(v_sm,1);
    
    % Remove NaN of variables
    idx_NaN = find(isnan(vmean_sm) | isnan(umean_sm));
    
    time(idx_NaN) = [];
    umean_sm(idx_NaN) = [];
    vmean_sm(idx_NaN) = [];   
    lat(idx_NaN) = [];
    lon(idx_NaN) = [];
    vx(idx_NaN) = [];
    vy(idx_NaN) = [];
    
    %......................................................................
    % Part 2 : Data interpolation on a reference time
    %......................................................................
    % Reference time
    time_ref = min(time):dt:max(time);
    
    % Data interpolation
    lat_int = interp1(time,lat,time_ref);
    lon_int = interp1(time,lon,time_ref);    
    vx_int = interp1(time,vx,time_ref);
    vy_int = interp1(time,vy,time_ref);
    
    u_sm_int = interp1(time,umean_sm,time_ref);
    v_sm_int = interp1(time,vmean_sm,time_ref);    
    
    if plot_out == 1
        %......................................................................
        % Part 3 : Stickplots
        %......................................................................
        % Colors
        blue = 	[0, 0.4470, 0.7410];
        orange = 	[0.8500, 0.3250, 0.0980];
        
        % Limits
        lw = min(lon_int)-.5; le = max(lon_int)+.5;
        ls = min(lat_int)-.2; ln = max(lat_int)+.2;
        
        h = figure('units','normalized','outerposition',[0 0 1 1]);
        
        grid on;box on; hold on;
        m_proj('mercator','lon', [lw le],'lat',[ls ln]);
        hold on;
        load golCoast.dat;
        h1=m_line(golCoast(:,1),golCoast(:,2),'color','k','linewidth',0.5);
        load griddedBathyGOL;
        col = 0.8* [1 1 1];
        [cs1,h3]=m_contour(XI,YI,ZI,[-120:10:0],'Color',col);
        set(h3,'linestyle','-');
        hold on;
        [cs2,h4]=m_contour(XI,YI,ZI,[-4000:200:0],'Color',col);
        set(h4,'linestyle','-');
        hold on;
        m_plot(lon_int,lat_int,'-k','Markersize',.5);
        hold on;
        h6 = m_quiver(lon_int,lat_int,vx_int,vy_int,'Color',blue,'Linewidth',1.5);
        set(h6,'AutoScale','on', 'AutoScaleFactor', 2)
        hold on;
        m_grid('box','fancy','tickdir','in','linestyle','none');
        hold on;
        h7 = m_quiver(lon_int,lat_int,u_sm_int,v_sm_int,'Color',orange,'Linewidth',1.5);
        set(h7,'AutoScale','on', 'AutoScaleFactor', 2)
        legend([h6,h7],{'Glider DACs','ADCP_S_M DACs'})
        xlabel('Longitude','fontsize',12,'fontweight','bold');
        title('U Shear Method .vs U Depth Average Current');
        
        cd(dpath);
        saveas(h,'Stickplots_DAC_from_methods.jpeg');
    end

    %......................................................................
    % Part 4 : Get the statistics
    % mean, std, rmsd, correlation
    %......................................................................     
    stats_vel{ii}.vx_usm = allstats(vx_int,u_sm_int);
    stats_vel{ii}.vy_vsm = allstats(vy_int,v_sm_int);

    %......................................................................
    % Part 5 : Allocate Outputs
    %......................................................................
    i_vel{ii}.time = time_ref;
    i_vel{ii}.lat = lat_int;
    i_vel{ii}.lon = lon_int;    
    i_vel{ii}.vx = vx_int;
    i_vel{ii}.vy = vy_int;
    i_vel{ii}.u_sm = u_sm_int;
    i_vel{ii}.v_sm = v_sm_int;
    
end

%..........................................................................
% Part 5 : Mean velocities on the deployment (Taylor's diagram)
% Use depth average current (DAC such as reference
%..........................................................................
% Convert cell to matrix
struct_stats = cell2mat(stats_vel);

% Concatenate data
for j = 1:length(struct_stats)
    vx_usm_conc(1:length(struct_stats(j).vx_usm),j) = struct_stats(j).vx_usm(:,2);
    vy_vsm_conc(1:length(struct_stats(j).vy_vsm),j) = struct_stats(j).vy_vsm(:,2); 
    vx_conc(1:length(struct_stats(j).vx_usm),j) = struct_stats(j).vx_usm(:,1); 
    vy_conc(1:length(struct_stats(j).vy_vsm),j) = struct_stats(j).vy_vsm(:,1);   
end

% Calculate an average over the entire deployment
% Standard deviation
STD_vx = nanmean(vx_conc(2,:));
STD_vy = nanmean(vy_conc(2,:));
STD_usm = nanmean(vx_usm_conc(2,:));
STD_vsm = nanmean(vy_vsm_conc(2,:));
STDs = [STD_vx STD_vy STD_usm STD_vsm];

% RMSD
RMS_vx = nanmean(vx_conc(3,:));
RMS_vy = nanmean(vy_conc(3,:));
RMS_usm = nanmean(vx_usm_conc(3,:));
RMS_vsm = nanmean(vy_vsm_conc(3,:));
RMSs = [RMS_vx RMS_vy RMS_usm RMS_vsm];

% Correlation
COR_vx = nanmean(vx_conc(4,:));
COR_vy = nanmean(vy_conc(4,:));
COR_usm = nanmean(vx_usm_conc(4,:));
COR_vsm = nanmean(vy_vsm_conc(4,:));
CORs = [COR_vx COR_vy COR_usm COR_vsm];

% Plot Taylor's diagram
STDsm = [STDs(1:2) STDs(3) STDs(4)];
RMSsm = [RMSs(1:2) RMSs(3) RMSs(4)];
CORsm = [CORs(1:2) CORs(3) CORs(4)];

h1 = figure;

set(h1,'Units','Normalized','Outerposition',[0 0 .6 .6]);

[hp ht axl] = taylordiag_v2(STDsm,RMSsm,CORsm);
title('Shear method');
ax = gca;
ax.FontSize = 12;

saveas(h1,'Taylor_diagram.jpeg');

end
