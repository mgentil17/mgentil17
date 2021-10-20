function [i_glidstruct,i_curstruct,i_turbstruct] = glider_casts(i_glidstruct,i_curstruct,i_turbstruct,thresh_cast,filt_order)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Depths Comparisons and Alignment
% Glider Casts Indexing
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%==========================================================================
% Define Variables
%==========================================================================
cur = i_curstruct;          % current
turb = i_turbstruct;        % turbidity
glid = i_glidstruct;        % glider

time = glid.sci(:,1);
lat = glid.nav(:,2);
P = glid.sci(:,3);          % pressure
dpth_nav = glid.nav(:,8);
dpth_adcp = cur.ser.dpth;

%==========================================================================
% Depth Sensors Comparisons
%==========================================================================
% Derived depth of science sensor
dpth_sci = sw_dpth(P,lat);

% Colors
pink = rgb('DeepPink');
green = rgb('YellowGreen');
blue = rgb('SteelBlue');

% Plot
figure('units','normalized','outerposition',[0 0 1 1]);
plot(time,-dpth_nav,'o','Color',blue,'MarkerSize',5);
hold on; grid minor;
plot(time,-dpth_sci,'--','Color',pink,'Linewidth',2);
plot(time,-dpth_adcp,'-','Color',green,'Linewidth',2);
legend('nav_d_p_t_h', 'ctd_d_p_t_h', 'adcp_d_p_t_h');
ylabel('Depth [m]','Fontweight','bold');
xlim([time(1) time(1000)]);
ax = gca;
ax.FontSize = 12;
set(ax,'TickDir','out');
datetick('x','dd HH:MM','keepticks','keeplimits');

% Compute RMSE between science and navigation sensors
RMSE_dpth = sqrt(nanmean(dpth_nav - dpth_sci).^2);                         % Root Mean Squared Error
fprintf('###### %d m: RMSE between the pressure sensors\n #######',RMSE_dpth);

if RMSE_dpth < 0.5
    % Alignment of science depth sensor on navigation depths sensors
    disp('##### Depth alignment on the navigation sensor #####');
    disp('##### Correction of science pressure #####');
    glid.sci(:,3) = sw_pres(dpth_nav,lat);
else
    disp('##### Error : you have to apply an offset #####');
end

%==========================================================================
% Casts Index
% flag 1 : downcasts
% flag 2 : upcasts
% flag 3 : surfacings
%==========================================================================
depth = dpth_nav;
[cast] = index_cast_v2(depth,time,thresh_cast,filt_order);

% Remove the last cell of matrixes to have the same size as the cast vector
glid.Dac(end,:) = [];
glid.nav(end,:) = [];
glid.sci(end,:) = [];

cur.ser.BTvel1(end,:) = [];
cur.ser.BTvel2(end,:) = [];
cur.ser.BTvel3(end,:) = [];
cur.ser.BTvel4(end,:) = [];
cur.ser.sal(end,:) = [];     turb.ser.sal(end,:) = [];
cur.ser.temp(end,:) = [];    turb.ser.temp(end,:) = [];
cur.ser.dpth(end,:) = [];    turb.ser.dpth(end,:) = [];
cur.ser.head(end,:) = [];    turb.ser.head(end,:) = [];
cur.ser.ptch(end,:) = [];    turb.ser.ptch(end,:) = [];
cur.ser.roll(end,:) = [];    turb.ser.roll(end,:) = [];
cur.ser.time(end,:) = [];    turb.ser.time(end,:) = [];

cur.ser.alt(end,:) = [];   
turb.ser.altB1(end,:) = [];
turb.ser.altB2(end,:) = [];
turb.ser.altB3(end,:) = [];
turb.ser.altB4(end,:) = [];

cur.pro.vel1(:,end) = [];    turb.pro.BI1(:,end) = [];
cur.pro.vel2(:,end) = [];    turb.pro.BI2(:,end) = [];
cur.pro.vel3(:,end) = [];    turb.pro.BI3(:,end) = [];
cur.pro.vel4(:,end) = [];    turb.pro.BI4(:,end) = [];

cur.pro.time(:,end) = [];  

cur.pro.Dpth(:,end) = [];   
turb.pro.Dpth1(:,end) = [];
turb.pro.Dpth2(:,end) = [];
turb.pro.Dpth3(:,end) = [];
turb.pro.Dpth4(:,end) = [];

time(end) = [];
dpth_nav(end) = [];
dpth_adcp(end) = [];
depth(end) = [];

%==========================================================================
% ADCP Depth
%==========================================================================
% Apply an offset on the depth ADCP sensor
dpth_offset = dpth_nav - dpth_adcp;

ind_downcast = find(cast.data==1);                                         % find index of downcasts
ind_upcast = find(cast.data==2);                                           % find index of upcasts
ind_surfcast = find(cast.data==3);                                         % find index of surfcasts

% On depth current series
cur.ser.dpth(ind_downcast) = dpth_adcp(ind_downcast) + dpth_offset(ind_downcast);
cur.ser.dpth(ind_upcast) = dpth_adcp(ind_upcast) + dpth_offset(ind_upcast);
cur.ser.dpth(ind_surfcast) = dpth_adcp(ind_surfcast) + dpth_offset(ind_surfcast);

% Plot
figure;
plot(time,-depth,'o','Color',pink,'MarkerSize',5);
hold on;
plot(time,-cur.ser.dpth,'-','Color',green,'Linewidth',2);
legend('sci_d_p_t_h', 'adcp_d_p_t_h + offset');
ylabel('Depth [m]','Fontweight','bold');
xlim([time(1) time(1000)]);
ax = gca;
ax.FontSize = 12;
set(ax,'TickDir','out');
datetick('x','dd HH:MM','keepticks','keeplimits');

% On depth profiles
for i = 1:length(dpth_offset')
    cur.pro.Dpth(:,i) = cur.pro.Dpth(:,i) + dpth_offset(i);
    
    % Apply the correction on turbidity  profiles
    turb.pro.Dpth1(:,i) = turb.pro.Dpth1(:,i) + dpth_offset(i);
    turb.pro.Dpth2(:,i) = turb.pro.Dpth2(:,i) + dpth_offset(i);
    turb.pro.Dpth3(:,i) = turb.pro.Dpth3(:,i) + dpth_offset(i);
    turb.pro.Dpth4(:,i) = turb.pro.Dpth4(:,i) + dpth_offset(i);
end

% Apply the correction on turbidity serie 
turb.ser.dpth = cur.ser.dpth;

%==========================================================================
% Allocate Outputs
%==========================================================================
i_glidstruct = glid;
i_glidstruct.cast = cast.data;
i_curstruct = cur;
i_turbstruct = turb;

end