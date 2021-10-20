function [proj,h] = dist_glider_v0(DAC,ref_time,dpth,startlon,startlat,endlon,endlat,proj_angle,plot_out,dpath)

%--------------------------------------------------------------------------
% Estimate distance along a projection 
%--------------------------------------------------------------------------

% Pre-allocate Outputs
LG = NaN*ref_time;
LT = NaN*ref_time;
DIST = NaN*ref_time;
DISTproj = NaN*ref_time;

% Define variables
time = DAC(:,1);
lat = DAC(:,2);             % gps latitude
lon = DAC(:,3);             % gps longitude
vx = DAC(:,4);
vy = DAC(:,5);

% Keep only unique position data
[lat, index] = unique(lat,'stable');
lon = lon(index);
vx = vx(index);
vy = vy(index);
time = time(index);

% Deployment caracteristics
nbsurfacing = num2str(length(vx));
display(['Number of FixGPS = ' nbsurfacing])

nbprofile = num2str(length(ref_time));
display(['Number of raw Profile = ' nbprofile])

% Interpolate time FixGPS to time each glider profile
latint = interp1(time,lat,ref_time);
lonint = interp1(time,lon,ref_time);

% Allocate lat and lon glider profile
LG = lonint;
LT = latint;

% Roughly the angle between projection line and current direction
% Current direction has been taken as parallel to the coast
az_proj = azimuth(startlat, startlon, endlat, endlon);

% Maximum spacing in degrees between interpolated positions
maxdiff = .0000085;
% Translates to roughly every 1 meter.
[proj_lat,proj_lon] = interpm([startlat; endlat],[startlon; endlon],maxdiff, 'gc');

if plot_out == 1
    % Plot
    h = figure;
    subplot 211;
    plot(startlon, startlat, 's'); hold on; grid minor;
    xlabel('Longitude (°)','Fontweight','bold');
    ylabel('Latitude (°)','Fontweight','bold');
    plot(endlon, endlat, 's');
    plot(proj_lon, proj_lat,'.');
    plot(lon, lat, '.r');
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    axis equal;
    legend('start point','end point','projection','real position');
end

% gcxgc calculated intersection point of lines starting at given
% lat/lon/azimuth angle pairs
% Projection angle is taken as perpendicular to the user-defined cross
% section
Proj_Angle = proj_angle;
% Proj_Angle = az_proj-proj_angle;

proj_lat(1) = startlat;
proj_lon(1) = startlon;

nlat = NaN(size(LT));
nlon = NaN(size(LG));

tic
wb = waitbar(0,'Projecting profiles ...');
for i = 1:size(LT)
    waitbar(i/length(LT))
    [newlat, newlon] = gcxgc(LT(i),LG(i),Proj_Angle,proj_lat(i),proj_lon(i),az_proj);
    nlat(i) = newlat(2);
    nlon(i) = newlon(2);
end
t=toc;
fprintf('0 of 2- Project data... [%2.2f sec]\n', t);
close(wb)

% % if plot_out == 1
% %     subplot 312;
% %     plot(proj_lon,proj_lat,'Linewidth',2);
% %     xlabel('longitude in degrees');
% %     ylabel('latitude in degrees');
% %     hold on; grid minor;
% %     plot(lon,lat,'r','Linewidth',2);
% %     plot(nlon',nlat','--k','Linewidth',2);
% %     axis equal;
% %     legend('proj lat/lon','lat/lon (real)','nlat/nlon (new)');
% % end

% Allocate dist
LGproj = nlon;
LGproj(isnan(LG)) = NaN;            %apply mask
LTproj = nlat;
LTproj(isnan(LT)) = NaN;            %apply mask

% Projet lat/lon surfacing data
clear nlon nlat
for i=1:length(lat)
    if (~isnan(lat(i)) && lat(i) ~= 0 && ~isnan(lon(i)) && lon(i) ~= 0)
        [newlat, newlon] = gcxgc(lat(i),lon(i),Proj_Angle,proj_lat(i),proj_lon(i),az_proj);
        nlat(i) = newlat(2);
        nlon(i) = newlon(2);
    end
end

LGvproj = nlon;
LTvproj = nlat;

% Calculate distance to Rhone RIVER MOUTH for each surfacing + projected_distance
DISTvx = NaN*lat;
DISTvxproj = NaN*lat;
for i = 1:length(lat)
    DISTvx(i) = deg2km(distance([startlat, startlon], [lat(i), lon(i)]));
    DISTvxproj(i) = deg2km(distance([startlat, startlon], [LTvproj(i), LGvproj(i)]));
end

wb = waitbar(0,'Calculate distances ...');
for i = 1:length(latint)
    waitbar(i/length(latint))
    DIST(i) = deg2km(distance([startlat, startlon], [LT(i), LG(i)]));
    DISTproj(i) = deg2km(distance([startlat, startlon], [LTproj(i), LGproj(i)]));
end
close(wb)

if plot_out == 1
    subplot 212;
%     plot(DIST, 'r','Linewidth',2);
    hold on; grid minor;
    plot(DISTproj, '-b','Linewidth',2);
    xlabel('Number of measures','Fontweight','bold');
    ylabel('Distance (km)','Fontweight','bold');
    legend('Projected Distance', 'Location', 'Best');
%     legend('Real distance', 'Projected Distance', 'Location', 'Best');     
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    cd(dpath);
    saveas(h,'data_projection.jpeg');
end

% Allocate Outputs
DIST(isnan(LGproj))=NaN;
DISTproj(isnan(LGproj))=NaN;
proj.time = time;
proj.lat = LTvproj;
proj.lon = LGvproj;
proj.vx = vx;
proj.vy = vy;
proj.distvx = DISTvx;
proj.distvxproj = DISTvxproj;
proj.latint = LTproj;
proj.lonint = LGproj;
proj.dist = DIST;
proj.distproj = DISTproj;
proj.az_proj = az_proj;

end
