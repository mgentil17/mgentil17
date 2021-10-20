function [nav, out_nav] = filt_nav(nav,mag_declination,lat_O,lon_O)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Cleaning glider navigation data
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

no_data = 1:length(nav.data(:,1));

%==========================================================================
% Magnetic Declination correction
%==========================================================================
nav.data(:,12) = nav.data(:,12) + mag_declination;          % heading

%==========================================================================
% Remove Coordinates Outliers
%==========================================================================
i_bad = find(nav.data(:,2)<lat_O(1) | nav.data(:,2)>lat_O(2) ...
    | nav.data(:,4)<lat_O(1) | nav.data(:,4)>lat_O(2) ...
    | nav.data(:,6)<lat_O(1) | nav.data(:,6)>lat_O(2) ...
    | nav.data(:,3)<lon_O(1) | nav.data(:,3)>lon_O(2) ...
    | nav.data(:,5)<lon_O(1) | nav.data(:,5)>lon_O(2) ...
    | nav.data(:,7)<lon_O(1) | nav.data(:,7)>lon_O(2));

fprintf('##### %d outliers latitude and longitude are removed\n #####',length(i_bad));
nav.data(i_bad,:) = [];

%==========================================================================
% Kill NaN
%==========================================================================
ii_bad = find(isnan(nav.data(:,1)) | isnan(nav.data(:,6)) | ...
    isnan(nav.data(:,7)) | isnan(nav.data(:,14)));

fprintf('##### %d Remove NaN values\n #####',length(ii_bad));
nav.data(ii_bad,:) = [];

out_nav = length(i_bad) + length(ii_bad);

end