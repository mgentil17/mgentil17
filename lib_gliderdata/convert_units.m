function [nav, sci] = convert_units(glid_nav,glid_sci,sinceEden,SEC,r2d,CX)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Unit Conversion for glider data (science and navigation
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%==========================================================================
% Time units conversion
%==========================================================================
% Navigation data
fromEden = glid_nav.data(:,1)./SEC;
glid_nav.data(:,1) = sinceEden + fromEden;

% Science data
fromEden = glid_sci.data(:,1)./SEC;
glid_sci.data(:,1) = sinceEden + fromEden;

%==========================================================================
% Position units conversion
%==========================================================================
% Fix GPS
glid_nav.data(:,2) = floor(((glid_nav.data(:,2)/100-floor(glid_nav.data(:,2)/100)) ...
    *(10000/60))*100)/10000 + floor(glid_nav.data(:,2)/100);               % ltgps
glid_nav.data(:,3) = floor(((glid_nav.data(:,3)/100-floor(glid_nav.data(:,3)/100)) ...
    *(10000/60))*100)/10000 + floor(glid_nav.data(:,3)/100);               % lggps

% Waypoints
glid_nav.data(:,4) = floor(((glid_nav.data(:,4)/100-floor(glid_nav.data(:,4)/100)) ...
    *(10000/60))*100)/10000 + floor(glid_nav.data(:,4)/100);               % ltwpt
glid_nav.data(:,5) = floor(((glid_nav.data(:,5)/100-floor(glid_nav.data(:,5)/100)) ...
    *(10000/60))*100)/10000 + floor(glid_nav.data(:,5)/100);               % lgwpt

% Dead-reckoning
glid_nav.data(:,6) = floor(((glid_nav.data(:,6)/100-floor(glid_nav.data(:,6)/100)) ...
    *(10000/60))*100)/10000 + floor(glid_nav.data(:,6)/100);               % lt
glid_nav.data(:,7) = floor(((glid_nav.data(:,7)/100-floor(glid_nav.data(:,7)/100)) ...
    *(10000/60))*100)/10000 + floor(glid_nav.data(:,7)/100);               % lg

%==========================================================================
% Attitude units conversion
%==========================================================================
% Radian to degrees
glid_nav.data(:,11) = r2d * glid_nav.data(:,11);                           % pitch
glid_nav.data(:,12) = r2d * glid_nav.data(:,12);                           % heading
glid_nav.data(:,13) = r2d * glid_nav.data(:,13);                           % roll

%==========================================================================
% Scence units conversion
%==========================================================================
% Conductivity ratio
glid_sci.data(:,4) = glid_sci.data(:,4).*CX; 

% Pressure
glid_sci.data(:,3) = glid_sci.data(:,3).*CX; 

%==========================================================================
% Allocate Outputs
%==========================================================================
nav = glid_nav;
sci = glid_sci;

end