function [X,Xbis,Y,Z,lat,Col,noS,n] = config_plot(glidmat_pro,Variable)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Configuration to Plot Data
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Glider Variable
g = glidmat_pro;
v = Variable;

% Number of sections
names = fieldnames(g.position);
noS = length(names);

% Get the configuration
for i = 1:length(names)
    
    % Position Variables
    Xtmp = g.position.(names{i}).X;             % distance
    Xbistmp = g.position.(names{i}).time;       % time
    Ytmp = g.position.(names{i}).Y;             % depth
    lattmp = g.position.(names{i}).lat;         % latitude
    
    % Identification of the variable to plot and color to use
    if v == 1
        Ztmp = g.hydrology.(names{i}).temp;            % temperature
        Col = cmocean('thermal');
        n = 'T [°C]';
    elseif v == 2
        Ztmp = g.hydrology.(names{i}).SA;              % absolute salinity
        Col = cmocean('haline');
        n = 'SA [g kg^-^1]';
    elseif v == 3
        Ztmp = g.hydrology.(names{i}).sigma0;          % density anomaly
        Col = cmocean('dense');
        n = '\sigma\theta [kg m^3]';
    elseif v == 4
        Ztmp = g.hydrology.(names{i}).N2;              % Brunt-Vaisala frequency
        Col = cmocean('amp');
        n = 'N^2 [s^-^1]';
    elseif v == 5
        Ztmp = g.hydrodynamique.(names{i}).U;          % horizontal current magnitude
        Col = cmocean('speed');
        n = 'U_m_a_g [m s^-^1]';
    elseif v == 6
        Ztmp = g.hydrodynamique.(names{i}).u;          % east-west current
        Col = cmocean('balance');
        n = 'u [m s^-^1]';
    elseif v == 7
        Ztmp = g.hydrodynamique.(names{i}).v;          % north-south current
        Col = cmocean('balance');
        n = 'v [m s^-^1]';
    elseif v == 8
        Ztmp = g.hydrodynamique.(names{i}).Ri;         % Richardson number
        Col = cmocean('ice');
        n = 'Ri';
    elseif v == 9
        Ztmp = g.hydrology.(names{i}).bbp;             % optical backscatter
        Col = cmocean('curl');
        n = 'b_b_p [m^-^1]'; 
    elseif v == 10
        Ztmp = g.hydrology.(names{i}).BI;              % backscatter index
        Col = cmocean('deep');
        n = 'BI [dB]';
    elseif v == 11
        Ztmp = g.hydrology.(names{i}).SSC;             % SPM concentration
        Col = cmocean('curl');
        n = 'SSC [mg L^-^1]';
    elseif v == 12
        Ztmp = g.hydrology.(names{i}).chla;            % chlorophyll-a
        Col = cmocean('algae');
        n = 'Chla [\mug.L^-^1]';
    end
    
    % time condition
    tArray = nanmean(Xbistmp);
    if tArray(1) > tArray(2)
        Xtmp = flip(Xtmp,2);
        Xbistmp = flip(Xbistmp,2);
        Ytmp = flip(Ytmp,2);
        Ztmp = flip(Ztmp,2);
        lattmp = flip(lattmp,2);
    end
    
    % Store variables 
    eval(['X.section' num2str(i) '=Xtmp;']);
    eval(['lat.section' num2str(i) '=lattmp;']);
    eval(['Xbis.section' num2str(i) '=Xbistmp;']);
    eval(['Y.section' num2str(i) '=Ytmp;']); 
    eval(['Z.section' num2str(i) '=Ztmp;']);
end

end