function [proj_glid] = distglider_v0(glidstruct,startlon,startlat,endlon,endlat,proj_angle,plot_out,dpath)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Estimate distance along a projection on glider sections
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Number of glider sections
names = fieldnames(glidstruct);

for ii = 1:length(names)
    % Select data
    nav = glidstruct.(names{ii}).nav;
    DAC = glidstruct.(names{ii}).Dac;
    ref_time = nav(:,1);
    dpth = nav(:,end);
    
    % Projected distance
    [proj,h] = dist_glider_v0(DAC,ref_time,dpth,startlon,startlat,endlon,endlat,proj_angle,plot_out,dpath);
    
    % Store projected distance in glider section
    eval(['proj_glid.' names{ii} '=proj;']);
end

end