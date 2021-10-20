function [sm_vel_filt, DAC_cur] = smooth_vel(sm_vel,constrain,max_dpth,nwinc,nwinr,plot_out,dpath)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% Smooth velocity data 
% Using smooth2a function
% Smooths a 2D matrix using a mean filter over a user-defined rectangle. 
% Ignores and preserves NaNs.
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

names = fieldnames(sm_vel);                       % number of sections

for ii = 1:length(names)
    sm_cur = sm_vel.(names{ii});  
    num = num2str(ii);
    
    %======================================================================
    % DATA PREPARATION
    %======================================================================
    % Initialize time and position vectors
    time = NaN(1,length(sm_cur));
    lat = NaN(1,length(sm_cur));                            % one value at each profile
    lon = NaN(1,length(sm_cur));                            % one value at each profile
    dpth = NaN(1,length(sm_cur));

    % Depth average current
    % Initialize vector
    vx = NaN(1,length(sm_cur));
    vy = NaN(1,length(sm_cur));

    % Currents from shear method
    % Initialize matrix
    u_sm = NaN(max_dpth,length(sm_cur));                    % east-west component
    v_sm = NaN(max_dpth,length(sm_cur));                    % north-south component
    ustd_sm = NaN(max_dpth,length(sm_cur));                 % uncertainties                 
    vstd_sm = NaN(max_dpth,length(sm_cur));                 % uncertainties              
    z_sm = NaN(max_dpth,length(sm_cur));                    % depth
    
    ui_std = NaN(max_dpth,length(sm_cur)); 
    vi_std = NaN(max_dpth,length(sm_cur));  
    unbstd = NaN(max_dpth,length(sm_cur)); 
    vnbstd = NaN(max_dpth,length(sm_cur)); 

    % Initialize matrix of bottom track
    ub = NaN(max_dpth,length(sm_cur));
    vb = NaN(max_dpth,length(sm_cur));
    
    % Backscatter index
    % Initialize matrix
    BI_dm = NaN(max_dpth,length(sm_cur));
    BI_sm = NaN(max_dpth,length(sm_cur));
    BIstd_dm = NaN(max_dpth,length(sm_cur));
    BIstd_sm = NaN(max_dpth,length(sm_cur));
    
    %======================================================================
    % FLAG BY NAN PROFILES WITHOUT CONSTRAIN VALUE 
    %======================================================================
    for j = 1:length(sm_cur)
        if ~any(sm_cur{j}.uoffset) || ~any(sm_cur{j}.voffset)
            sm_cur{j}.uc_adjust = NaN(size(sm_cur{j}.vx));
            sm_cur{j}.vc_adjust = NaN(size(sm_cur{j}.vx));
        end
    end

    % Convert cell to matrix
    struct_sm = cell2mat(sm_cur);                       % shear currents data
    
    for j = 1:length(sm_cur)
        
        %==================================================================
        % EXTRACT VELOCITY DATA FROM BOTH METHODS
        %==================================================================
        time(1:length(struct_sm(j).min_time),j) = struct_sm(j).min_time;
        lat(1:length(struct_sm(j).lat),j) = struct_sm(j).lat;
        lon(1:length(struct_sm(j).lon),j) = struct_sm(j).lon;
        dpth(1:length(struct_sm(j).edges'),j) = struct_sm(j).edges';
        
        % Depth Average Current
        vx(1:length(struct_sm(j).vx(1)),j) = struct_sm(j).vx(1);
        vy(1:length(struct_sm(j).vy(1)),j) = struct_sm(j).vy(1);
        
        % Shear method
        u_sm(1:length(struct_sm(j).uc_adjust),j) = struct_sm(j).uc_adjust;
        v_sm(1:length(struct_sm(j).vc_adjust),j) = struct_sm(j).vc_adjust;
        ustd_sm(1:length(struct_sm(j).uc_adjust_std),j) = struct_sm(j).uc_adjust_std;
        vstd_sm(1:length(struct_sm(j).vc_adjust_std),j) = struct_sm(j).vc_adjust_std;
        z_sm(1:length(struct_sm(j).edges),j) = struct_sm(j).edges;
        
        ui_std(1:length(struct_sm(j).uc_std),j) = struct_sm(j).uc_std;
        vi_std(1:length(struct_sm(j).vc_std),j) = struct_sm(j).vc_std;
        
        % Constrain
        ucst(1:length(struct_sm(j).ucst),j) = struct_sm(j).ucst;
        vcst(1:length(struct_sm(j).vcst),j) = struct_sm(j).vcst;
        ucst_std(1:length(struct_sm(j).ucst_std),j) = struct_sm(j).ucst_std;
        vcst_std(1:length(struct_sm(j).vcst_std),j) = struct_sm(j).vcst_std;
        
        % Bottom track
        ub(1:length(struct_sm(j).ub(:,1)),j) = struct_sm(j).ub(:,1);
        vb(1:length(struct_sm(j).vb(:,1)),j) = struct_sm(j).vb(:,1);
        
        % Backscatter index
        BI_sm(1:length(struct_sm(j).BI_mean),j) = struct_sm(j).BI_mean;
        BIstd_sm(1:length(struct_sm(j).BI_std),j) = struct_sm(j).BI_std;
    end

    %======================================================================
    % SMOOTHING VELOCITY DATA
    %======================================================================
    % Shear method
    u_sm_filt = smooth2a(u_sm,nwinc,nwinr);
    v_sm_filt = smooth2a(v_sm,nwinc,nwinr);
    ustd_sm_filt = smooth2a(ustd_sm,nwinc,nwinr);
    vstd_sm_filt = smooth2a(vstd_sm,nwinc,nwinr);
    
    % Backscatter index
    BI_sm_filt = smooth2a(BI_sm,nwinc,nwinr);
    BIstd_sm_filt = smooth2a(BIstd_sm,nwinc,nwinr);
    
    if plot_out == 1
        %==================================================================
        % PLOT
        %==================================================================
        h = figure('Name','Raw and Filt Absolute velocities from Shear Method');
        subplot 221; contourf(u_sm); axis ij; cmocean('balance'); colorbar; caxis([-0.5 0.5]); ylim([0 60]); title('Raw uc (m s^-^1)');ylabel('Depth (bin of 2 meters)','Fontweight','bold');
        subplot 222; contourf(v_sm); axis ij; cmocean('balance'); colorbar; caxis([-0.5 0.5]); ylim([0 60]); title('Raw vc (m s^-^1)');
        subplot 223; contourf(u_sm_filt); axis ij; cmocean('balance'); colorbar; caxis([-0.5 0.5]); ylim([0 60]); title('Smooth uc (m s^-^1)'); xlabel('Profiles number','Fontweight','bold'); ylabel('Depth (bin of 2 meters)','Fontweight','bold');
        subplot 224; contourf(v_sm_filt); axis ij; cmocean('balance'); colorbar; caxis([-0.5 0.5]); ylim([0 60]); title('Smooth vc (m s^-^1)'); xlabel('Profiles number','Fontweight','bold');
        cd(dpath);
        saveas(h,'Smooth_shear_velocities.jpeg');
    end
    
    %======================================================================
    % ALLOCATE OUTPUTS
    %======================================================================
    if constrain == 1
        cst = 'Bottom track';
    elseif constrain == 2
        cst = 'Depth Average Current';
    end
    
    sm_vel_filt.(names{ii}).time = time;
    sm_vel_filt.(names{ii}).dpth = z_sm;
    sm_vel_filt.(names{ii}).u = u_sm_filt;
    sm_vel_filt.(names{ii}).v = v_sm_filt;
    sm_vel_filt.(names{ii}).ustd = ustd_sm_filt;
    sm_vel_filt.(names{ii}).vstd = vstd_sm_filt;
    sm_vel_filt.(names{ii}).cst = cst;
    sm_vel_filt.(names{ii}).BI = BI_sm_filt;
    sm_vel_filt.(names{ii}).BIstd = BIstd_sm_filt;
    
    DAC_cur.(names{ii}).time = time;
    DAC_cur.(names{ii}).lat = lat;
    DAC_cur.(names{ii}).lon = lon;
    DAC_cur.(names{ii}).vx = vx;
    DAC_cur.(names{ii}).vy = vy;
    
end