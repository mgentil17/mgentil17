function [sm_vel] = cur_shear_method_MC_v0(dataset,constrain,alt_lim,c_glider_data,SizeStack,BinStack,nwin,nMC,stddev,plot_out)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% PART I:
%
% Stack ADCP data according to the depth to reduce the uncertainty
% associated with a single ping
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

d = dataset;
names = fieldnames(d);

tic
wbar = waitbar(0,'Data Stacking, it is a long process, please wait...');

for m = 1:length(names)
    waitbar(m/(length(names)))
    
    %======================================================================
    % DEFINE VARIABLES
    %======================================================================
    % dataset
    gnav = d.(names{m}).nav;
    gDAC = d.(names{m}).Dac;
    Guwm = d.(names{m}).Guwm;
    curser = d.(names{m}).curser;
    turbser = d.(names{m}).turbser;
    curpro = d.(names{m}).curpro;
    turbpro = d.(names{m}).turbpro;
    cfg = d.(names{m}).cfg_adcp;
    
    % Gliders's navigation
    time = gnav(:,1);
    lat  = gnav(:,2);
    lon  = gnav(:,3);
    pitch = gnav(:,5);
    heading = gnav(:,6);
    dpth_ser = gnav(:,8);
    vx = gDAC(:,4);
    vy = gDAC(:,5);
    
    % Velocities components (ug,vg,wg) from glider's flight model
    ug = Guwm(:,5);
    vg = Guwm(:,6);
    wg = Guwm(:,4);
    Hg = Guwm(:,8);                         % horizontal velocity
    
    % Velocities components (ub,vb,wb) from glider's bottom track
    ub = curser(:,1);
    vb = curser(:,2);
    wb = curser(:,3);
    
    % Velocities components (um,vm,wm) from glider's ADCP profiles
    um0 = curpro.vel1';
    vm0 = curpro.vel2';
    wm0 = curpro.vel3';
    
    % Turbidity data
    BI = turbpro.BI3';
    
    % Depths of each ADCP bins along the trajectory
    dpthcur_pro = curpro.Dpth';
    dpthturb_pro = turbpro.Dpth3';    
    
    % Profiles setup from glider's configuration
    nbins = cfg.WN;
    BinSize = cfg.BinSize;
    Bin1Mid = cfg.Bin1Mid;
    
    %======================================================================
    % DOWNCASTS IDENTIFICATION
    %======================================================================
    % Identifies idexes of start and end of downcast
    % from the depth vector
    ind_value = find(~isnan(dpth_ser));
    ind = find(diff(ind_value)>1);
    nend = ind_value(ind);
    clearvars ind_value ind
    
    ind_value = find(isnan(dpth_ser));
    ind = find(diff(ind_value)>1);
    nstart = ind_value(ind)+1;
    
    if nend(1) <= nstart(1); nend = nend(2:end); end
    if nend(end) <= nstart(end); nstart = nstart(1:end-1); end
    
    % Downcasts only
    tmp_start = nstart;
    tmp_end = nend;    
    for ii = 1:length(nstart)
        if nend(ii) - nstart(ii) <= 10
            tmp_start(ii) = [];
            tmp_end(ii) = [];
        end
    end
    nstart = tmp_start;
    nend = tmp_end;
    clearvars tmp_start tmp_end
    
    %======================================================================
    % DATA STACKING
    %======================================================================
    for ii = 1:length(nstart)
      
        %------------------------------------------------------------------
        % Profile's variables
        %------------------------------------------------------------------
        p.time = time(nstart(ii):nend(ii));
        p.lat = lat(nstart(ii):nend(ii));
        p.lon = lon(nstart(ii):nend(ii));
        p.pitch = pitch(nstart(ii):nend(ii));
        p.heading = heading(nstart(ii):nend(ii));
        p.dpth_ser= dpth_ser(nstart(ii):nend(ii));
        p.vx = vx(nstart(ii):nend(ii));
        p.vy = vy(nstart(ii):nend(ii));
        p.ug = ug(nstart(ii):nend(ii));
        p.vg = vg(nstart(ii):nend(ii));
        p.wg = wg(nstart(ii):nend(ii));
        p.Hg = Hg(nstart(ii):nend(ii));
        p.ub = ub(nstart(ii):nend(ii));
        p.vb = vb(nstart(ii):nend(ii));
        p.wb = wb(nstart(ii):nend(ii));
        p.um0 = um0(nstart(ii):nend(ii),:);
        p.vm0 = vm0(nstart(ii):nend(ii),:);
        p.wm0 = wm0(nstart(ii):nend(ii),:);
        p.BI = BI(nstart(ii):nend(ii),:);
        p.dpthcur_pro = dpthcur_pro(nstart(ii):nend(ii),:);
        p.dpthturb_pro = dpthturb_pro(nstart(ii):nend(ii),:);
        p.nbins = nbins;
        p.BinSize = BinSize;
        p.Bin1Mid = Bin1Mid;
        p.ndata = length(p.time);
        min_time = min(p.time);
        
        %------------------------------------------------------------------
        % Filter bad data
        %------------------------------------------------------------------
        % Bottom depth
        altmax = p.dpth_ser(end) + alt_lim;
        
        % Turbidity 
        idbad = p.dpthturb_pro >= altmax;
        p.BI(idbad==1) = NaN;
        p.dpthturb_pro(idbad==1) = NaN;
        clearvars idbad
        
        % Currents
        idbad = p.dpthcur_pro >= altmax;
        p.um0(idbad==1) = NaN;
        p.vm0(idbad==1) = NaN;
        p.wm0(idbad==1) = NaN;
        p.dpthcur_pro(idbad==1) = NaN;
        clearvars idbad altmax
        
        %==================================================================
        % INITIATE MONTE CARLO SIMULATIONS AT EACH PROFILE
        %==================================================================
        for imc = 1:nMC

            %--------------------------------------------------------------
            % Turbidity data
            %--------------------------------------------------------------
            % Initialize structure
            BImean  = NaN(size(p.BI));
            zstackturb = NaN(size(p.BI));
            
            % Stack and moving window 'nwin' on turbidity data
            for k = nwin+1:p.ndata-nwin
                zmin  = min(min(p.dpthturb_pro(k-nwin:k+nwin,:)));
                zmax  = max(max(p.dpthturb_pro(k-nwin:k+nwin,:)));
                if any(zmin) | any(zmax)
                    edgesturb = floor(zmin):SizeStack:ceil(zmax)+SizeStack;
                    y = discretize(p.dpthturb_pro(k-nwin:k+nwin,:),edgesturb);
                    yR = reshape(y,[],1);
                    biR = reshape(p.BI(k-nwin:k+nwin,:),[],1);
                    idbad = find(isnan(yR));
                    yR(idbad) = []; biR(idbad) = [];
                    ans1 = accumarray(yR,biR,[],@nanmedian); BImean(k,1:size(ans1,1)) = ans1;
                    zstackturb(k,1:size(ans1,1)) = edgesturb(1:size(ans1,1));
                end
            end
            
            clearvars idx_nan k zmin zmax ans1 y yR idbad edgesturb biR
            
            %--------------------------------------------------------------
            % Currents data
            %--------------------------------------------------------------
            % Generate vectors of random inputs for each components
            % um,vm ~ Error Normally Distributed
            r1 = randn(size(p.um0));
            um = ( r1 .* stddev ) + p.um0;
            vm = ( r1 .* stddev ) + p.vm0;
            
            % Initialize structure
            um_mean = NaN(size(um));
            vm_mean = NaN(size(vm));
            um_std  = NaN(size(um));
            vm_std  = NaN(size(vm));
            zstackcur  = NaN(size(vm));
                        
            % Stack and moving window 'nwin' on currents data
            for k = nwin+1:p.ndata-nwin
                zmin  = min(min(p.dpthcur_pro(k-nwin:k+nwin,:)));
                zmax  = max(max(p.dpthcur_pro(k-nwin:k+nwin,:)));
                if any(zmin) | any(zmax)
                    edgescur = floor(zmin):SizeStack:ceil(zmax)+SizeStack;
                    y = discretize(p.dpthcur_pro(k-nwin:k+nwin,:),edgescur);
                    yR = reshape(y,[],1);
                    umR = reshape(um(k-nwin:k+nwin,:),[],1);
                    vmR = reshape(vm(k-nwin:k+nwin,:),[],1);
                    idbad = find(isnan(yR));
                    yR(idbad) = []; umR(idbad) = []; vmR(idbad) = [];
                    ans = accumarray(yR,umR,[],@nanmedian); um_mean(k,1:size(ans,1)) = ans;
                    ans = accumarray(yR,vmR,[],@nanmedian); vm_mean(k,1:size(ans,1)) = ans;
                    ans = accumarray(yR,umR,[],@nanstd); um_std(k,1:size(ans,1)) = ans;
                    ans = accumarray(yR,vmR,[],@nanstd); vm_std(k,1:size(ans,1)) = ans;
                    zstackcur(k,1:size(ans,1)) = edgescur(1:size(ans,1));
                end
            end
            
            % Removal of singleton on um_mean ou vm_mean (NaN)
            idx_nan = find(isnan(um_mean) | isnan(vm_mean));
            um_mean(idx_nan) = NaN;
            vm_mean(idx_nan) = NaN;
            um_std(idx_nan) = NaN;
            vm_std(idx_nan) = NaN;
            
            clearvars idx_nan k zmin zmax ans1 ans3 y yR idbad ...
                edgescur umR vmR
            
            %**************************************************************
            % Plot
            if plot_out & ii == 10 & imc == 1
                figure('Name','Selected set of consecutive raw vs. smoothed and stacked profiles');
                subplot(1,2,1);
                for i=1:p.ndata
                    hold on;
                    plot(um(i,:),-p.dpthcur_pro(i,:),'-x');
                    title('Raw E-W current');xlabel('speed (m/s)');ylabel('depth (m)')
                end
                subplot(1,2,2);
                for i=1:nwin+1:p.ndata-nwin
                    hold on;
                    plot(um_mean(i,:),-p.dpthcur_pro(i,:),'-o');
                    title('Stacked & Smoothed E-W current');xlabel('speed (m/s)');ylabel('depth (m)')
                end
                
                figure('Name','Selected set of consecutive raw vs. smoothed and stacked profiles');
                subplot(1,2,1);
                for i=1:p.ndata
                    hold on; plot(vm(i,:),-p.dpthcur_pro(i,:),'-x');
                    title('Raw N-S velocity current');xlabel('speed (m/s)');;ylabel('depth (m)')
                end
                subplot(1,2,2);
                for i=1:nwin+1:p.ndata-nwin
                    hold on; plot(vm_mean(i,:),-p.dpthcur_pro(i,:),'-o');
                    title('Stacked & Smoothed N-S current');xlabel('speed (m/s)');;ylabel('depth (m)')
                end
            end
            %**************************************************************
            
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% PART II:
%
% Shear method - The vertical derivatives of the current profiles 
% (according to each component) are stacked per layer to build a profile
% which is then vertically integrated to find the profiles of
% each component of the current (with one constant). The profiles are
% then compared and adjusted with a constrain.
%
% i) calculate single-ensemble shear by vertically differentiating ADCP % 
% velocity profiles; 
% ii) grid resulting shear estimates in depth space; 
% iii) vertically integrate shear to yield baroclinic velocity profile
% iv) integration constant is required to add barotropic velocity component
% velocity referencing using bottom track
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

            %--------------------------------------------------------------
            % Turbidy data
            %--------------------------------------------------------------
            % Define data
            bi = BImean;
            
            % Bad value
            id_bad = find(bi == -Inf | bi == Inf);
            bi(id_bad) = NaN;
            
            % Control plots of a consecutive profiles for turbidity components
            zd = zstackturb(:,1:nbins);
            
            % Stacking of profile data per user-defined layers (BinStack)
            edgesturb = min(min(zd)):BinStack:max(max(zd)+BinStack);
            y = discretize(zd,edgesturb);
            yR = reshape(y,[],1);
            biR = reshape(bi,[],1);
            
            % Removal of missing data (NaN)
            idx_nan = find(isnan(yR));
            biR(idx_nan) = [];
            yR(idx_nan) = [];
            
            BImed(:,imc) = accumarray(yR,biR,[],@nanmean);
            edgesturb = edgesturb(1:end-1);
            
            clearvars id_bad zd y yR idx_nan
            
            %--------------------------------------------------------------
            % Currents data
            %--------------------------------------------------------------
            % Derivative currents
            [dummy,g_u_prof] = gradient(um_mean);
            [dummy,g_v_prof] = gradient(vm_mean);
            [dummy,g_z_prof] = gradient(zstackcur);
            dum = g_u_prof./g_z_prof;
            dvm = g_v_prof./g_z_prof;
            
            % Bad value
            id_bad = find(dum == -Inf | dum == Inf);
            dum(id_bad) = NaN;
            id_bad = find(dvm == -Inf | dvm == Inf);
            dvm(id_bad) = NaN;
            
            %**************************************************************
            % Plot
            if plot_out & ii == 10 & imc == 1
                % Control plots of a consecutive derivative profiles for both velocity components
                zd = zstackcur(:,1:nbins);
                figure('Name','Shear method : Selected set consecutive derivatives');
                subplot(1,2,1);
                for i=1:p.ndata
                    hold on; plot(dum(i,:),-zd(i,:),'-x');xlabel('speed (m/s)');ylabel('depth (m)')
                    title('Consecutive E-W derivatives');
                end
                subplot(1,2,2);
                for i=1:p.ndata
                    hold on; plot(dvm(i,:),-zd(i,:),'-o');xlabel('speed (m/s)');ylabel('depth (m)')
                    title('Consecutive N-S derivatives');
                end
            end
            %**************************************************************

            % Stacking of profile data per user-defined layers (BinStack)
            zd = zstackcur(:,1:nbins);
            edges = min(min(zd)):BinStack:max(max(zd)+BinStack);
            y = discretize(zd,edges);
            yR = reshape(y,[],1);
            dumR = reshape(dum,[],1);
            dvmR = reshape(dvm,[],1);
            
            % Removal of missing data (NaN)
            idx_nan = find(isnan(yR));
            dumR(idx_nan) = [];
            dvmR(idx_nan) = [];
            yR(idx_nan) = [];
            
            dummed = accumarray(yR,dumR,[],@nanmedian);
            dvmmed = accumarray(yR,dvmR,[],@nanmedian);
            dumstd = accumarray(yR,dumR,[],@nanstd);
            dvmstd = accumarray(yR,dvmR,[],@nanstd);
            num = max(numel(dummed));
            
            % Determines the interpolation limit
            id_u = find(~isnan(dummed),1,'last');
            id_v = find(~isnan(dvmmed),1,'last');
            
            if id_u >= id_v
                id = id_u;
            else
                id = id_v;
            end
            
            % Derivative interpolation on the defined limit
            dummed = inpaint_nans(dummed(1:id));
            dvmmed = inpaint_nans(dvmmed(1:id));
            dumstd = inpaint_nans(dumstd(1:id));
            dvmstd = inpaint_nans(dvmstd(1:id));
            
            % Keeps the size of the vector from stacking process
            dummed(end+1:num) = NaN;
            dvmmed(end+1:num) = NaN;
            dumstd(end+1:num) = NaN;
            dvmstd(end+1:num) = NaN;
            
            % Vertical integration from the surface
            ui(:,imc) = cumtrapz(dummed)*BinStack;
            vi(:,imc) = cumtrapz(dvmmed)*BinStack;
            
            %**************************************************************
            % Plot
            if plot_out & ii == 10 & imc == 1
                figure('Name','Shear method: Stacked derivatives and reconstructed current');
                subplot(1,2,1), plot(dummed,-edges(1:end-1),'r'); grid minor;
                title('Derivative E-W current');xlabel('du/dz (1/s)');ylabel('depth (m)')
                subplot(1,2,2); plot(ui,-edges(1:end-1),'b'); grid minor;
                title('Reconstructed E-W current');xlabel('speed (m/s)');ylabel('depth (m)')
                
                figure('Name','Shear method: Stacked derivatives and reconstructed current');
                subplot(1,2,1); plot(dvmmed,-edges(1:end-1),'r'); grid minor;
                title('Derivative N-S current');xlabel('dv/dz (1/s)');ylabel('depth (m)')
                subplot(1,2,2); plot(vi,-edges(1:end-1),'b'); grid minor;
                title('Reconstructed N-S current');xlabel('speed (m/s)');ylabel('depth (m)')
            end
            %**************************************************************
            
            clearvars id_bad id_u id_v zd y yR idx_nan
 
        end
        
        % Calculate summary statistics
        ui_mean = nanmean(ui,2);
%         ui_std  = nanstd(ui,2);
        ui_std  = nanstd(ui,0,2);        
        vi_mean = nanmean(vi,2);
%         vi_std  = nanstd(vi,2);
        vi_std  = nanstd(vi,0,2);
        BI_mean = nanmean(BImed,2);
%         BI_std = nanstd(BImed,2); 
        BI_std = nanstd(BImed,0,2);         
        
        clearvars ui vi BImed BIstd
        
        %******************************************************************
        % Navigation data and Bottom track data
        %******************************************************************
        timenav = c_glider_data.data(:,1);
        vxp = c_glider_data.data(:,8);                                % clean glider data
        vyp = c_glider_data.data(:,9);                                % clean glider data
        latp = c_glider_data.data(:,2);
        lonp = c_glider_data.data(:,3);
        
        % Kill NaN
        idx_NaN = find(isnan(vxp) | isnan(vyp));
        timenav(idx_NaN) = [];
        vxp(idx_NaN) = [];
        vyp(idx_NaN) = [];
        latp(idx_NaN) = [];
        lonp(idx_NaN) = [];
        
        % Find the depth average current corresponding to each absolute
        % velocity profile
        % Use the time vector
        diff_time = abs(timenav - min_time);
        [~,idx] = min(diff_time);
        
        vxp = vxp(idx);
        vyp = vyp(idx);
        latp = latp(idx);
        lonp = lonp(idx);
        timenav= timenav(idx);
        
        % bottom track
        ubp = repmat(p.ub,1,nbins);
        vbp = repmat(p.vb,1,nbins);
        
        if constrain == 1
            %..............................................................
            % Apply the Bottom Track Constrain
            %..............................................................
            cst = 'Bottom track';
            
            % Compute speed profiles using the bottom track (referenced
            % profiles) = near bottom current
            ucst = p.um0 - ubp;
            vcst = p.vm0 - vbp;
            
            % Removal of singleton on unb ou vnb (NaN)
            idx_nan = find(isnan(ucst) | isnan(vcst));
            ucst(idx_nan) = NaN;
            vcst(idx_nan) = NaN;
            
            % Stacking of referenced profiles
            y = discretize(p.dpthcur_pro,edges);
            yR = reshape(y,[],1);
            ucstR = reshape(ucst,[],1);
            vcstR = reshape(vcst,[],1);
            
            % Removal of missing data (NaN)
            idx_nan = find(isnan(yR));
            ucstR(idx_nan) = [];
            vcstR(idx_nan) = [];
            yR(idx_nan) = [];
            
            ucstmed = accumarray(yR,ucstR,[],@nanmedian);
            vcstmed = accumarray(yR,vcstR,[],@nanmedian);
            ucststd = accumarray(yR,ucstR,[],@nanstd);
            vcststd = accumarray(yR,vcstR,[],@nanstd);
            
            % Filter bad bottom track data
            tmp = nanstd(ucstmed);
            if tmp >= 0.13
                idbad = find(abs(ucstmed)>1.5*tmp);
                ucstmed(idbad) = NaN;
                ucststd(idbad) = NaN;
            end
            tmp = nanstd(vcstmed);
            if tmp >= 0.13
                idbad = find(abs(vcstmed)>1.5*tmp);
                vcstmed(idbad) = NaN;
                vcststd(idbad) = NaN;
            end
            
            % Smooth the bottom track constrain
            ucstmed = movmedian(ucstmed,3);
            vcstmed = movmedian(vcstmed,3);
            ucststd = movmedian(ucststd,3);
            vcststd = movmedian(vcststd,3);
            
            % Offset
            uoffset = nanmedian(ui_mean-ucstmed);
            voffset = nanmedian(vi_mean-vcstmed);
            
        elseif constrain == 2
            %..............................................................
            % Apply the DAC Track Constrain
            %..............................................................
            cst = 'Depth Average Current';
            
            % Create a size vector of
            ucstmed = ones(size(ui_mean))*vxp;
            vcstmed = ones(size(ui_mean))*vyp;
            ucststd = ones(size(ui_mean))*0.05;             % Uncertainty of GPS 5 cm.s-1 (Merckelbach et al., 2008)
            vcststd = ones(size(ui_mean))*0.05;             % Uncertainty of GPS 5 cm.s-1 (Merckelbach et al., 2008)
            
            % Offset
            uoffset = nanmedian(ui_mean) - nanmedian(ucstmed);
            voffset = nanmedian(vi_mean) - nanmedian(vcstmed);
        end
        
        %******************************************************************
        % Plot
        if plot_out & ii == 10
            figure('Name','Shear Method: baroclinic and referenced profiles')
            subplot(1,2,1), plot(ui_mean,-edges(2:end),'r','LineWidth',3); hold on; plot(ucstmed,-edges(2:end),'-r');
            plot(ucstmed-ucststd,-edges(2:end),'--r');plot(ucstmed+ucststd,-edges(2:end),'--r');
            grid minor;
            title('Baroclinic and referenced E-W current ');
            xlabel('speed (m/s)');ylabel('depth (m)')
            subplot(1,2,2); plot(vi_mean,-edges(2:end),'b','LineWidth',3); hold on; plot(vcstmed,-edges(2:end),'-b');
            plot(vcstmed-vcststd,-edges(2:end),'--b');plot(vcstmed+vcststd,-edges(2:end),'--b');
            grid minor;
            title('Baroclinic and referenced N-S current');
            xlabel('speed (m/s)');ylabel('depth (m)')
        end
        %******************************************************************
        
        % Shifts the velocity profiles with respect to the bottom track
        % current profile
        ui_adjust = ui_mean - uoffset;
        vi_adjust = vi_mean - voffset;
        ui_adjust_std = nanmedian(sqrt(ui_std.^2 + ucststd.^2));
        vi_adjust_std = nanmedian(sqrt(vi_std.^2 + vcststd.^2));
        
        % Create a vector
        ui_adjust_std = ones(size(ui_adjust))*ui_adjust_std;
        vi_adjust_std = ones(size(vi_adjust))*vi_adjust_std;
    
        edges = edges(1:end-1);
        
        %******************************************************************
        % Plot
        if plot_out & ii == 10
            figure('Name','Shear Method: Adjusted profiles')
            subplot(1,2,1), plot(ui_adjust,-edges,'r','LineWidth',3); hold on; plot(ucstmed,-edges,':k','LineWidth',2);
            plot(ucstmed-ucststd,-edges,'--k');plot(ucstmed+ucststd,-edges,'--k');
            grid minor;
            xlabel('speed (m/s)');ylabel('depth (m)');axis tight
            title('Adjusted and referenced E-W current ')
            subplot(1,2,2); plot(vi_adjust,-edges,'b','LineWidth',3); hold on; plot(vcstmed,-edges,':k','LineWidth',2);
            plot(vcstmed-vcststd,-edges,'--k');plot(vcstmed+vcststd,-edges,'--k');
            grid minor;
            title('Adjusted and referenced N-S current');
            xlabel('speed (m/s)');ylabel('depth (m)');axis tight
        end
        %******************************************************************
        
        clearvars ucstR vcstR zdepth idx_nan id_bad zd y yR
        
        %==================================================================
        % INTERP DATA
        %==================================================================
        ref_dpth = edges';
        i_BI = interp1(edgesturb',BI_mean,ref_dpth);
        i_BIstd = interp1(edgesturb',BI_std,ref_dpth);
        
        %==================================================================
        % ALLOCATE OUTPUTS PROFILES
        %==================================================================
        sm{ii}.timedac = timenav;
        sm{ii}.lat = latp;
        sm{ii}.lon = lonp;
        sm{ii}.vx = vxp;
        sm{ii}.vy = vyp;
        sm{ii}.min_time = min_time;
        sm{ii}.edges = edges;
        sm{ii}.um0 = p.um0;
        sm{ii}.vm0 = p.vm0;
        sm{ii}.dummed = dummed;
        sm{ii}.dvmmed = dvmmed;
        sm{ii}.uc = ui_mean;
        sm{ii}.vc = vi_mean;
        sm{ii}.uc_std = ui_std;
        sm{ii}.vc_std = vi_std;
        sm{ii}.ub = ubp;
        sm{ii}.vb = vbp;
        sm{ii}.cst = cst;
        sm{ii}.ucst = ucstmed;
        sm{ii}.vcst = vcstmed;
        sm{ii}.ucst_std = ucststd;
        sm{ii}.vcst_std = vcststd;
        sm{ii}.uoffset = uoffset;
        sm{ii}.voffset = voffset;
        sm{ii}.uc_adjust = ui_adjust;
        sm{ii}.vc_adjust = vi_adjust;
        sm{ii}.uc_adjust_std = ui_adjust_std;
        sm{ii}.vc_adjust_std = vi_adjust_std;
        sm{ii}.BI_mean = i_BI;
        sm{ii}.BI_std = i_BIstd;

        clearvars p timenav latp lonp vxp vyp min_time ui_mean vi_mean ui_std ...
            vi_std ubp vbp ucst vcstb ucstmed vcstmed ucststd vcststd uoffset ...
            voffset ui_adjust vi_adjust ui_adjust_std vi_adjust_std ...
            i_BI i_BIstd
    end
    
    %======================================================================
    % ALLOCATE OUTPUTS SECTIONS
    %======================================================================
    nstart_prof.(names{m}) = nstart;
    nend_prof.(names{m}) = nend;
    sm_vel.(names{m}) = sm;
    
    clearvars sm nstart nend
end
t = toc;
fprintf('Monte-Carlo simulations on a section... [%2.2f sec]\n', t);
close(wbar);
                     
end