function [dm_vel] = cur_direct_method_MC_v0(dataset,constrain,alt_lim,c_glider_data,SizeStack,BinStack,nwin,nMC,stddev)

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

for i = 1:length(names)
    waitbar(i/(length(names)))
    
    %======================================================================
    % DEFINE VARIABLES
    %======================================================================
    % dataset
    gnav = d.(names{i}).nav;
    gDAC = d.(names{i}).Dac;
    Guwm = d.(names{i}).Guwm;
    curser = d.(names{i}).curser;
    turbser = d.(names{i}).turbser;
    curpro = d.(names{i}).curpro;
    turbpro = d.(names{i}).turbpro;
    cfg = d.(names{i}).cfg_adcp;
    
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
            
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% PART II:
%
% Direct method - Difference in measured ADCP speeds and static speed
% of the glider determined by the flight model
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

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
            % Velocities components (ug,vg,wg) from glider's flight model
            ugp = repmat(p.ug,1,nbins);
            vgp = repmat(p.vg,1,nbins);

            % Compute absolute velocities
            uc = um_mean - ugp;                    % east-west component
            vc = vm_mean - vgp;                    % north-south component

            % Stacking of profile data per user-defined layers (BinStack)
            edges = min(min(zstackcur)-BinStack):BinStack:max(max(zstackcur)+BinStack);
            y = discretize(zstackcur,edges);
            yR = reshape(y,[],1);
            ucR = reshape(uc,[],1);
            vcR = reshape(vc,[],1);
            
            % Removal of missing data (NaN)
            idx_nan = find(isnan(yR));
            ucR(idx_nan) = [];
            vcR(idx_nan) = [];
            yR(idx_nan) = [];
            
            umean(:,imc) = accumarray(yR,ucR,[],@nanmean);
            vmean(:,imc) = accumarray(yR,vcR,[],@nanmean);
            
            clearvars y yR ucR vcR idx_nan
        end
    
        % Calculate summary statistics
        u_mean = nanmean(umean,2);
        u_std  = nanstd(umean,2);
        v_mean = nanmean(vmean,2);
        v_std  = nanstd(vmean,2);
        BI_mean = nanmean(BImed,2);
        BI_std = nanstd(BImed,2); 
        
        clearvars umean vmean BImed 
        
% %         %******************************************************************
% %         % Navigation data and Bottom track data
% %         %******************************************************************
% %         timenav = c_glider_data.data(:,1);
% %         vxp = c_glider_data.data(:,8);                                % clean glider data
% %         vyp = c_glider_data.data(:,9);                                % clean glider data
% %         latp = c_glider_data.data(:,2);
% %         lonp = c_glider_data.data(:,3);
% %         
% %         % Kill NaN
% %         idx_NaN = find(isnan(vxp) | isnan(vyp));
% %         timenav(idx_NaN) = [];
% %         vxp(idx_NaN) = [];
% %         vyp(idx_NaN) = [];
% %         latp(idx_NaN) = [];
% %         lonp(idx_NaN) = [];
% %         
% %         % Find the depth average current corresponding to each absolute
% %         % velocity profile
% %         % Use the time vector
% %         diff_time = abs(timenav - min_time);
% %         [~,idx] = min(diff_time);
% %         
% %         vxp = vxp(idx);
% %         vyp = vyp(idx);
% %         latp = latp(idx);
% %         lonp = lonp(idx);
% %         timenav= timenav(idx);
% %         
% %         % bottom track
% %         ubp = repmat(p.ub,1,nbins);
% %         vbp = repmat(p.vb,1,nbins);
% %         
% %          if constrain == 1
% %             %**************************************************************
% %             % Apply the Bottom Track Constrain
% %             %**************************************************************
% %             cst = 'Bottom track';
% %             
% %             % Compute speed profiles using the bottom track (referenced
% %             % profiles) = near bottom current
% %             ucst = p.um0 - ubp;
% %             vcst = p.vm0 - vbp;
% %             
% %             % Removal of singleton on unb ou vnb (NaN)
% %             idx_nan = find(isnan(ucst) | isnan(vcst));
% %             ucst(idx_nan) = NaN;
% %             vcst(idx_nan) = NaN;
% %             
% %             % Stacking of referenced profiles
% %             y = discretize(p.dpthcur_pro,edges);
% %             yR = reshape(y,[],1);
% %             ucstR = reshape(ucst,[],1);
% %             vcstR = reshape(vcst,[],1);
% %             
% %             % Removal of missing data (NaN)
% %             idx_nan = find(isnan(yR));
% %             ucstR(idx_nan) = [];
% %             vcstR(idx_nan) = [];
% %             yR(idx_nan) = [];
% %             
% %             ucstmed = accumarray(yR,ucstR,[],@nanmedian);
% %             vcstmed = accumarray(yR,vcstR,[],@nanmedian);
% %             ucststd = accumarray(yR,ucstR,[],@nanstd);
% %             vcststd = accumarray(yR,vcstR,[],@nanstd);
% %             
% %             % Filter bad bottom track data
% %             tmp = nanstd(ucstmed);
% %             if tmp >= 0.13
% %                 idbad = find(abs(ucstmed)>1.5*tmp);
% %                 ucstmed(idbad) = NaN;
% %                 ucststd(idbad) = NaN;
% %             end
% %             tmp = nanstd(vcstmed);
% %             if tmp >= 0.13
% %                 idbad = find(abs(vcstmed)>1.5*tmp);
% %                 vcstmed(idbad) = NaN;
% %                 vcststd(idbad) = NaN;
% %             end
% %             
% %             % Smooth the bottom track constrain
% %             ucstmed = movmedian(ucstmed,3);
% %             vcstmed = movmedian(vcstmed,3);
% %             ucststd = movmedian(ucststd,3);
% %             vcststd = movmedian(vcststd,3);
% %             
% %             % Offset
% %             uoffset = nanmedian(u_mean-ucstmed);
% %             voffset = nanmedian(v_mean-vcstmed);
% %             
% %          elseif constrain == 2
% %             %**************************************************************
% %             % Apply the DAC Track Constrain
% %             %**************************************************************
% %             cst = 'Depth Average Current';
% %             
% %             % Create a size vector of
% %             ucstmed = ones(size(u_mean))*vxp;
% %             vcstmed = ones(size(u_mean))*vyp;
% %             ucststd = ones(size(u_mean))*0.05;             % Uncertainty of GPS 5 cm.s-1 (Merckelbach et al., 2008)
% %             vcststd = ones(size(u_mean))*0.05;             % Uncertainty of GPS 5 cm.s-1 (Merckelbach et al., 2008)
% %             
% %             % Offset
% %             uoffset = nanmedian(u_mean) - nanmedian(ucstmed);
% %             voffset = nanmedian(v_mean) - nanmedian(vcstmed);
% %          end
% %             
% %         % Shifts the velocity profiles with respect to the bottom track
% %         % current profile or Depth average current
% %         umean_adjust = u_mean - uoffset;
% %         vmean_adjust = v_mean - voffset;
% %         umean_adjust_std = nanmedian(sqrt(u_std.^2 + ucststd.^2));
% %         vmean_adjust_std = nanmedian(sqrt(v_std.^2 + vcststd.^2));
        
        umean_adjust = u_mean;
        vmean_adjust = v_mean;
        umean_adjust_std = u_std;
        vmean_adjust_std = v_std;
        
%         % Create a vector
%         umean_adjust_std = ones(size(umean_adjust))*umean_adjust_std;
%         vmean_adjust_std = ones(size(vmean_adjust))*vmean_adjust_std;
    
        edges = edges(1:end-1);
        
        clearvars ucstR vcstR zdepth idx_nan id_bad zd y yR
        
        %==================================================================
        % INTERP DATA
        %==================================================================
        ref_dpth = edges';
        i_BI= interp1(edgesturb',BI_mean,ref_dpth);
        i_BIstd = interp1(edgesturb',BI_std,ref_dpth);
            
        %==================================================================
        % ALLOCATE OUTPUTS PROFILES
        %==================================================================
%         dm{ii}.timedac = timenav;
%         dm{ii}.lat = latp;
%         dm{ii}.lon = lonp;
%         dm{ii}.vx = vxp;
%         dm{ii}.vy = vyp;
        dm{ii}.min_time = min_time;
        dm{ii}.edges = edges;
        dm{ii}.um0 = p.um0;
        dm{ii}.vm0 = p.vm0;
        dm{ii}.ug = p.ug;
        dm{ii}.vg = p.vg;
        dm{ii}.uc = u_mean;
        dm{ii}.vc = v_mean;
        dm{ii}.uc_std = u_std;
        dm{ii}.vc_std = v_std;
%         dm{ii}.ub = ubp;
%         dm{ii}.vb = vbp;
%         dm{ii}.cst = cst;
%         dm{ii}.ucst = ucstmed;
%         dm{ii}.vcst = vcstmed;
%         dm{ii}.ucst_std = ucststd;
%         dm{ii}.vcst_std = vcststd;
%         dm{ii}.uoffset = uoffset;
%         dm{ii}.voffset = voffset;
        dm{ii}.uc_adjust = umean_adjust;
        dm{ii}.vc_adjust = vmean_adjust;
        dm{ii}.uc_adjust_std = umean_adjust_std;
        dm{ii}.vc_adjust_std = vmean_adjust_std;
        dm{ii}.BI_mean = i_BI;
        dm{ii}.BI_std = i_BIstd;

        clearvars p timenav latp lonp vxp vyp min_time u_mean v_mean u_std ...
            v_std ubp vbp ucst vcstb ucstmed vcstmed ucststd vcststd uoffset ...
            voffset umean_adjust vmean_adjust umean_adjust_std vmean_adjust_std ...
            i_BI i_BIstd
    end
    
    %======================================================================
    % ALLOCATE OUTPUTS SECTIONS
    %======================================================================
    dm_vel.(names{i}) = dm;
    
    clearvars dm nstart nend
    
end
t = toc;
fprintf('Monte-Carlo simulations on a section... [%2.2f sec]\n', t);
close(wbar);

end