function [glidmat_pro] = der_spatial_param_v0(glidmat,deploy)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Derivation of spatial parameters
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

g = glidmat;
names = fieldnames(glidmat);

wbar = waitbar(0,'Derived parameters, please wait...');

for i = 1:length(names)
    waitbar(i/(length(names)))
    
    %======================================================================
    % Filter bad value
    %======================================================================
    % Remove NaN profiles
    idbad = find(isnan(nanmedian(g.(names{i}).CT)) | ...
        isnan(nanmedian(g.(names{i}).u)));
    
    if ~any(idbad)
        s = g.(names{i});
        namesvar = fieldnames(s);
        for n = 1:length(namesvar)
            tmps = s.(namesvar{n});
            tmps(:,idbad) = [];
            S.(namesvar{n}) = tmps;
        end
        clearvars namesvar tmps s
    else
        S = g.(names{i});
    end
    
    % Replace by NaN bad value of bb due to the surfacings of glider
    bb = S.bb;
    bbp = S.bbp;
    ssc = S.SSC;
    bbtmp = bb;
    bbptmp = bbp;
    ssctmp = ssc;
    X = S.X;
    Y = S.Y;
    
    for ii = 1:size(X,2)
        if deploy{1} == 'matugli_2016'
            if X(1,ii) > 15
                for jj = 1:size(X,1)
                    if Y(jj,ii)<=30                     % depth limit
                        if bbtmp(jj,ii) >= 0.001
                            bbtmp(jj,ii) = NaN;
                            bbptmp(jj,ii) = NaN;
                            ssctmp(jj,ii) = NaN;
                        end
                    end
                end
            end
        elseif deploy{1} == 'matugli_2017'
            if X(1,ii) > 22
                for jj = 1:size(X,1)
                    if Y(jj,ii)<=50                     % depth limit
                        if bbtmp(jj,ii) >= 0.0013
                            bbtmp(jj,ii) = NaN;
                            bbptmp(jj,ii) = NaN;
                            ssctmp(jj,ii) = NaN;
                        end
                    end
                end
            end
        elseif deploy{1} == 'matugli_2018'
            if X(1,ii) > 10
                for jj = 1:size(X,1)
                    if Y(jj,ii)<=30                     % depth limit
                        if bbtmp(jj,ii) >= 0.6
                            bbtmp(jj,ii) = NaN;
                            bbptmp(jj,ii) = NaN;
                            ssctmp(jj,ii) = NaN;
                        end
                    end
                end
            end
        end
    end
    
    % Extrapolate NaN values
    idlim = 50;
    bbitot = [];
    bbpitot = [];
    sscitot = [];
    for ii = 1:size(X,2)
        idbad = find(isnan(bbtmp(1:idlim,ii)));        
        if any(idbad)
            bbi = inpaint_nans(bbtmp(1:idlim,ii),4);
            bbpi = inpaint_nans(bbptmp(1:idlim,ii),4);
            ssci = inpaint_nans(ssctmp(1:idlim,ii),4);
            
            num = max(numel(bbptmp(:,ii)));
            bbi(end+1:num) = bbtmp(idlim+1:end,ii);
            bbpi(end+1:num) = bbptmp(idlim+1:end,ii);
            ssci(end+1:num) = ssctmp(idlim+1:end,ii); 
        else
            bbi = bbtmp(:,ii);
            bbpi = bbptmp(:,ii);
            ssci = ssctmp(:,ii);
        end
        bbitot = [bbitot,bbi];
        bbpitot = [bbpitot,bbpi];
        sscitot = [sscitot,ssci];
    end
    S.bb = bbitot;
    S.bbp = bbpitot;
    S.SSC = sscitot;
    
    %======================================================================
    % Compute Brunt Vaisala Frequency
    %======================================================================
    % Define variables
    SA = S.SA;
    CT = S.CT;
    p = S.P;
    lat = S.lat;
    
    % Limit of pressure (124 m)
    id = find(round(p) == 125);
    p(id) = NaN;
    
    % Brunt Vaisala frequency
    [N2, ~] = gsw_Nsquared(SA,CT,p,lat);
    
    clearvars SA CT p lat id
    
    % Adjust matrix on N2 size
    s = S;
    namesvar = fieldnames(s);
    
    for n = 1:length(namesvar)-1
        tmps = s.(namesvar{n});
        tmps = tmps(1:end-1,:);
        S.(namesvar{n}) = tmps;
    end
    clearvars namesvar tmps s n
    
    %======================================================================
    % Compute Richardson Number
    %======================================================================
    % Define variables
    u = S.u;             % along-shelf component
    v = S.v;             % cross-shef component
    dpth = S.Y;
    
    % Compute the magnitude of total current
    U = sqrt(u.^2 + v.^2);
    
    % Shear stress 
    dUm = diff2(U')./diff2(dpth');        % Centered derivative
    dUm = dUm';
    
    % Richardson Number
    Ri = N2 ./ (dUm.^2);

    % Allocate Outputs
    S.N2 = N2;
    S.Ri = Ri;
    S.U = U;
    S.dUm = dUm;
    
    clearvars u v dpth U dUm Ri
    
    %======================================================================
    % Compute Relative Geostrophic current
    %======================================================================
    % Delta distance
    dist = S.X;
    delt_dist = diff(dist, 1, 2);                   % difference between rows
    delt_dist = delt_dist.*1000;                    % [m]
    
    % Adjust matrix on delta distance size
    s = S;
    namesvar = fieldnames(s);
    for n = 1:length(namesvar)
        tmps = s.(namesvar{n});
        tmps = tmps(:,1:end-1);
        S.(namesvar{n}) = tmps;
    end  
    botdpth = S.botdpth;
    S = rmfield(S,namesvar(end-4));
    clearvars namesvar tmps s n
    
    % Define variables
    dist = S.X;
    dist_vg = S.X + ((delt_dist/2)./1000);   % [km]
    dpth = S.Y;                 % same that preessure in coastal area
    lat = S.lat;
    SA = S.SA;
    CT = S.CT;
    u = S.u;
    
    % Coriolis parameter
    omega = 7.27e-5;
    fcoriol = 2 * omega * sin(pi/180 * (lat));

    % Smooth physical parameters
    SA_sm = NaN(size(SA));
    CT_sm = NaN(size(CT));
    for ii = 1:size(SA,1)
        SA_sm(ii,:) = mean_av(SA(ii,:),5);            % every 1500 m
        CT_sm(ii,:) = mean_av(CT(ii,:),5);        
    end
    
    % Define the reference depth to compute geopotential anomalies at each
    % profiles
    % We choose as a reference 5m below the surface, where each profile
    % contains a value
    ref_p = 5;
    
    % Geopotential anomalies
    for n = 1:size(dist,2)
        delt_phi(:,n) = gsw_geo_strf_dyn_height(SA_sm(:,n),CT_sm(:,n),dpth(:,n),ref_p);   
    end
    delt_phi = smooth2a(delt_phi,5,5);
    
    % Geostrophic current
    [row,col] = size(delt_phi);
    vg = NaN(row,col-1);
    for ii = 2:col
        for jj = 1:row
            vg(jj,ii) = -(1/(delt_dist(jj,ii)*fcoriol(jj,ii)).*([delt_phi(jj,ii)-delt_phi(jj,ii-1)])); % m.s-1
        end
    end   

    % Barotropic adjustment to retrieve absolute geostrophic velocity
    vgmean = nanmean(vg);
    vgmed = nanmedian(vg);
    umean = nanmean(u);
    umed = nanmedian(u);
    
    vgxmean = umean-vgmean;     
    vgxmed = umed-vgmed; 

    vgxcorr_mean = NaN(size(vg));
    vgxcorr_med = NaN(size(vg));    
    for n = 1:size(vg,1)
        vgxcorr_mean(n,:) = vg(n,:)+vgxmean;
        vgxcorr_med(n,:) = vg(n,:)+vgxmed;        
    end
    S.Vg = vgxcorr_med;
    Ri = S.Ri;

    %=======================================================================
    % SMOOTH DATA CURRENT
    %=======================================================================
    S.u = smooth2a(S.u,5,2);
    S.v = smooth2a(S.v,5,2);
    S.U = smooth2a(S.U,5,2);
    S.Vg = smooth2a(S.Vg,5,2);
    
   %=======================================================================
   % ALLOCATE OUTPUTS 
   %=======================================================================
   % Position data
   glidmat_pro.position.(names{i}).X = S.X;
   glidmat_pro.position.(names{i}).Y = S.Y;
   glidmat_pro.position.(names{i}).time = S.time;   
   glidmat_pro.position.(names{i}).lat = S.lat;   
   glidmat_pro.position.(names{i}).lon = S.lon;   
   glidmat_pro.position.(names{i}).distVg = dist_vg;   
   glidmat_pro.position.(names{i}).refpVg = ref_p;
   glidmat_pro.position.(names{i}).botdpth = botdpth;
  
   % Hydrological data
   glidmat_pro.hydrology.(names{i}).temp = S.temp;
   glidmat_pro.hydrology.(names{i}).cond = S.cond;
   glidmat_pro.hydrology.(names{i}).chla = S.chla;
   glidmat_pro.hydrology.(names{i}).bb = S.bb;
   glidmat_pro.hydrology.(names{i}).bbp = S.bbp;
   glidmat_pro.hydrology.(names{i}).SSC = S.SSC;
   glidmat_pro.hydrology.(names{i}).p = S.P;
   glidmat_pro.hydrology.(names{i}).SP = S.SP;
   glidmat_pro.hydrology.(names{i}).dens = S.dens;
   glidmat_pro.hydrology.(names{i}).TP = S.TP;
   glidmat_pro.hydrology.(names{i}).SA = S.SA;
   glidmat_pro.hydrology.(names{i}).CT = S.CT;
   glidmat_pro.hydrology.(names{i}).sigma0 = S.sigma0;
   glidmat_pro.hydrology.(names{i}).BI = S.BI;
   glidmat_pro.hydrology.(names{i}).BIstd = S.BIstd;   
   glidmat_pro.hydrology.(names{i}).N2 = S.N2;   
   if deploy{1} == 'matugli_2016' | deploy{1} == 'matugli_2017'
      glidmat_pro.hydrology.(names{i}).cdom = S.cdom;
   end
    
   % Hydrodynamical data
   glidmat_pro.hydrodynamique.(names{i}).U = S.U;
   glidmat_pro.hydrodynamique.(names{i}).u = S.u;
   glidmat_pro.hydrodynamique.(names{i}).v = S.v;
   glidmat_pro.hydrodynamique.(names{i}).ustd = S.ustd;
   glidmat_pro.hydrodynamique.(names{i}).vstd = S.vstd;
   glidmat_pro.hydrodynamique.(names{i}).vx = S.vx;
   glidmat_pro.hydrodynamique.(names{i}).vy = S.vy;
   glidmat_pro.hydrodynamique.(names{i}).Vg = S.Vg;
   glidmat_pro.hydrodynamique.(names{i}).dUm = S.dUm;
   glidmat_pro.hydrodynamique.(names{i}).Ri = Ri;
   
   clearvars -except names g i glidmat deploy dpath glidmat_pro wbar
end
close(wbar);

end