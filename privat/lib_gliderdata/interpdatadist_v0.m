function [glidmat] = interpdatadist_v0(proj,i_curstruct_vf,i_sections,diststep,dpthstep,alt_lim,deploy)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Create a regular spatial grid for glider data
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

s = i_sections;
c = i_curstruct_vf;
names = fieldnames(proj);

wbar = waitbar(0,'Spatial Interpolation, please wait...');

for i = 1:length(names)
    waitbar(i/(length(names)))
    
    %======================================================================
    % Define variables
    %======================================================================
    % Projected data
    lat = proj.(names{i}).latint;
    lon = proj.(names{i}).lonint;
    dist_proj = proj.(names{i}).distproj;
    
    % Navigation data
    time = s.(names{i}).nav(:,1);
    dpth_nav = s.(names{i}).nav(:,end);      
    
    % Science data
    temp = s.(names{i}).sci(:,2);
    cond = s.(names{i}).sci(:,4);
    chla = s.(names{i}).sci(:,5);
    bb = s.(names{i}).sci(:,6);
    if deploy{1} == 'matugli_2016' | deploy{1} == 'matugli_2017'
        cdom = s.(names{i}).sci(:,7);
    end
    
    % Depth average current data
    vx = s.(names{i}).Dac(:,end-1);
    vy = s.(names{i}).Dac(:,end);
    
    % Derived parameters
    P = s.(names{i}).param(:,2);
    SP = s.(names{i}).param(:,3);
    dens = s.(names{i}).param(:,4);
    TP = s.(names{i}).param(:,5);
    SA = s.(names{i}).param(:,6);
    CT = s.(names{i}).param(:,7);
    sigma0 = s.(names{i}).param(:,8);  
    bbp = s.(names{i}).param(:,end-1);
    ssc = s.(names{i}).param(:,end);    

    % ADCP current data
    dpth_ac = c.(names{i}).dpth;
    u = c.(names{i}).u;
    v = c.(names{i}).v;
    ustd = c.(names{i}).ustd;
    vstd = c.(names{i}).vstd;
    time_ac = c.(names{i}).time;
    
    % ADCP turbidity data
    BI = c.(names{i}).BI;
    BIstd = c.(names{i}).BIstd;
    
    %======================================================================
    % Extrapolate BI data at the surface
    %======================================================================
    BItmp = [];
    num = size(BI,1);
   
    for m = 1:size(BI,2)
        idlim = find(~isnan(BI(:,m)),1,'first');
        if idlim > 1
            tmp_BI = BI(:,m);
            idend = find(~isnan(tmp_BI),1,'last');
            tmp = inpaint_nans(tmp_BI(1:idend));
            % Keeps the size of the vector from stacking process
            tmp(end+1:num) = NaN;
            % Output
            BItmp = [BItmp,tmp];
            clearvars tmp_BI
        else
            tmp_BI = BI(:,m);
            BItmp = [BItmp,tmp_BI];
            clearvars tmp_BI
        end
    end
    BI = []; BI = BItmp;
    
    %======================================================================
    % Keep only unique data
    %======================================================================
    % Reference distance
    [distvec, indexvec] = unique(dist_proj,'stable');             % for vector
    
    % Remove duplicates for vector
    lat = lat(indexvec);
    lon = lon(indexvec);
    time = time(indexvec);
    dpth_nav = dpth_nav(indexvec);
    temp = temp(indexvec);
    cond = cond(indexvec);
    chla = chla(indexvec);
    bb = bb(indexvec);
    vx = vx(indexvec);
    vy = vy(indexvec);
    P = P(indexvec);
    SP = SP(indexvec);
    dens = dens(indexvec);
    TP = TP(indexvec);
    SA = SA(indexvec);
    CT = CT(indexvec);
    sigma0 = sigma0(indexvec);
    bbp = bbp(indexvec);
    ssc = ssc(indexvec);
    if deploy{1} == 'matugli_2016' | deploy{1} == 'matugli_2017'
        cdom = cdom(indexvec);
    end
    
    %======================================================================
    % Remove NaN values
    %======================================================================
    ind_bad = find(isnan(lat) | isnan(distvec) | isnan(dpth_nav));

    lat(ind_bad)= [];
    lon(ind_bad) = [];
    distvec(ind_bad) = [];
    time(ind_bad) = [];
    dpth_nav(ind_bad) = [];
    temp(ind_bad) = [];
    cond(ind_bad) = [];
    chla(ind_bad) = [];
    bb(ind_bad) = [];
    vx(ind_bad) = [];
    vy(ind_bad) = [];
    P(ind_bad) = []; 
    SP(ind_bad) = [];
    dens(ind_bad) = [];
    TP(ind_bad) = [];
    SA(ind_bad) = [];
    CT(ind_bad) = [];
    sigma0(ind_bad) = [];
    bbp(ind_bad) = [];
    ssc(ind_bad) = [];
    if deploy{1} == 'matugli_2016' | deploy{1} == 'matugli_2017'
        cdom(ind_bad) = [];
    end
    
    %======================================================================
    % Interp navigation distance to adcp distance
    %======================================================================
    distac = interp1(time,distvec,time_ac);
    
    % Remove NaN values
    idbad = find(isnan(distac));
    
    dpth_ac(:,idbad) = [];
    u(:,idbad) = [];
    v(:,idbad) = [];
    ustd(:,idbad) = [];
    vstd(:,idbad) = [];
    BI(:,idbad) = []; 
    BIstd(:,idbad) = []; 
    time_ac(idbad) = [];
    distac(idbad) = [];
    
    % =====================================================================
    % Interpolation on a new mesh grid by rows using linear interpolation 
    % for sparse data
    % =====================================================================
    dist = distvec;
    xref = round(min(min(dist))):diststep:round(max(max(dist)));
    yref = 1:dpthstep:200;
    [X,Y] = meshgrid(xref,yref);
        
    % Projected data
    i_lat = griddata(distvec,dpth_nav,lat,X,Y);
    i_lon = griddata(distvec,dpth_nav,lon,X,Y);
    
    % Nav data
    i_time = griddata(distvec,dpth_nav,time,X,Y);
    
    % Science data
    i_temp = griddata(distvec,dpth_nav,temp,X,Y);
    i_cond = griddata(distvec,dpth_nav,cond,X,Y);
    i_chla = griddata(distvec,dpth_nav,chla,X,Y);
    i_bb = griddata(distvec,dpth_nav,bb,X,Y);
    i_vx = griddata(distvec,dpth_nav,vx,X,Y);
    i_vy = griddata(distvec,dpth_nav,vy,X,Y);
    if deploy{1} == 'matugli_2016' | deploy{1} == 'matugli_2017'
        i_cdom = griddata(distvec,dpth_nav,cdom,X,Y);
    end
    
    % Derived param data
    i_P = griddata(distvec,dpth_nav,P,X,Y);    
    i_SP = griddata(distvec,dpth_nav,SP,X,Y);
    i_dens = griddata(distvec,dpth_nav,dens,X,Y);
    i_TP = griddata(distvec,dpth_nav,TP,X,Y);
    i_SA = griddata(distvec,dpth_nav,SA,X,Y);
    i_CT = griddata(distvec,dpth_nav,CT,X,Y);
    i_sigma0 = griddata(distvec,dpth_nav,sigma0,X,Y);
    i_bbp = griddata(distvec,dpth_nav,bbp,X,Y);
    i_ssc = griddata(distvec,dpth_nav,ssc,X,Y);    
    
    i_bbp = smooth2a(i_bbp,5,1);
    i_ssc = smooth2a(i_ssc,5,1);
    
    % Interpolate current and turbidity data
    % Griddata doesn't work (noisy result), so we do a linear
    % interpolation, first on the vertical, secondly on the horizontal
    yref = yref';
    if X(1,1) > X(1,end)
        xref = [max(X(1,:)) : -diststep : min(X(1,:))];
        n = 1;
    else
        xref = [min(X(1,:)) : diststep : max(X(1,:))];
        n = 0;
    end
    
    % Reshape distac in matrix
    distac = distac.*ones(size(u));
    
    if n==1
       distac = flip(distac,2);
       u = flip(u,2);
       v = flip(v,2);
       ustd = flip(ustd,2);
       vstd = flip(vstd,2);
       BI = flip(BI,2);
       BIstd = flip(BIstd,2);
       dpth_ac = flip(dpth_ac,2);
       xref = flip(xref);
    end
    
    % Vertical interpolation
    i_u = [];
    i_v = [];
    i_ustd = [];
    i_vstd = [];
    i_BI = [];
    i_BIstd = [];
    
    for ii = 1:size(u,2)
        
        iok = find(~isnan(u(:,ii)));
        id = find(~isnan(dpth_ac(:,ii)));
        if ~ any(iok) | length(iok) ==1
            tmpu = NaN(size(yref));
            tmpv = NaN(size(yref));
            tmpustd = NaN(size(yref));
            tmpvstd = NaN(size(yref));
        elseif length(iok)<= length(id)
            tmpu = interp1(dpth_ac(iok,ii),u(iok,ii),yref);
            tmpv = interp1(dpth_ac(iok,ii),v(iok,ii),yref);
            tmpustd = interp1(dpth_ac(iok,ii),ustd(iok,ii),yref);
            tmpvstd = interp1(dpth_ac(iok,ii),vstd(iok,ii),yref);
        elseif length(iok)>= length(id)
            tmpu = interp1(dpth_ac(id,ii),u(id,ii),yref);
            tmpv = interp1(dpth_ac(id,ii),v(id,ii),yref);
            tmpustd = interp1(dpth_ac(id,ii),ustd(id,ii),yref);
            tmpvstd = interp1(dpth_ac(id,ii),vstd(id,ii),yref);            
        end
        i_u = [i_u,tmpu];
        i_v = [i_v,tmpv];
        i_ustd = [i_ustd,tmpustd];
        i_vstd = [i_vstd,tmpvstd];
        
        clearvars iok
        
        iok = find(~isnan(BI(:,ii)));
        if ~any(iok) | length(iok) == 1
            tmpBI = NaN*(ones(size(yref)));
            tmpBIstd = NaN*(ones(size(yref)));
        elseif length(iok)<= length(id) & any(iok)
            tmpBI = interp1(dpth_ac(iok,ii),BI(iok,ii),yref);
            tmpBIstd = interp1(dpth_ac(iok,ii),BIstd(iok,ii),yref);
        elseif length(iok)>= length(id)
            tmpBI = interp1(dpth_ac(id,ii),BI(id,ii),yref);
            tmpBIstd = interp1(dpth_ac(id,ii),BIstd(id,ii),yref);
        end
        i_BI = [i_BI,tmpBI];
        i_BIstd = [i_BIstd,tmpBIstd];
        
        clearvars iok
    end
    
    % Horizontal interpolation
    uii = [];
    vii = [];
    ustdii = [];
    vstdii = [];
    BIii = [];
    BIstdii = [];
    
    for ii = 1:size(i_u,1)
        
        iok = find(~isnan(i_u(ii,:)));
        if ~ any(iok) | length(iok) ==1
            tmpu = NaN(size(xref));
            tmpv = NaN(size(xref));
            tmpustd = NaN(size(xref));
            tmpvstd = NaN(size(xref));
        else
            tmpu = interp1(distac(ii,iok),i_u(ii,iok),xref);
            tmpv = interp1(distac(ii,iok),i_v(ii,iok),xref);
            tmpustd = interp1(distac(ii,iok),i_ustd(ii,iok),xref);
            tmpvstd = interp1(distac(ii,iok),i_vstd(ii,iok),xref);
        end
        uii = [uii;tmpu];
        vii = [vii;tmpv];
        ustdii = [ustdii;tmpustd];
        vstdii = [vstdii;tmpvstd];
        
        clearvars iok
        
        iok = find(~isnan(i_BI(ii,:)));
        if ~ any(iok) | length(iok) ==1
            tmpBI = NaN(size(xref));
            tmpBIstd = NaN(size(xref));
        else
            tmpBI = interp1(distac(ii,iok),i_BI(ii,iok),xref);
            tmpBIstd = interp1(distac(ii,iok),i_BIstd(ii,iok),xref);
        end
        BIii = [BIii;tmpBI];
        BIstdii = [BIstdii;tmpBIstd];
        
        clearvars iok
    end
    
    %======================================================================
    % Derivation of the bathymetrie at each section
    %======================================================================
    % Define variables
    cast =  s.(names{i}).cast;
    dpth_dup = s.(names{i}).nav(:,end);               % dupplicate variable
    dist_dup = proj.(names{i}).distproj;              % dupplicate variable
    
    % Remove glider surfcasts
    ind_bad = find(cast == 3);
    cast(ind_bad) = [];
    dpth_dup(ind_bad) = [];
    dist_dup(ind_bad) = [];
    
    % Find bottom depth at each surfcast
    diff_profile = diff(cast);
    start_up = find(diff_profile==1);                       % find the index corresponding to the beginning of upcast profile
    ind_depth_max = start_up;
    dpth_max = dpth_dup(ind_depth_max);                     % maximum depth of the glider
    bottom_depth = dpth_max + alt_lim;                      % bottom depth
    
    % Interpolate data bottom depth
    distbot = dist_dup(ind_depth_max);                            % distance of bottom depth
    ind_bad = find(isnan(distbot) | isnan(bottom_depth));
    distbot(ind_bad) = [];
    bottom_depth(ind_bad) = [];
    [distbot, index] = unique(distbot,'stable');
    i_botdpth = interp1(distbot',bottom_depth(index)',X(1,:));
    
    %======================================================================
    % Allocate Outputs
    %======================================================================
    glidmat.(names{i}).X = X;
    glidmat.(names{i}).Y = Y;
    glidmat.(names{i}).time = i_time;
    glidmat.(names{i}).lat = i_lat;
    glidmat.(names{i}).lon = i_lon;
    glidmat.(names{i}).temp = i_temp;
    glidmat.(names{i}).cond = i_cond;
    glidmat.(names{i}).chla = i_chla;
    glidmat.(names{i}).bb = i_bb;
    glidmat.(names{i}).vx = i_vx;
    glidmat.(names{i}).vy = i_vy;
    glidmat.(names{i}).P = i_P;    
    glidmat.(names{i}).SP = i_SP;
    glidmat.(names{i}).dens = i_dens;
    glidmat.(names{i}).TP = i_TP;
    glidmat.(names{i}).SA = i_SA;
    glidmat.(names{i}).CT = i_CT;
    glidmat.(names{i}).sigma0 = i_sigma0;
    glidmat.(names{i}).bbp = i_bbp;
    glidmat.(names{i}).SSC = i_ssc;    
    glidmat.(names{i}).BI = BIii;
    glidmat.(names{i}).BIstd = BIstdii;
    glidmat.(names{i}).u = uii;
    glidmat.(names{i}).v = vii;
    glidmat.(names{i}).ustd = ustdii;
    glidmat.(names{i}).vstd = vstdii;
    if deploy{1} == 'matugli_2016' | deploy{1} == 'matugli_2017'
        glidmat.(names{i}).cdom = i_cdom;
    end
    glidmat.(names{i}).botdpth = i_botdpth;
    
    %======================================================================
    % CLEAN-UP DATA ACCORDING TO THE BATHYMETRY
    %======================================================================
    % Distance limit
    if deploy{1} == 'matugli_2016' | deploy{1} == 'matugli_2017'
        Xlim = 35;
    else
        Xlim = 60;
    end
    idbad = find(X(1,:)>Xlim | isnan(i_botdpth));
    st = glidmat.(names{i});
    namesvar = fieldnames(st);
    for ii = 1:length(namesvar)
        tmps = st.(namesvar{ii});
        tmps(:,idbad) = [];
        S.(namesvar{ii}) = tmps;
    end
    clearvars namesvar tmps st idbad ii
    
    % Replace by NaN value below the maximum depth
    namesvar = fieldnames(S);
    y = S.Y;
    nlength = 1:size(y,1);
    nlength = nlength';
    botdpth = round(S.botdpth);
    
    idbad = [];
    for ii = 1:size(y,2)
        idtmp = find(y(:,ii)>botdpth(ii));
        ni = max(numel(nlength),numel(idtmp));
        idtmp(end+1:ni) = NaN;
        idbad(:,ii) = idtmp;
        clearvars idtmp
    end
    clearvars ii

    for ii = 1:length(namesvar)-1
        tmps = S.(namesvar{ii});
        for ni = 1:size(tmps,2)
            id = idbad(:,ni);
            find(isnan(id)); id(ans) = [];
            tmps(id,ni) = NaN;
        end
        S2.(namesvar{ii}) = tmps;
        clearvars tmps;
    end
      
    % Outputs final
    S2.botdpth =  botdpth;
    glidmat.(names{i}) = S2;
        
    clearvars -except s c names i proj i_curstruct_vf i_sections ...
        diststep dpthstep alt_lim deploy dpath i_botdpth wbar waitbar glidmat
    
end
close(wbar);

end