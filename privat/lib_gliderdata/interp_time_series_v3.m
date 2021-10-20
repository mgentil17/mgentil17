function [i_glidstruct,i_curstruct,i_turbstruct] = interp_time_series_v3(nav,sci,opt,curstruct,turbstruct,time_sync,nb_smooth)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Synchronization of all data on a reference time
% and
% Interpolation of data on a reference time step
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Variables
Dac = [nav.data(:,1:3) nav.data(:,8:9)];                          % depth average current
navx = [nav.data(:,1) nav.data(:,6:7) nav.data(:,10:end-2)];        % navigation data
scix = sci.data;
optx = opt.data;
cs = curstruct.ser;
cp = curstruct.pro;
ts = turbstruct.ser;
tp = turbstruct.pro;

% Kill NaN of Dac time series
idbad = find(isnan(Dac(:,end-1)) | isnan(Dac(:,end)));
Dac(idbad,:) = [];

% Transpose acoustic matrix profiles
cp.vel1 = cp.vel1';
cp.vel2 = cp.vel2';
cp.vel3 = cp.vel3';
cp.vel4 = cp.vel4';
cp.Dpth = cp.Dpth';
cp.time = cp.time';

tp.BI1 = tp.BI1';
tp.BI2 = tp.BI2';
tp.BI3 = tp.BI3';
tp.BI4 = tp.BI4';
tp.Dpth1 = tp.Dpth1';
tp.Dpth2 = tp.Dpth2';
tp.Dpth3 = tp.Dpth3';
tp.Dpth4 = tp.Dpth4';

%==========================================================================
% Synchronization and Interpolation of glider data
%==========================================================================
% Define the variables numbers of each matrices
[~,col_DAC] = size(Dac);
[~,col_nav] = size(navx);
[~,col_sci] = size(scix);
[~,col_opt] = size(optx);

% Depth average current
for i = 2:col_DAC
    Dac_sm(:,i) = medfilt1(Dac(:,i),nb_smooth);                  % smoothing of variables using median moving windows
    iDac(:,i) = interp1(Dac(:,1),Dac_sm(:,i),time_sync);         % temporal interpolation
end

% Navigation
for i = 2:col_nav                                                          
    nav_sm(:,i) = medfilt1(navx(:,i),nb_smooth);                  
    inav(:,i) = interp1(navx(:,1),nav_sm(:,i),time_sync);                
end

% Physical science
for i = 2:col_sci                                                          
    sci_sm(:,i) = medfilt1(scix(:,i),nb_smooth);                  
    isci(:,i) = interp1(scix(:,1),sci_sm(:,i),time_sync);                
end

% Bio-iptical science
for i = 2:col_opt                                                          
    opt_sm(:,i) = medfilt1(optx(:,i),nb_smooth);                  
    iopt(:,i) = interp1(optx(:,1),opt_sm(:,i),time_sync);                
end

% Add the reference time vector at each matrix
iDac = [time_sync' iDac(:,2:end)];
inav = [time_sync' inav(:,2:end)];
isci = [time_sync' isci(:,2:end)];
iopt = [time_sync' iopt(:,2:end)];

%==========================================================================
% Synchronization and Interpolation of ADCP data
% Use linear interpolation for sparse data
%==========================================================================
% Define variables 
x = cs.time;
xi = ts.time;
X = time_sync';

% Find non NaN values of current time-series
iok = find(~isnan(x));

% Time-series interpolation
cur.ser.BTvel1 = interp1(x(iok),cs.BTvel1(iok),X);
cur.ser.BTvel2 = interp1(x(iok),cs.BTvel2(iok),X);
cur.ser.BTvel3 = interp1(x(iok),cs.BTvel3(iok),X);
cur.ser.BTvel4 = interp1(x(iok),cs.BTvel4(iok),X);
cur.ser.alt = interp1(x(iok),cs.alt(iok),X);
cur.ser.sal = interp1(x(iok),cs.sal(iok),X);
cur.ser.temp = interp1(x(iok),cs.temp(iok),X);
cur.ser.dpth = interp1(x(iok),cs.dpth(iok),X);
cur.ser.head = interp1(x(iok),cs.head(iok),X);
cur.ser.ptch = interp1(x(iok),cs.ptch(iok),X);
cur.ser.roll = interp1(x(iok),cs.roll(iok),X);
cur.ser.time = X;

% Find non NaN values of turbidity time-series
iok = find(~isnan(xi));

turb.ser.sal = interp1(xi(iok),ts.sal(iok),X);
turb.ser.temp = interp1(xi(iok),ts.temp(iok),X);
turb.ser.dpth = interp1(xi(iok),ts.dpth(iok),X);
turb.ser.head = interp1(xi(iok),ts.head(iok),X);
turb.ser.ptch = interp1(xi(iok),ts.ptch(iok),X);
turb.ser.roll = interp1(xi(iok),ts.roll(iok),X);
turb.ser.altB1 = interp1(x(iok),ts.altB1(iok),X);
turb.ser.altB2 = interp1(x(iok),ts.altB2(iok),X);
turb.ser.altB3 = interp1(x(iok),ts.altB3(iok),X);
turb.ser.altB4 = interp1(x(iok),ts.altB4(iok),X);
turb.ser.time = X;

% Profiles interpolation
% Time interpolation on each of the bins
[nbins, ~] = size(cp.vel1);

for ii = 1:nbins
    cur.pro.vel1(ii,:) = interp1(x(iok),cp.vel1(ii,iok),X);
    cur.pro.vel2(ii,:) = interp1(x(iok),cp.vel2(ii,iok),X);
    cur.pro.vel3(ii,:) = interp1(x(iok),cp.vel3(ii,iok),X);
    cur.pro.vel4(ii,:) = interp1(x(iok),cp.vel4(ii,iok),X);
    cur.pro.Dpth(ii,:) = interp1(x(iok),cp.Dpth(ii,iok),X);
    cur.pro.time(ii,:) = interp1(x(iok),cp.time(ii,iok),X);        
    
    turb.pro.BI1(ii,:) = interp1(x(iok),tp.BI1(ii,iok),X);
    turb.pro.BI2(ii,:) = interp1(x(iok),tp.BI2(ii,iok),X);
    turb.pro.BI3(ii,:) = interp1(x(iok),tp.BI3(ii,iok),X);
    turb.pro.BI4(ii,:) = interp1(x(iok),tp.BI4(ii,iok),X);
    turb.pro.Dpth1(ii,:) = interp1(x(iok),tp.Dpth1(ii,iok),X);
    turb.pro.Dpth2(ii,:) = interp1(x(iok),tp.Dpth2(ii,iok),X);
    turb.pro.Dpth3(ii,:) = interp1(x(iok),tp.Dpth3(ii,iok),X);
    turb.pro.Dpth4(ii,:) = interp1(x(iok),tp.Dpth4(ii,iok),X);   
end

%==========================================================================
% Allocate Outputs
%==========================================================================
% Merge Physical and Bio-optical data
i_sci = [isci iopt(:,2:end)];

i_glidstruct.Dac = iDac;
i_glidstruct.nav = inav;
i_glidstruct.sci = i_sci;

i_curstruct.ser = cur.ser;
i_curstruct.pro = cur.pro;
i_curstruct.cfg = curstruct.cfg;

i_turbstruct.ser = turb.ser;
i_turbstruct.pro = turb.pro;
i_turbstruct.cfg = turbstruct.cfg;

end