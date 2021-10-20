function [C_fileout, T_fileout] = QA_QC_adcp_v2(r_adcp,time_lag_trick,ptch_offset,Kc,cnt_tresh,vel_tresh,ofs_BT)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Control quality of ADCP data. The velocity data still represents the sum
% of the ocean velocity observed and the glider velocity
%
% Mathieu Gentil
% 01/21
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

struct = r_adcp;

%%
if time_lag_trick
  disp('Time lag trick on depth. 5s. delay');
end

disp(['Hardware pitch offset : ' num2str(ptch_offset) ' degrees'])
disp(['Count to dB conversion factor Kc = ' num2str(Kc)])
disp(['Correlation treshold : ' num2str(cnt_tresh) ' counts'])
% disp(['Decibel limit treshold : ' num2str(db_tresh) ' decibels'])
disp(['Velocity limit treshold : ' num2str(vel_tresh) ' m.s-1'])
disp(['Bottom track offset limit : ' num2str(ofs_BT) ' m.s-1'])

%%
%==========================================================================
% INITIALIZE STRUCTURE
%==========================================================================
% Define variables
series   = struct.series;
profiles = struct.profiles;
cfg      = struct.cfg;

%..........................................................................
% Turbidity data
%..........................................................................
% Initialize pro structure
Tpro.BI1  = [];
Tpro.BI2  = [];
Tpro.BI3  = [];
Tpro.BI4  = [];
Tpro.Dpth1  = [];
Tpro.Dpth2  = [];
Tpro.Dpth3  = [];
Tpro.Dpth4  = [];

%Inititalize series structure
Tser.sal  = [];
Tser.temp  = [];
Tser.time  = [];
Tser.dpth  = [];
Tser.head  = [];
Tser.ptch  = [];
Tser.roll  = [];
Tser.altB1 = [];
Tser.altB2 = [];
Tser.altB3 = [];
Tser.altB4 = [];
Tser.bot_depth = [];

%..........................................................................
% Current data
%..........................................................................
% Initialize pro structure
Cpro.vel1 = [];
Cpro.vel2 = [];
Cpro.vel3 = [];
Cpro.vel4 = [];
Cpro.Dpth = [];  

% Initialize ser structure
Cser.BTvel1   = [];
Cser.BTvel2   = [];
Cser.BTvel3   = [];
Cser.BTvel4   = [];
Cser.alt = [];
Cser.sal    = [];
Cser.temp  = [];
Cser.time = [];
Cser.dpth = [];
Cser.head = [];
Cser.ptch = [];
Cser.roll = [];
Cser.bot_depth = [];

%%
%==========================================================================
% DATA
%==========================================================================
theta   = cfg.theta;
R       = cfg.Bin1Mid + [1:cfg.WN]*cfg.BinSize;
T       = series.temp;
S       = series.sal;
f       = cfg.f;
EA1cnt  = profiles.EAcntB1;
EA2cnt  = profiles.EAcntB2;
EA3cnt  = profiles.EAcntB3;
EA4cnt  = profiles.EAcntB4;
C1cnt   = profiles.corrB1;
C2cnt   = profiles.corrB2;
C3cnt   = profiles.corrB3;
C4cnt   = profiles.corrB4;
vel1  = profiles.vel1;
vel2  = profiles.vel2;
vel3  = profiles.vel3;
vel4  = profiles.vel4;
BT1dpth = series.BTDepthB1;
BT2dpth = series.BTDepthB2;
BT3dpth = series.BTDepthB3;
BT4dpth = series.BTDepthB4;
BTvel1  = series.BTvel1;
BTvel2  = series.BTvel2;
BTvel3  = series.BTvel3;
BTvel4  = series.BTvel4;
h_sync  = series.head;
p_sync  = series.ptch + ptch_offset;
r_sync  = series.roll;
d_sync  = series.dpth;
t_exp   = series.time;

%%
%==========================================================================
% QA/QC OF ADCP DATA
%==========================================================================
%% For our deployments replace by NaN temperature = 21
idbad = find(round(T) == 21);

if any(idbad)
    T(idbad) = NaN;
end

%% Filter surface data while glider communicating
isurf=find(d_sync<1);                    %find surface data (inf to 1m)
disp(['Filter surface data (<1m) while glider communicating :', num2str(length(isurf)), ' data were removed']);
if ~isempty(isurf)
    d_sync(isurf)=NaN;
    t_exp(isurf)=NaN;
    p_sync(isurf)=NaN;
    r_sync(isurf)=NaN;
    h_sync(isurf)=NaN;
    T(isurf)=NaN;
    S(isurf)=NaN;
    
    BT1dpth(isurf)=NaN;
    BT2dpth(isurf)=NaN;
    BT3dpth(isurf)=NaN;
    BT4dpth(isurf)=NaN;
    
    BTvel1(isurf)=NaN;
    BTvel2(isurf)=NaN;
    BTvel3(isurf)=NaN;
    BTvel4(isurf)=NaN;
    
    C1cnt(isurf,:)=NaN;
    C2cnt(isurf,:)=NaN;
    C3cnt(isurf,:)=NaN;
    C4cnt(isurf,:)=NaN;
    
    EA1cnt(isurf,:)=NaN;
    EA2cnt(isurf,:)=NaN;
    EA3cnt(isurf,:)=NaN;
    EA4cnt(isurf,:)=NaN;
    
    vel1(isurf,:)=NaN;
    vel2(isurf,:)=NaN;
    vel3(isurf,:)=NaN;
    vel4(isurf,:)=NaN;
end

%% TIME LAG TRICK (P.CAUCHY)  
  if time_lag_trick
    tmpdpth = d_sync;
    tmpt = t_exp + 5/(60*60*12);	% 5s. delay (default 5/(...))
    d_sync = interp1(t_exp, tmpdpth, tmpt);
  end
    
%% CORRECT THE REAL DEPTH OF EACH CELL ONLY FOR TURBIDITY DATA
T_sync = repmat(t_exp,1 , cfg.WN);% - t_exp(1);
% Max bottom's depth
DMax  = max(d_sync + cfg.Bin1Mid + cfg.BinSize*cfg.WN);

% Correct tilt, pitch and roll effects.
% Projection on the real vertical
disp(['Correct tilt, pitch and roll effects']);
Rv = R/cosd(theta);
altB1 = BT1dpth/cosd(theta).*cosd(theta+r_sync).*cosd(p_sync);
altB2 = BT2dpth/cosd(theta).*cosd(theta-r_sync).*cosd(p_sync);
altB3 = BT3dpth/cosd(theta).*cosd(theta+p_sync).*cosd(r_sync);
altB4 = BT4dpth/cosd(theta).*cosd(theta-p_sync).*cosd(r_sync);
R1corr = repmat(Rv, length(d_sync), 1) .* repmat(cosd(theta+r_sync).*cosd(p_sync), 1, cfg.WN); %TEST AU 20180213 - earth coordinates (?)
R2corr = repmat(Rv, length(d_sync), 1) .* repmat(cosd(theta-r_sync).*cosd(p_sync), 1, cfg.WN);
R3corr = repmat(Rv, length(d_sync), 1) .* repmat(cosd(theta+p_sync).*cosd(r_sync), 1, cfg.WN);
R4corr = repmat(Rv, length(d_sync), 1) .* repmat(cosd(theta-p_sync).*cosd(r_sync), 1, cfg.WN);
D1 = repmat(d_sync, 1, cfg.WN) + R1corr;	% Real depth of each cell
D2 = repmat(d_sync, 1, cfg.WN) + R2corr;	% Real depth of each cell
D3 = repmat(d_sync, 1, cfg.WN) + R3corr;	% Real depth of each cell
D4 = repmat(d_sync, 1, cfg.WN) + R4corr;	% Real depth of each cell

% Max bottom's depth
D1Min = max(d_sync + altB1);  if isnan(D1Min); D1Min=1000;end
D2Min = max(d_sync + altB2);  if isnan(D2Min); D2Min=1000;end
D3Min = max(d_sync + altB3);  if isnan(D3Min); D3Min=1000;end
D4Min = max(d_sync + altB4);  if isnan(D4Min); D4Min=1000;end
  
%% Backscatter Index estimation
disp(['Backscatter Index estimation in dB']);
BI1 = backscatter_corr(EA1cnt, Kc, Rv, T, S, f, d_sync, 'trdi');
BI2 = backscatter_corr(EA2cnt, Kc, Rv, T, S, f, d_sync, 'trdi');
BI3 = backscatter_corr(EA3cnt, Kc, Rv, T, S, f, d_sync, 'trdi');
BI4 = backscatter_corr(EA4cnt, Kc, Rv, T, S, f, d_sync, 'trdi');

% % %----------------------------------------------------------------------
% % % raw backscatter data
% % raw_ac.ser.sal    = S;
% % raw_ac.ser.temp  = T;
% % raw_ac.ser.time = t_exp;
% % raw_ac.ser.dpth = d_sync;
% % raw_ac.ser.head = h_sync;
% % raw_ac.ser.ptch = p_sync;
% % raw_ac.ser.roll = r_sync;
% % raw_ac.ser.altB1 = altB1;
% % raw_ac.ser.altB2 = altB2;
% % raw_ac.ser.altB3 = altB3;
% % raw_ac.ser.altB4 = altB4;
% % raw_ac.ser.bot_depth = DMax;
% % 
% % raw_ac.pro.BI1  = BI1;
% % raw_ac.pro.BI2  = BI2;
% % raw_ac.pro.BI3  = BI3;
% % raw_ac.pro.BI4  = BI4;
% % raw_ac.pro.Dpth1 = D1;
% % raw_ac.pro.Dpth2 = D2;
% % raw_ac.pro.Dpth3 = D3;
% % raw_ac.pro.Dpth4 = D4;
% % 
% % raw_ac.cfg=cfg;
% % 
% % save('raw_ac.mat','raw_ac');
% % %----------------------------------------------------------------------

%% Correlation treshold
% Apply correlation criterium to BI
[mask1, scorr1] = correlation_treshold(cnt_tresh, C1cnt, R1corr, altB1);
D1mask                   = NaN*D1;
D1mask(find(D1 < D1Min)) = 1;
BI1_c           = BI1.*mask1.*D1mask;
[mask2, scorr2] = correlation_treshold(cnt_tresh, C2cnt, R2corr, altB2);
D2mask                   = NaN*D2;
D2mask(find(D2 < D2Min)) = 1;
BI2_c           = BI2.*mask2.*D2mask;
[mask3, scorr3] = correlation_treshold(cnt_tresh, C3cnt, R3corr, altB3);
D3mask                   = NaN*D3;
D3mask(find(D3 < D3Min)) = 1;
BI3_c           = BI3.*mask3.*D3mask;
[mask4, scorr4] = correlation_treshold(cnt_tresh, C4cnt, R4corr, altB4);
D4mask                   = NaN*D4;
D4mask(find(D4 < D4Min)) = 1;
BI4_c           = BI4.*mask4.*D4mask;

tmp = C1cnt(:); tmp(isnan(tmp)) = []; id = find(tmp<=64); 
disp(['Apply correlation criterium to BI on the beam 1 :', num2str(length(id)), ' data were removed on ', num2str(length(tmp)), ' data'])
tmp = C2cnt(:); tmp(isnan(tmp)) = []; id = find(tmp<=64); 
disp(['Apply correlation criterium to BI on the beam 2 :', num2str(length(id)), ' data were removed on ', num2str(length(tmp)), ' data'])
tmp = C3cnt(:); tmp(isnan(tmp)) = []; id = find(tmp<=64); 
disp(['Apply correlation criterium to BI on the beam 3 :', num2str(length(id)), ' data were removed on ', num2str(length(tmp)), ' data'])
tmp = C4cnt(:); tmp(isnan(tmp)) = []; id = find(tmp<=64); 
disp(['Apply correlation criterium to BI on the beam 4 :', num2str(length(id)), ' data were removed on ', num2str(length(tmp)), ' data'])

%% FOR VELOCITY DATA DERIVES A SINGLE DEPTH PROFILE FOR THE 4 ADCP BEAM
% The beam mapping replaced the velocity cells of the different beams 
% in space
% Along a beam
% We choose the adcp beam 3 like a reference
% Indeed, it faces forward during a downcast and gives the most realistic
% measurement
BTdpth = BT3dpth;
alt = BTdpth/cosd(theta).*cosd(theta+p_sync).*cosd(r_sync);

Bin1dist = cfg.Bin1Mid;                                                % distance from the middle of the 1st cell of the adcp [m]
BinSize = cfg.BinSize;                                                 % [m]
WN = cfg.WN;                                                           % number of bin

Dpth = d_sync' + Bin1dist;
for j = 2:WN'
    Dpth(j,:) = Dpth(j-1,:)+BinSize;
end
Dpth = Dpth';  

%% APPLY A FILTER DEPTH ON VELOCITY DATA
% Replace by NaN velocity under the bottom depth

% Compute the bottom depth
bot_depth = d_sync + alt;

% Correction of velocity
for i = 1:length(bot_depth)
    idx = find(Dpth(i,:) > bot_depth(i));
    if idx > 1
        vel1(i,idx) = NaN;
        vel2(i,idx) = NaN;
        vel3(i,idx) = NaN;
        vel4(i,idx) = NaN;
        
        BTvel1(idx) = NaN;
        BTvel2(idx) = NaN;
        BTvel3(idx) = NaN;
        BTvel4(idx) = NaN;
    end
end

%% Correlation treshold on velocity
DMin = max(bot_depth);  if isnan(DMin); DMin=1000;end
Dmask                   = NaN*Dpth;
Dmask(find(Dpth < DMin)) = 1;

[mask1, scorr1] = correlation_treshold(cnt_tresh, C1cnt, Dpth, alt);
vel1_c           = vel1.*mask1.*Dmask;
[mask2, scorr2] = correlation_treshold(cnt_tresh, C2cnt, Dpth, alt);
vel2_c           = vel2.*mask2.*Dmask;
[mask3, scorr3] = correlation_treshold(cnt_tresh, C3cnt, Dpth, alt);
vel3_c           = vel3.*mask3.*Dmask;
[mask4, scorr4] = correlation_treshold(cnt_tresh, C4cnt, Dpth, alt);
vel4_c           = vel4.*mask4.*Dmask;

%% Velocity treshold
% On series
idbad = find(BTvel1>vel_tresh | BTvel2>vel_tresh | BTvel3>vel_tresh ...
    | BTvel4>vel_tresh | BTvel1<-vel_tresh | BTvel2<-vel_tresh |...
    BTvel3<-vel_tresh | BTvel4<-vel_tresh);

if ~isempty(idbad)
    BT1dpth(idbad)=NaN;
    BT2dpth(idbad)=NaN;
    BT3dpth(idbad)=NaN;
    BT4dpth(idbad)=NaN;
    
    BTvel1(idbad)=NaN;
    BTvel2(idbad)=NaN;
    BTvel3(idbad)=NaN;
    BTvel4(idbad)=NaN;
end

% On profiles
for ii = 1:size(vel1_c,1)
    idbad = find(vel1_c(ii,:)>vel_tresh | vel2_c(ii,:)>vel_tresh | ...
        vel3_c(ii,:)>vel_tresh | vel4_c(ii,:)>vel_tresh | ...
        vel1_c(ii,:)<-vel_tresh | vel2_c(ii,:)<-vel_tresh | ...
        vel3_c(ii,:)<-vel_tresh | vel4_c(ii,:)<-vel_tresh);
    
    vel1_c(ii,idbad) = NaN;
    vel2_c(ii,idbad) = NaN;
    vel3_c(ii,idbad) = NaN;
    vel4_c(ii,idbad) = NaN;
end

% %% Bottom track offset
% % Attention fonctionne pour matugli 2018 à vérifier pour 2016/2017
% 
% ubt_raw = BTvel1;     % east-west component
% vbt_raw = BTvel2;     % north-south component
% u_raw = nanmedian(vel1_c,2);
% v_raw = nanmedian(vel2_c,2);
% 
% % % figure;
% % % subplot 211; plot(u_raw); hold on; plot(ubt_raw);
% % % legend('uadcp','ubt'); ylim([-1 1]);
% % % title('Raw data');
% % % subplot 212; plot(v_raw); hold on; plot(vbt_raw);
% % % legend('vadcp','vbt'); ylim([-1 1]);
% % % title('Raw data');
% 
% % Filter velocity data
% tmp_filt = movmin(u_raw,200,'omitnan');
% u_filt = movmax(tmp_filt,200,'omitnan');
% tmp_filt = movmin(ubt_raw,200,'omitnan');
% ubt_filt = movmax(tmp_filt,200,'omitnan');
% 
% tmp_filt = movmin(v_raw,200,'omitnan');
% v_filt = movmax(tmp_filt,200,'omitnan');
% tmp_filt = movmin(vbt_raw,200,'omitnan');
% vbt_filt = movmax(tmp_filt,200,'omitnan');
% 
% % Compute the offset
% uoffset = u_filt - ubt_filt;
% voffset = v_filt - vbt_filt;
% uoffset = nanmedian(uoffset);
% voffset = nanmedian(voffset);
% 
% % Apply the offset
% if abs(uoffset) >= ofs_BT
%     disp('### Caution ###: application of an offset on the east-west component of the Bottom track, give priority to the DAC constraint for current correction.');
%     uref = nanmedian(u_filt);
%     for ii = 1:length(ubt_raw)
%         if u_filt(ii) >= uref
%             ubt_corr(ii) = ubt_raw(ii) + uoffset;
%         elseif u_filt(ii) >= uref
%             ubt_corr(ii) = ubt_raw(ii) - uoffset;
%         elseif u_filt(ii) == NaN
%             ubt_corr(ii) = NaN;
%         end
%     end
% else
%     ubt_corr = ubt_raw;    
% end
% BTvel1_c = ubt_corr;
% 
% if abs(voffset) >= ofs_BT
%     disp('### Caution ###: application of an offset on the north-south component of the Bottom track, give priority to the DAC constraint for current correction.');
%     vref = nanmedian(v_filt);
%     for ii = 1:length(vbt_raw)
%         if v_filt(ii) >= vref
%             vbt_corr(ii) = vbt_raw(ii) + voffset;
%         elseif u_filt(ii) >= vref
%             vbt_corr(ii) = vbt_raw(ii) - voffset;
%         elseif v_filt(ii) == NaN
%             vbt_corr(ii) = NaN;
%         end
%     end
% else
%     vbt_corr = vbt_raw;
% end
% BTvel2_c = vbt_corr;

% % figure;
% % subplot 211; plot(u_raw); hold on; plot(ubt_raw); plot(ubt_corr);
% % legend('uadcp','ubt','ubtcorr'); ylim([-1 1]);
% % title('Raw data');
% % subplot 212; plot(v_raw); hold on; plot(vbt_raw); plot(vbt_corr)
% % legend('vadcp','vbt','vbtcorr'); ylim([-1 1]);
% % title('Raw data');

%%
%==========================================================================
% ALLOCATE OUTPUTS
%==========================================================================

% Turbidity
Tser.sal    = S;
Tser.temp  = T;
Tser.time = t_exp;
Tser.dpth = d_sync;
Tser.head = h_sync;
Tser.ptch = p_sync;
Tser.roll = r_sync;
Tser.altB1 = altB1;
Tser.altB2 = altB2;
Tser.altB3 = altB3;
Tser.altB4 = altB4;
Tser.bot_depth = DMax;

Tpro.BI1  = BI1_c;
Tpro.BI2  = BI2_c;
Tpro.BI3  = BI3_c;
Tpro.BI4  = BI4_c;
Tpro.Dpth1 = D1;
Tpro.Dpth2 = D2;
Tpro.Dpth3 = D3;
Tpro.Dpth4 = D4;

T_fileout.cfg=cfg;
T_fileout.pro=Tpro;
T_fileout.ser=Tser;

% Current
Cser.BTvel1   = BTvel1;
Cser.BTvel2   = BTvel2;
Cser.BTvel3   = BTvel3;
Cser.BTvel4   = BTvel4;
Cser.alt = alt;
Cser.sal    = S;
Cser.temp  = T;
Cser.time = t_exp;
Cser.dpth = d_sync;
Cser.head = h_sync;
Cser.ptch = p_sync;
Cser.roll = r_sync;
Cser.bot_depth = bot_depth;

Cpro.vel1 = vel1_c;
Cpro.vel2 = vel2_c;
Cpro.vel3 = vel3_c;
Cpro.vel4 = vel4_c;
Cpro.Dpth = Dpth;
Cpro.time = T_sync;

C_fileout.cfg=cfg;
C_fileout.pro=Cpro;
C_fileout.ser=Cser;

end