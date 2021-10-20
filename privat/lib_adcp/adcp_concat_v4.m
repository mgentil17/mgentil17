function [r_adcp_concat] = adcp_concat_v4(adressadcp, WN_ref, BinSize_ref, sal_val)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% Select and concatenate adcp data in a structure
%
% version 4.0
% 10/2019
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%


%==========================================================================
%% LOAD ADCP DATA
%==========================================================================
root = adressadcp;                                              % adcp file path
cd(root);
file_adcp = dir('*.PD0');

for i=1:length(file_adcp)
  fnames{i} = file_adcp(i).name;
end

% %----------------------------------------------------------------------
% % test au 21/07/20
% if deploy{1} == 'matugli_2018'
% %     fnames = fnames(1:87);                    % (ADCP C1)
% %     fnames = fnames(88:124);                    % (ADCP C2)
% %     fnames = fnames(125:277);                   % stormy conditions (ADCP C3)    
%     fnames = fnames(278:424);                    % (ADCP C4)
% 
% end
% %----------------------------------------------------------------------

% INITIALIZATION
% Time series
series.time = []; % Julian days
series.svel = []; % m/s
series.dpth = []; % m
series.head = []; % deg
series.ptch = []; % deg
series.roll = []; % deg
series.sal  = []; % PSU
series.temp = []; % Celcius deg
series.pres = []; % dbar
% Bottom track
% No tilt corrected
series.BTDepthB1 = [];        % Bottom depth beam1 (m)
series.BTDepthB2 = [];        % Bottom depth beam2 (m)
series.BTDepthB3 = [];        % Bottom depth beam3 (m)
series.BTDepthB4 = [];        % Bottom depth beam4 (m)
% Here in Earth coordinates, tilts corrected.
series.BTvel1 = [];        % Bottom track velocities (m/s)
series.BTvel2 = [];        % Bottom track velocities (m/s)
series.BTvel3 = [];        % Bottom track velocities (m/s)
series.BTvel4 = [];        % Bottom track velocities (m/s)

% Profiles
% Velocities
% Here in Earth coordinates, tilts corrected.
profiles.vel1 = [];
profiles.vel2 = [];
profiles.vel3 = [];
profiles.vel4 = [];
% Received level
profiles.EAcntB1 = [];   	% Received level beam1 (counts)
profiles.EAcntB2 = [];	% Received level beam2 (counts)
profiles.EAcntB3 = [];        % Received level beam3 (counts)
profiles.EAcntB4 = [];        % Received level beam4 (counts)
% Correlation
profiles.corrB1 = [];        % Correlation beam1
profiles.corrB2 = [];        % Correlation beam2
profiles.corrB3 = [];        % Correlation beam3
profiles.corrB4 = [];        % Correlation beam4

wb = waitbar(0,'Concatenating pd0 to create section...');

for i=1:length(fnames)
    waitbar(i/(length(fnames)))
    [c, s, p] = load_adcp(root, fnames{i});
    cfg=c;
    
    % If more than one configuration of the adcp is used during the deployment
    if c.BinSize > BinSize_ref || c.BinSize < BinSize_ref
        warning('Multiple ADCP configurations');
        disp('Aligning data to the reference configuration defined in Glider_ADCP_define_param.m');
      
        % Calcule des profondeurs
        dpth = c.BinSize:BinSize_ref:WN_ref;
        dpth_ref = 2:BinSize_ref:WN_ref;
    
        % Interpolation des données de profiles
        for ii = 1:size(p.vel1,1)
            pi.vel1(ii,:) = interp1(dpth,p.vel1(ii,:),dpth_ref);
            pi.vel2(ii,:) = interp1(dpth,p.vel2(ii,:),dpth_ref);
            pi.vel3(ii,:) = interp1(dpth,p.vel3(ii,:),dpth_ref);
            pi.vel4(ii,:) = interp1(dpth,p.vel4(ii,:),dpth_ref);
            pi.EAcntB1(ii,:) = interp1(dpth,p.EAcntB1(ii,:),dpth_ref);
            pi.EAcntB2(ii,:) = interp1(dpth,p.EAcntB2(ii,:),dpth_ref);
            pi.EAcntB3(ii,:) = interp1(dpth,p.EAcntB3(ii,:),dpth_ref);
            pi.EAcntB4(ii,:) = interp1(dpth,p.EAcntB4(ii,:),dpth_ref);
            pi.corrB1(ii,:) = interp1(dpth,p.corrB1(ii,:),dpth_ref);
            pi.corrB2(ii,:) = interp1(dpth,p.corrB2(ii,:),dpth_ref);
            pi.corrB3(ii,:) = interp1(dpth,p.corrB3(ii,:),dpth_ref);
            pi.corrB4(ii,:) = interp1(dpth,p.corrB4(ii,:),dpth_ref);
        end
        
        % Allocate in variables
        p.vel1 = pi.vel1;
        p.vel2 = pi.vel2;
        p.vel3 = pi.vel3;
        p.vel4 = pi.vel4;
        p.EAcntB1 = pi.EAcntB1;
        p.EAcntB2 = pi.EAcntB2;
        p.EAcntB3 = pi.EAcntB3;
        p.EAcntB4 = pi.EAcntB4;
        p.corrB1 = pi.corrB1;
        p.corrB2 = pi.corrB2;
        p.corrB3 = pi.corrB3;
        p.corrB4 = pi.corrB4;   
        
        % Configuration
        cfg.WN = WN_ref;
        cfg.BinSize = BinSize_ref;
        
        clear pi
    end
    
    % Conditions on the size of profiles (different ADCP settings along a
    % delpoyment) i.e. number of cells
    if size(p.vel1,2)<size(profiles.vel1,2)
        p.vel1(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        p.vel2(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        p.vel3(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        p.vel4(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        
        p.EAcntB1(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        p.EAcntB2(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        p.EAcntB3(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        p.EAcntB4(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        
        p.corrB1(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        p.corrB2(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        p.corrB3(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        p.corrB4(:,size(p.vel1,2):size(profiles.vel1,2))=NaN;
        
    elseif size(profiles.vel1,2)>0 && size(p.vel1,2)>size(profiles.vel1,2)
        profiles.vel1(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        profiles.vel2(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        profiles.vel3(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        profiles.vel4(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        
        profiles.EAcntB1(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        profiles.EAcntB2(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        profiles.EAcntB3(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        profiles.EAcntB4(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        
        profiles.corrB1(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        profiles.corrB2(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        profiles.corrB3(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        profiles.corrB4(:,size(profiles.vel1,2):size(p.vel1,2))=NaN;
        
    end
    
    % Time series
    series.time = [series.time; s.time]; % Julian days
    series.svel = [series.svel; s.svel]; % m/s
    series.dpth = [series.dpth; s.dpth]; % m
    series.head = [series.head; s.head]; % deg
    series.ptch = [series.ptch; s.ptch]; % deg
    series.roll = [series.roll; s.roll]; % deg
    series.sal  = [series.sal; s.sal]; % PSU
    series.temp = [series.temp; s.temp]; % Celcius deg
    series.pres = [series.pres; s.pres]; % dbar
    % Bottom track
    % No tilt corrected
    series.BTDepthB1 = [series.BTDepthB1; s.BTDepthB1];        % Bottom depth beam1 (m)
    series.BTDepthB2 = [series.BTDepthB2; s.BTDepthB2];        % Bottom depth beam2 (m)
    series.BTDepthB3 = [series.BTDepthB3; s.BTDepthB3];        % Bottom depth beam3 (m)
    series.BTDepthB4 = [series.BTDepthB4; s.BTDepthB4];        % Bottom depth beam4 (m)
    % Tilts corrected ?
    series.BTvel1 = [series.BTvel1; s.BTvel1];        % Bottom track velocities (m/s)
    series.BTvel2 = [series.BTvel2; s.BTvel2];        % Bottom track velocities (m/s)
    series.BTvel3 = [series.BTvel3; s.BTvel3];        % Bottom track velocities (m/s)
    series.BTvel4 = [series.BTvel4; s.BTvel4];        % Bottom track velocities (m/s)
    
    % Velocities
    % Tilts corrected ?
    profiles.vel1 = [profiles.vel1; p.vel1];
    profiles.vel2 = [profiles.vel2; p.vel2];
    profiles.vel3 = [profiles.vel3; p.vel3];
    profiles.vel4 = [profiles.vel4; p.vel4];
    % Received level
    profiles.EAcntB1 = [profiles.EAcntB1; p.EAcntB1];        % Received level beam1 (counts)
    profiles.EAcntB2 = [profiles.EAcntB2; p.EAcntB2];        % Received level beam2 (counts)
    profiles.EAcntB3 = [profiles.EAcntB3; p.EAcntB3];        % Received level beam3 (counts)
    profiles.EAcntB4 = [profiles.EAcntB4; p.EAcntB4];        % Received level beam4 (counts)
    % Correlation
    profiles.corrB1 = [profiles.corrB1; p.corrB1];        % Correlation beam1
    profiles.corrB2 = [profiles.corrB2; p.corrB2];        % Correlation beam2
    profiles.corrB3 = [profiles.corrB3; p.corrB3];        % Correlation beam3
    profiles.corrB4 = [profiles.corrB4; p.corrB4];        % Correlation beam4
    
    fclose all;
end
close(wb)

% fileout = ['SECTION_ns.mat'];
% plot(series.time, -series.dpth, 'k.'); grid on; datetick('x')
% save(fileout, 'cfg', 'series', 'profiles')

% Filter bad values of temperature (21°C)
idbad = find(series.temp == 21);
series.temp(idbad) = NaN;

% Sets the salinity at 38
series.sal = sal_val.*ones(size(series.sal));
%-------------------------------------------------

%==========================================================================
%% ALLOCATE OUTPUT
%==========================================================================
r_adcp_concat = struct('cfg',cfg,'series',series,'profiles',profiles);


end