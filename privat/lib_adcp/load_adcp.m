function [cfg, series, profiles] = load_adcp(dir, fname)

[fid,message] = fopen([dir '/' fname],'r','l');

[fd,vd,veld,cmd,ead,pgd,status,btd] = rdread(fid);

% Configuration
cfg.f	    = fd(1);	% Sampling frequency (kHz)
cfg.up      = fd(2);	% Orientation (upward if 1)
cfg.theta   = fd(3);	% Beam angle (Deg)
cfg.WN      = fd(4);	% Number of depth cells
cfg.WP      = fd(5);	% Pings per ensemble
cfg.BinSize = fd(6); 	% Bin size (m)
cfg.WF      = fd(7); 	% Blank after transmit
switch fd(8)
  case 0
    cfg.coord = 'beam';
  case 1
    cfg.coord = 'instrument';
  case 2
    cfg.coord = 'ship';
  case 3
    cfg.coord = 'earth';
end
cfg.tilts   = fd(9);
cfg.S3beams = fd(10);
cfg.BMapp   = fd(11);
cfg.Bin1Mid = fd(12);	% Middle of first bin
clear fd

% overall time and depth
% Otime = squeeze(vd(:, :, 1));	% Julian days
% Odpth = squeeze(vd(:, :, 3));	% m
% plot(Otime, Odpth, 'k.')

% Time series
series.time = squeeze(vd(:, :, 1));	% Julian days
series.svel = squeeze(vd(:, :, 2));	% m/s
series.dpth = squeeze(vd(:, :, 3));	% m
series.head = squeeze(vd(:, :, 4));	% deg	
series.ptch = squeeze(vd(:, :, 5));	% deg
series.roll = squeeze(vd(:, :, 6));	% deg
series.sal  = squeeze(vd(:, :, 7));	% PSU
series.temp = squeeze(vd(:, :, 8));	% Celcius deg
series.pres = squeeze(vd(:, :, 9));	% dbar
  
% Velocities
% Tilts corrected ?
profiles.vel1 = squeeze(veld(:, :, 1));
profiles.vel2 = squeeze(veld(:, :, 2));
profiles.vel3 = squeeze(veld(:, :, 3));
profiles.vel4 = squeeze(veld(:, :, 4));
  
% Received level
profiles.EAcntB1 = squeeze(ead(:, :, 1));	% Received level beam1 (counts)
profiles.EAcntB2 = squeeze(ead(:, :, 2));	% Received level beam2 (counts)
profiles.EAcntB3 = squeeze(ead(:, :, 3));	% Received level beam3 (counts)
profiles.EAcntB4 = squeeze(ead(:, :, 4));	% Received level beam4 (counts)
  
% Correlation
profiles.corrB1 = squeeze(cmd(:, :, 1));	% Correlation beam1
profiles.corrB2 = squeeze(cmd(:, :, 2));	% Correlation beam2
profiles.corrB3 = squeeze(cmd(:, :, 3));	% Correlation beam3
profiles.corrB4 = squeeze(cmd(:, :, 4));	% Correlation beam4
  
% Bottom track
% Not tilt corrected
series.BTDepthB1 = squeeze(btd(:, :, 1));	% Bottom depth beam1 (m)
series.BTDepthB2 = squeeze(btd(:, :, 2));	% Bottom depth beam2 (m)
series.BTDepthB3 = squeeze(btd(:, :, 3));	% Bottom depth beam3 (m)
series.BTDepthB4 = squeeze(btd(:, :, 4));	% Bottom depth beam4 (m)
  
% Here in Earth coordinates, tilts corrected.
series.BTvel1 = squeeze(btd(:, :, 5));	% Bottom track velocities (m/s)
series.BTvel2 = squeeze(btd(:, :, 6));	% Bottom track velocities (m/s)
series.BTvel3 = squeeze(btd(:, :, 7));	% Bottom track velocities (m/s)
series.BTvel4 = squeeze(btd(:, :, 8));	% Bottom track velocities (m/s)


% plot(series.time, series.dpth, 'k.')
% grid on

end
