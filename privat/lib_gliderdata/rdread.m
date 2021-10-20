function varargout = rdread(fid)
%RDREAD Read RDI BB data.
%  Adapted from RDREAD(FID)

%  Christian Mertens, IfM Kiel

% rewind to beginning of file and read header to get the number of bytes
% in each ensemble, the number of data types, and the address offsets
status = fseek(fid,0,'bof');
[nbytes,dtype,offset] = rdhead(fid);
% disp([' raw data has ',int2str(nbytes),'+2 bytes per ensemble'])
ntypes = length(dtype);

% get the number of ensembles from file size; each ensemble has nbytes
% plus two bytes for the checksum
status = fseek(fid,0,'eof');
m = floor(ftell(fid)/(nbytes + 2));
status = fseek(fid,0,'bof');

% number of bins is the offset difference (minus 2 bytes for the ID code)
% between velocity data and correlation magnitude devided by 4 beams
% times 2 bytes
n = (offset(4) - offset(3) - 2)/(2*4);

% data parameters
scale = [NaN,NaN,0.001,1,1,1,1,0.001];
precision = {'','','int16','uint8','uint8','uint8','uint8',''};
%varid = [0 128 256 512 768 1024 1536]; % 0     : fixed leader
                                        % 128   : Variable leader
                                        % 256   : Velocity
                                        % 512   : Correlation magnitude
                                        % 768   : Echo Intensity
                                        % 1024  : Percent good
                                        % 1280  : Status data
                                        % 1536  : Bottom Track
varid = [0,128,256,512,768,1024,1280,1536];
bad = scale.*[NaN,NaN,-32768,0,NaN,NaN,1,-32768];

% initialize output variables
for k = 1:length(varid)
  if k == 1
    varargout{k} = NaN*ones(1,1,12);
  elseif k == 2
    varargout{k} = NaN*ones(m,1,9);
  elseif (k >= 3 & k <= 7)
    varargout{k} = NaN*ones(m,n,4);
  elseif k == 8
    varargout{k} = NaN*ones(m,1,16);
  end
end

% read fixed leader data
status = fseek(fid,offset(1)+2,'bof');
varargout{1} = rdflead(fid);

icheck=0;

for i = 1:m
  % read ensemble to verify the checksum
  status = fseek(fid,(i-1)*(nbytes+2),'bof');
  buffer = fread(fid,nbytes,'uint8');
  checksum = fread(fid,1,'uint16');

  % read ensemble if checksum is ok
  if checksum == rem(sum(buffer),2^16);
    for kk = 2:length(dtype)
      k = dtype(kk);
      % set file pointer to beginning of data
      status = fseek(fid,(i-1)*(nbytes+2)+offset(kk)+2,'bof');
      switch varid(k)
        case varid(2)
          % variable leader data
          varargout{k}(i,1,:) = rdvlead(fid);
        case varid(8)
          % bottom track data
          varargout{k}(i,1,:) = rdbtrack(fid);
        otherwise
          % velocity, correlation, echo intensity, percent-good or status data
          a = fread(fid,4*n,precision{k});
          varargout{k}(i,:,:) = scale(k)*reshape(a,4,n)';
      end
    end
  else
   icheck=icheck+1;
  end
end

if icheck > m*0.01
 disp([' WARNING  found ',int2str(icheck),' ensembles with bad checksum '])
end

% check for bad values
i = find(varargout{3} == bad(3));
varargout{3}(i) = NaN;
i = find(varargout{4} == bad(4));
varargout{4}(i) = NaN;
% bottom track
i = find(varargout{8}(:,1,5) == bad(8));
varargout{8}(i,1,1:8) = NaN;
i = find(varargout{8}(:,1,6) == bad(8));
varargout{8}(i,1,1:8) = NaN;
i = find(varargout{8}(:,1,7) == bad(8));
varargout{8}(i,1,1:8) = NaN;
i = find(varargout{8}(:,1,8) == bad(8));
varargout{8}(i,1,1:8) = NaN;



%-------------------------------------------------------------------------------

function [nbytes,dtype,offset] = rdhead(fid)
%RDHEAD Read the header data from a raw ADCP data file.
%  [NBYTES,DTYPE,OFFSET] = RDHEAD(FID)

hid = 127;  % header identification byte
sid = 127;  % data source identification byte

% get file position pointer
fpos = ftell(fid);

% check header and data source identification bytes
[id,n] = fread(fid,2,'uint8');
if  (n < 2 | feof(fid))
  error('Unexpected end of file.')
end
if (id(1) ~= hid | id(2) ~= sid)
  error('Header identification byte not found.')
end

% read the number of bytes
nbytes = fread(fid,1,'uint16');

% skip spare byte
fseek(fid,1,'cof');

% read the number of data types
ndt = fread(fid,1,'uint8');
 if ndt >= 9; ndt=8; end;        %%% DT bug fix 2009-01-07

% read address offsets for data types
offset = fread(fid,ndt,'uint16');

% read variable identifiers
%varid = [0 128 256 512 768 1024 1536]; % 0 	: fixed leader
					% 128	: Variable leader
					% 256 	: Velocity
					% 512 	: Correlation magnitude
					% 768 	: Echo Intensity
					% 1024 	: Percent good
					% 1280	: Status data
					% 1536	: Bottom Track
varid = [0 128 256 512 768 1024 1280 1536];
for i = 1:ndt
  fseek(fid,fpos+offset(i),'bof');
  id = fread(fid,1,'uint16');
  dtype(i) = find(id == varid);
end

% rewind to the beginning of the ensemble
fseek(fid,fpos,'bof');


%-------------------------------------------------------------------------------

function fl = rdflead(fid);
%RDFLEAD Read the fixed leader data from a raw ADCP data file.
%  FL = RDFLEAD(FID)

fseek(fid,2,'cof');
% System configuration
cfg = fread(fid,2,'uchar');
% Interpretation.
R = mod(cfg(1), 8);	% Frequency
if R==0
  fl(1) = 75;	% (kHz)
elseif R==1
  fl(1) = 150;	% (kHz)
elseif R==2
  fl(1) = 300;	% (kHz)
elseif R==3
  fl(1) = 600;	% (kHz)
elseif R==4
  fl(1) = 1200;	% (kHz)
elseif R==5
  fl(1) = 2400;	% (kHz)
else
  disp('Error decoding frequency')
end
fl(2) = fix(cfg(1)/ 2^8);	% Uplooking if 1, Downlooking otherwise
Th = mod(cfg(2), 4);		% Beam angle
if Th==0
  fl(3) = 15;
elseif Th==1
  fl(3) = 20;
elseif Th==2
  fl(3) = 30;
else
  disp('Error decoding beam angle')
end

fseek(fid,3,'cof');

% number of depth cells
fl(4) = fread(fid,1,'uint8');

% pings per ensemble, depth cell length in cm, blank after transmit
fl(5:7) = fread(fid,3,'uint16');
fl(6) = 0.01*fl(6);
fl(7) = 0.01*fl(7);

fseek(fid,9,'cof');
c = fread(fid,1,'uint8');
% Interpretation
R = mod(c, 32);
if R<8
  fl(8) = 0;	% Beam coordinates
elseif R<16
  fl(8) = 1;	% Instrument coordinates
elseif R<24
  fl(8) = 2;	% Ship coordinates
else
  fl(8) = 3;	% Earth coordinates
end
R = mod(c, 8);
if R<4
  fl(9) = 0;	% Tilts not used in coordinates transformation
else
  fl(9) = 1;	% Tilts used in coordinates transformation
end
R = mod(c, 4);
if R<2
  fl(10) = 0; 	% 3-Beam solution not used
else
  fl(10) = 1; 	% 3-Beam solution used
end
R = mod(c, 2);
if R==1
  fl(11) = 1;	% Bin mapping used
else
  fl(11) = 0; 	% Bin mapping not used
end

fseek(fid,6,'cof');

% Bin 1 distance
fl(12) = 0.01*fread(fid,1,'ushort');



%-------------------------------------------------------------------------------

function vl = rdvlead(fid)
%PICLOD modified 20160325
%RDVLEAD Read the variable leader data from a raw ADCP data file.
%  VL = RDVLEAD(FID)

fseek(fid,2,'cof');

% time of ensemble
c = fread(fid,7,'uint8');
c(1)=y2k(c(1));
vl(1) = datenum(c(1),c(2),c(3),c(4),c(5),c(6)+c(7)/100);
fseek(fid,3,'cof');

% speed of sound (EC)
vl(2) = fread(fid,1,'uint16');

% Depth of transducer (ED)
vl(3) = 0.1*fread(fid,1,'uint16');

% heading (EH)
vl(4) = 0.01*fread(fid,1,'uint16');

% pitch (EP) and roll (ER)
vl(5:6) = 0.01*fread(fid,2,'int16');

% salinity (ES)
vl(7) = fread(fid,1,'uint16');

% temperature (ET)
vl(8) = 0.01*fread(fid,1,'int16');
fseek(fid,20,'cof');

% Pressure (P) in dbar
vl(9) = 0.001*fread(fid,1,'uint32');




%-------------------------------------------------------------------------------

function bt = rdbtrack(fid)
%RDBTRACK Read the bottom track data from a raw ADCP data file.
%  BT = RDBTRACK(FID)

%fseek(fid,14,'cof');
fseek(fid,7,'cof');
BTMode = fread(fid,1,'uint8');
fseek(fid,6,'cof');


% range
bt(1:4) = 0.01*fread(fid,4,'uint16');

% velocity
bt(5:8) = 0.001*fread(fid,4,'int16');

% correlation and magnitude
bt(9:16) = fread(fid,8,'uint8');

%-------------------------------------------------------------------------------

function d=y2k(d)
% fix date string
if d<80, d=2000+d; end
if d<100, d=1900+d; end


%-------------------------------------------------------------------------------

