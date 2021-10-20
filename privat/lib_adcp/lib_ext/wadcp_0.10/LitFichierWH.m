function [Variables]=LitFichierWH(fname,Nd,Nf)
%function [Variables]=LitFichierWH(fname,Nd,Nf)
% Fonction qui lit d'un coup les pings d'un fichier binaire WorkHorse
% en espérant (ce qui est d'ordinaire vÃ©rifiÃ©) que tous ont le mÃªme format.
% Louis Marié, le 09/03/06

fid=fopen(fname,'rb','ieee-le');
Tout=fread(fid,10,'*uint8');
fclose(fid);

if ((Tout(1)~=hex2dec('7f')) | (Tout(2)~=hex2dec('7f')))
    disp('Le fichier ne commence pas par un Binary Header valide');
    return;
end;


% On calcule le nombre d'octets par pings.
NBytes=double(ChampUInt16(Tout(3:4)))+2;

% On lit les pings demandés.
fid=fopen(fname,'rb','ieee-le');
fseek(fid,Nd*NBytes,'bof');
Tout=fread(fid,(Nf-Nd+1)*NBytes,'*uint8');
fclose(fid);

if (length(Tout)~=NBytes*floor(length(Tout)/NBytes))
    Tout=Tout(1:NBytes*floor(length(Tout)/NBytes));
end
Tout=reshape(Tout(1:end),NBytes,length(Tout)/NBytes);
ibad=find((Tout(1,:)~=hex2dec('7f')) | (Tout(2,:)~=hex2dec('7f')));
Tout(:,ibad)=[];

Variables=ParsePingWH(Tout);
return;



function [Variables]=ParsePingWH(Tout)
%function [Variables]=ParsePingWH(Tout)
% Fonction qui interprète toutes les colonnes de la matrice Tout comme des
% Pings sauvés par WH, en espérant (ce qui est d'ordinaire vérifié)
% que tous ont le même format. Louis Marié, le 09/03/06

% On s'offre le luxe de vérifier les checksums
chksum=mod(sum(Tout(1:end-2,:),1,'double'),65536);
if (chksum~=(double(ChampUInt16(Tout(end-1:end,:)))))
    fprintf('Le Checksum est faux sur l''un des pings.\n');
    return;
end;

% On lit les Binary Headers, pour avoir les offsets de début des variables.
Variables(1)=ParseVariableWH(Tout);
LesOffsets=sort(Variables(1).VariableData.DTOffsets(:,1));
LesOffsets(end+1)=size(Tout,1)-2;

% On lit les Fixed Leaders, pour avoir les nombres de cellules par pings.
Variables(2)=ParseVariableWH(Tout(LesOffsets(1)+1:LesOffsets(2),:));

for u=2:length(LesOffsets)-1
    Variables(u+1)=ParseVariableWH(Tout(LesOffsets(u)+1:LesOffsets(u+1),:));
end;

return;



function Variable=ParseVariableWH(Tout)
% function Variable=ParseVariableWH(Tout)
% Grosse fonction qui lit un tableau d'uint8,
% et qui interprÃ¨te ses colonnes comme des Ã©chantillons d'une
% variable sauvÃ©es par WH pour une sÃ©rie de pings.


Variable.VId=unique(ChampUInt16(Tout(1:2,:)));

if (1~=length(Variable.VId)) fprintf('Les variables ne sont pas dans le même ordre pour tous les pings!');return;end;

switch(Variable.VId)
    case hex2dec('7f7f')  % Header
              VD.NBytes  =ChampUInt16(Tout(3:4,:));
              NDT        =Tout(6,:);
              if (1~=length(unique(NDT))) fprintf('Il n''y a pas le même nombre de variables dans tous les pings!\n');return;end;
              for u=1:NDT(1);
                  VD.DTOffsets(u,:)=ChampUInt16(Tout((u-1)*2+7:(u-1)*2+8,:));
              end;
              
    case hex2dec('0000')  % Fixed Leader
              VD.FWVer       =double(Tout(3,:))+double(Tout(4,:))/100;            
              VD.SysConf     =ChampUInt16(Tout(5:6,:));
             [VD.Frequency,...
              VD.Convex,...
              VD.SnsConf,...
              VD.XdcrAtt,...
              VD.Up,...
              VD.BeamAngle,...
              VD.BConf]=decode_sysconfig(Tout(5,:),Tout(6,:));
              VD.RealSim     =Tout(7,:);
              VD.NBeams      =double(Tout(9,:));
              VD.NCells      =double(Tout(10,:));
              VD.NPings      =double(ChampUInt16(Tout(11:12,:)));
              VD.CellDepth   =double(ChampUInt16(Tout(13:14,:)))/100;
              VD.BlankDepth  =double(ChampUInt16(Tout(15:16,:)))/100;

              VD.SPMode      =Tout(17,:);
              VD.CorrThresh  =Tout(18,:);
              VD.NCodeReps   =Tout(19,:);
              VD.PGMin       =Tout(20,:);
              VD.EVMax       =double(ChampUInt16(Tout(21:22,:)))/1000;
              
              TppMin         =double(Tout(23,:));
              TppSec         =double(Tout(24,:));
              TppCent        =double(Tout(25,:));
              VD.Tpp         =60*TppMin+TppSec+TppCent/100;
              
              VD.CoordXform  =Tout(26,:);
              VD.HeadingAlgn =double(ChampInt16(Tout(27:28,:)))/100;
              VD.HeadingBias =double(ChampInt16(Tout(29:30,:)))/100;
              VD.SensorSrc   =Tout(31,:);
              VD.SensorAv    =Tout(32,:);
              VD.Bin1Dstnc   =double(ChampUInt16(Tout(33:34,:)))/100;
              VD.XmtLength   =double(ChampUInt16(Tout(35:36,:)))/100;
              VD.RefLayStart =Tout(37,:);
              VD.RefLayStop  =Tout(38,:);
              VD.FalseTgt    =Tout(39,:);
              VD.XmtLagDst   =double(ChampUInt16(Tout(41:42,:)))/100;
              if length(Tout)<43    %changed for BB ADCP
                    VD.CPUSN       =NaN(7,1);
                    VD.BandWidth   =NaN(2,1);
                    VD.SysPow      =NaN(1,1);
%               else
%                     VD.CPUSN       =Tout(43:50,:);
%                     VD.BandWidth   =ChampUInt16(Tout(51:52,:));
%                     VD.SysPow      =Tout(53,:);
              end
         
    case hex2dec('0080')  % Variable Leader
              ENum              =ChampUInt16(Tout(3:4,:));
              
              VD.Year           =double(Tout(58,:))*100+double(Tout(59,:));
              VD.Month          =double(Tout(6,:));
              VD.Day            =double(Tout(7,:));
              VD.Hour           =double(Tout(8,:));
              VD.Min            =double(Tout(9,:));
              
              Sec               =double(Tout(10,:));
              Cent              =double(Tout(11,:));
              VD.Sec            =Sec+Cent/100;
              
              ENumMSB           =Tout(12,:);
              VD.ENum           =65536*double(ENumMSB)+double(ENum);
              
              VD.BIT            =ChampUInt16(Tout(13:14,:));
              
              VD.SoundSpeed     =double(ChampUInt16(Tout(15:16,:)));
              VD.TransducerDepth=double(ChampUInt16(Tout(17:18,:)))/10;
              VD.Heading        =double(ChampUInt16(Tout(19:20,:)))/100;
              VD.Pitch          =double(ChampInt16(Tout(21:22,:)))/100;
              VD.Roll           =double(ChampInt16(Tout(23:24,:)))/100;
              VD.Salinity       =double(ChampUInt16(Tout(25:26,:)));
              VD.Temp           =double(ChampInt16(Tout(27:28,:)))/100;
              
              VD.MPTMin         =Tout(29,:);
              VD.MPTSec         =double(Tout(30,:))+double(Tout(31,:))/100;

              VD.HdgStd         =double(Tout(32,:));
              VD.PtchStd        =double(Tout(33,:)/10);
              VD.RollStd        =double(Tout(34,:)/10);
              
              VD.ADC0           =Tout(35,:);
              VD.ADC1           =Tout(36,:);
              VD.ADC2           =Tout(37,:);
              VD.ADC3           =Tout(38,:);
              VD.ADC4           =Tout(39,:);
              VD.ADC5           =Tout(40,:);
              VD.ADC6           =Tout(41,:);
              VD.ADC7           =Tout(42,:);
              
              VD.ESW            =ChampUInt32(Tout(43:46,:));
              
              VD.P0             =double(ChampInt32(Tout(49:52,:)))*10;
              VD.PStd           =double(ChampUInt32(Tout(53:56,:)))*10;
              
    case hex2dec('0100')  % Velocity Data
              VD=double(ChampInt16(Tout(3:end,:)))/1000;
              
    case hex2dec('0200')  % Correlation magnitude data.
              VD=double(Tout(3:end,:));
              
    case hex2dec('0300')  % Echo Intensity data.
              VD=double(Tout(3:end,:))*0.45;
              
    case hex2dec('0400') % Percent good Data.
              VD=double(Tout(3:end,:))/100;
     
    case hex2dec('0500') % ADCP status data.
              VD=Tout(3:end,:);

    case hex2dec('0600') % Bottom Track data.
              VD.BTPings     =ChampUInt16(Tout(3:4,:));
              VD.BTDelay     =ChampUInt16(Tout(5:6,:));
              VD.BTCorrMin   =Tout(7,:);
              VD.BTEvalMin   =Tout(8,:);
              VD.BTPGdMin    =double(Tout(9,:));
              VD.BTMode      =Tout(10,:);
              VD.BTErrVelMax =double(ChampUInt16(Tout(11:12,:)))/1000;
              VD.BTRange     =(double(Tout(78:81,:))*65536+double(ChampUInt16(Tout(17:24,:))))/100;
              VD.BTVel       =double(ChampInt16(Tout(25:32,:)))/1000;
              VD.BTCorr      =Tout(33:36,:);
              VD.BTEval      =Tout(37:40,:);
              VD.BTPGd       =Tout(41:44,:);
              
              VD.BTRefLayMin =double(ChampUInt16(Tout(45:46,:)))/10;
              VD.BTRefLayNear=double(ChampUInt16(Tout(47:48,:)))/10;
              VD.BTRefLayFar =double(ChampUInt16(Tout(49:50,:)))/10;
              VD.BTRefLayVel =double(ChampInt16(Tout(51:58,:)))/1000;
              VD.BTRefLayCor =Tout(59:62,:);
              VD.BTRefLayECI =double(Tout(63:66,:))*0.45;
              VD.BTRefLayPGd =Tout(67:70,:);
              
              VD.BTMaxDepth  =double(ChampUInt16(Tout(71:72,:)))/10;
              VD.BTRSSIAmp   =double(Tout(73:76,:))*0.45;
              VD.BTGain      =double(Tout(77,:));

    otherwise disp('Un type de variables non prÃ©vu dans LitVariable(...)');
end
Variable.VariableData=VD;
return;
            

% Quelques fonctions utilitaires.
function [v]=ChampUInt16(Tout)
v=uint16(Tout(1:2:end,:))+256*uint16(Tout(2:2:end,:));
return;

function [v]=ChampUInt32(Tout)
v=uint32(Tout(1:4:end,:))+256*uint32(Tout(2:4:end,:))+256*256*uint32(Tout(3:4:end,:))+256*256*256*uint32(Tout(4:4:end,:));
return;

function [v]=ChampInt16(Tout)
v=int32(Tout(1:2:end,:))+256*int32(Tout(2:2:end,:));
v=int16(mod(v+32768,65536)-32768);
return;

function [v]=ChampInt32(Tout)
v=double(Tout(1:4:end,:))+256*double(Tout(2:4:end,:))+256*256*double(Tout(3:4:end,:))+256*256*256*double(Tout(4:4:end,:));
v=mod(v+(2^31),2^32)-(2^31);
return;


function [frequency,convex,snsconf,xdcratt,up,angle,bconf] = decode_sysconfig(l, m)

freqs = [75, 150, 300, 600, 1200, 2400, 38];
freqbits = bitand(l, 7);
frequency = freqs(1 + freqbits);
convex    = bitand(bitshift(l,-3),1);
snsconf   = bitand(bitshift(l,-4),3);
xdcratt   = bitand(bitshift(l,-6),1);
up        = bitand(bitshift(l,-7),1);

angles    = [15, 20, 30, NaN];
anglebits = bitand(m, 3);
angle     = angles(1 + anglebits);
bconf     = bitand(bitshift(m,-4),15);
return;




