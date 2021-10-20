function [BI] = backscatter_corr(RL, Kc, Rv, T, S, f, P, method)
  % RL (counts)		: Received level
  % Kc (dB/counts)	: Count to dB conversion factor
  % Rv (m)		: Acoustic path = R/cos(theta)
  % T (°C)		: Temperature
  % S (PSU)		: Salinity
  % f (kHz)		: Frequency
  % P (m)		: Sensors depth
  % method		: BI estimation method deines or gostiaux
  % BI (dB)		: Backscatter index
  
  % Seawater absorption loss
  alpha_wdb = equ_att_son_garrison(f,P,T,S);% dB/m, seawater absorption
  a_wdb     = repmat(alpha_wdb, 1, length(Rv));
  TL_w      = 2*a_wdb .* repmat(Rv, size(a_wdb,1), 1);
  
  % Geometrical spreading
  tl_g = 20*log10(Rv);
  TL_g = repmat(tl_g, size(TL_w, 1), 1);
  
  % Noise level
  Er=min(RL(:));%P.Cauchy

  % Backscatter Indice
  if strcmp(method,'gostiaux')
    BI = 10*log10(10.^(Kc*RL/10) - 10^(Kc*Er/10)) + TL_w + TL_g;
  elseif strcmp(method,'deines')
    BI = Kc*(RL-Er) + TL_w + TL_g;
  elseif strcmp(method,'trdi')%Mulison (2017)
    BI = 10*log10(10.^(Kc*(RL-Er)/10)-1)+ TL_w + TL_g;
  else
    disp(['WARNING : Unknown Backscatter estimation method : ' method '.'])
    disp(' deines used instead')
    BI = Kc*(RL-Er) + TL_w + TL_g;
  end
end
