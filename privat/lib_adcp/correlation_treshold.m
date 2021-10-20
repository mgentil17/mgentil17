function [mask, scorr] = correlation_treshold(tresh, corr, R, alt)

mask  = NaN*corr;
scorr = NaN*corr;
for i=1:size(corr, 1)
  if ~isempty(find(R(i, :) < alt(i)))
    ib = max(find(R(i, :) < alt(i)));
  else
    ib = size(R, 2);
  end
  innan = find(~isnan(corr(i, :)));
  if length(innan)>10
    tmp         = filtfilt(ones(1, 4)/4, 1, corr(i, innan));
    scorr(i, :) = interp1(innan, tmp, 1:size(corr, 2));
  end
  itresh = find(scorr(i, 1:ib)>tresh);
  mask(i, itresh) = ones(1, length(itresh));
end
end
