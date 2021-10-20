%--------------------------------------------------------------------------
function y = diff2(x)
%DIFF2   Centered difference function .  
%        If X is a vector [x(1) x(2) ... x(n)]
%        then DIFF2(X) returns a vector of central differences between
%        elements [x(3)-x(1)  x(4)-x(2) ... x(n)-x(n-2)].  If X is a
%	     matrix, the differences are calculated down each column (i=1:m)

dn = 1;
[m,n] = size(x);
for i=1:m
    y(i,1) = x(i,2) - x(i,1);
    y(i,n) = x(i,n) - x(i,n-1);
    for j=2:n-1
      if isnan(x(i,j+dn)) == 0
        y(i,j) = x(i,j+dn) - x(i,j-dn);     
      else
        y(i,j) = x(i,j) - x(i,j-dn);
      end
    end
end

