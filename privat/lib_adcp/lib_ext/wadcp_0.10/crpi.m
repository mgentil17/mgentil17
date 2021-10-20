function y=crpi(v)
global mes BI_cal
y=sum(((v(1).*BI_cal+v(2))-mes).^2);
