function dph = rhs_K1dALL(ph, alpha, delta)
N=length(ph);
dph = 1/N*sum((1+delta*sin(ph'))./(1+delta*sin(ph)).*sin(ph'-ph+alpha) - sin(alpha), 2);

end