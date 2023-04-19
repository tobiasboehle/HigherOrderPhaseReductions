function dph = rhs_K1dALL(ph, alpha, delta)
% in the end, the rhs of a (1,infty) phase reduction is given by omega + K*dph, where dph is the output of this function.

N=length(ph);
dph = 1/N*sum((1+delta*sin(ph'))./(1+delta*sin(ph)).*sin(ph'-ph+alpha) - sin(alpha), 2);

end