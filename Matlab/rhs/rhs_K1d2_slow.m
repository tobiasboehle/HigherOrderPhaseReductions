function dph = rhs_K1d2_slow(ph, c, alpha)
N = length(ph);

mat = -sin(ph'-3*ph+alpha) + 2*sin(ph'-ph+alpha) - sin(ph+ph'+alpha) - sin(2*ph'-2*ph+alpha) + sin(-2*ph+alpha) + sin(2*ph'+alpha) - sin(alpha);

dph = 1/(4*N*c^2)*sum(mat,2);
    
end