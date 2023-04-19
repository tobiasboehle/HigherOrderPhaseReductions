function dph = rhs_K1d2_slow(ph, alpha)
% compute P^{1,2}, as defined in the paper.
% if N is large this function is slower than "rhs_K1d2", but the result is the same.

N = length(ph);

mat = -sin(ph'-3*ph+alpha) + 2*sin(ph'-ph+alpha) - sin(ph+ph'+alpha) - sin(2*ph'-2*ph+alpha) + sin(-2*ph+alpha) + sin(2*ph'+alpha) - sin(alpha);

dph = 1/(4*N)*sum(mat,2);
    
end