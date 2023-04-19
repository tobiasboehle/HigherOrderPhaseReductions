function dph = rhs_K1d0(ph, alpha)

%compute P^{1,0}, as defined in the paper.

Z = mean(exp(1i.*ph), 1);
op = abs(Z);
Psi = angle(Z(1));

dph = op*sin(Psi-ph+alpha)-sin(alpha);
    
end