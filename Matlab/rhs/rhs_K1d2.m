function dph = rhs_K1d2(ph, alpha)

%compute P^{1,2}, as defined in the paper.
% if N is large this function is faster than "rhs_K1d2_slow", but the result is the same.

Z = mean(exp((1:2).*1i.*ph), 1);
r = abs(Z);
ang = angle(Z);

dph = 1/4*(-r(1)*sin(ang(1)-3*ph+alpha) + 2*r(1)*sin(ang(1)-ph+alpha) - r(1)*sin(ph+ang(1)+alpha) - r(2)*sin(ang(2)-2*ph+alpha) + sin(-2*ph+alpha) + r(2)*sin(ang(2)+alpha) - sin(alpha));
    
end