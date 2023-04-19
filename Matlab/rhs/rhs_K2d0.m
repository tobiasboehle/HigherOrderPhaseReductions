function dph = rhs_K2d0(ph, m, alpha)
%compute P^{2,0}, as defined in the paper.
% if N is large this function is faster than "rhs_K2d0_slow", but the result is the same.

Z = mean(exp((1:2).*1i.*ph), 1);
r = abs(Z);
ang = angle(Z);

dph = (-1).*r(1).*sin(2.*alpha+(-1).*ph+ang(1))+r(1).^2.*sin(2.*alpha+( ...
  -2).*ph+2.*ang(1))+r(1).*r(2).*sin(ph+ang(1)+(-1).*ang(2));
dph = dph/2/m;

    
end