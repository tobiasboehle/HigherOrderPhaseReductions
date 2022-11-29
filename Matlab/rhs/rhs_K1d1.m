function dph = rhs_K1d1(ph, alpha)
    
Z = mean(exp((1:2).*1i.*ph), 1);
r = abs(Z);
ang = angle(Z);

dph = 1/2*(-r(1)*cos(-2*ph+ang(1)+alpha)-r(2)*cos(ang(2)-ph+alpha) + cos(-ph+alpha) + r(1)*cos(ang(1)+alpha));
    
end