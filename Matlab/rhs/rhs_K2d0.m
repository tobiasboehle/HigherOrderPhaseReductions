function dph = rhs_K2d0(ph, m, alpha)
    
Z = mean(exp((1:2).*1i.*ph), 1);
r = abs(Z);
ang = angle(Z);

% op = abs(Z);
% Psi = angle(Z(1));
% Theta = angle(Z(2));
% dph = 1/2/m*( -op(1)*sin(Psi-ph+2*alpha) + op(1)*op(2)*sin(ph-Theta+Psi) + op(1)^2*sin(2*Psi-2*ph+2*alpha) );

dph = (-1).*r(1).*sin(2.*alpha+(-1).*ph+ang(1))+r(1).^2.*sin(2.*alpha+( ...
  -2).*ph+2.*ang(1))+r(1).*r(2).*sin(ph+ang(1)+(-1).*ang(2));
dph = dph/2/m;

    
end