function dph = rhs_K2d1_slow(ph, m, omega, alpha)
N = length(ph);

fct0 =  1/(N*m);
R01 = fct0*sum(-cos(ph'-ph+alpha)+cos(alpha), 2);

c = 1;

fct1 = 1/(2*N*(m^2+omega^2));

mat1 = (-2).*omega.*cos(alpha+(-1).*ph)+(-3).*omega.*cos(alpha+ph)+4.* ...
  omega.*cos(alpha+transpose(ph))+2.*omega.*cos(alpha+(-2).*ph+ ...
  transpose(ph))+(-1).*omega.*cos(alpha+(-1).*ph+2.*transpose(ph))+ ...
  2.*m.*sin(alpha+(-1).*ph)+(-3).*m.*sin(alpha+ph)+4.*m.*sin(alpha+ ...
  transpose(ph))+(-2).*m.*sin(alpha+(-2).*ph+transpose(ph))+(-1).* ...
  m.*sin(alpha+(-1).*ph+2.*transpose(ph));
R11 = fct1*sum(mat1, 2);


grad0H = 1/N*sin(ph-ph'+alpha);
grad0H = grad0H - diag(diag(grad0H));
grad0H = grad0H + diag(-sum(grad0H,1));

grad1H = 1/N/c * (sin(ph)-sin(ph')).*sin(ph-ph'+alpha);
grad1H = grad1H - diag(diag(grad1H));
grad1H = grad1H + diag(-sum(grad1H,1));

%blue
dph1 = grad0H'*R11;

%green
dph2 = grad1H'*R01;

dph = dph1+dph2;
    
end