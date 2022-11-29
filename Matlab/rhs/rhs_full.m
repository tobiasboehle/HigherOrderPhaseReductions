function out = rhs_full(q, K, delta, m, omega, g, dg, alpha)
N = length(q)/2;
mat = reshape(q, [N,2]);
Rk = mat(:,1);
phi = mat(:,2);

Rl = Rk.';

lim_k = (1+delta*g(phi));
lim_l = lim_k.';

F = m*Rk.^2.*(Rk-1).*lim_k.^2;
G = 1/N*sum(Rl.*lim_l./lim_k.*cos(phi.'-phi+alpha)-Rk*cos(alpha) - delta*dg(phi).*(Rl.*lim_l./lim_k.^2.*sin(phi.'-phi+alpha)-Rk.*sin(alpha)./lim_k), 2);
H = 1/N*sum(Rl.*lim_l./(Rk.*lim_k).*sin(phi.'-phi+alpha)-sin(alpha), 2);

dR = F + K*G;
dp = omega + K*H;

out = [dR; dp];

end
