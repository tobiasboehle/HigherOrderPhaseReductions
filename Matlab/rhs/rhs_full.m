function out = rhs_full(q, K, delta, m, omega, g, dg, alpha)
%this function computes the rhs of the full/unreduced system.

% Inputs:
% q: (2N,1)-array containing R and phi
% K, delta, m, omega, alpha: parameters of the system
% g, dg: functions

%extract R and phi from q
N = length(q)/2;
mat = reshape(q, [N,2]);
Rk = mat(:,1); %this is R_k as in the paper
phi = mat(:,2);

%this is R_l as in the definition of the unreduced system in the paper.
Rl = Rk.';

lim_k = (1+delta*g(phi));
lim_l = lim_k.';

% functions F,G,H that define the rhs
F = m*Rk.^2.*(Rk-1).*lim_k.^2;
G = 1/N*sum(Rl.*lim_l./lim_k.*cos(phi.'-phi+alpha)-Rk*cos(alpha) - delta*dg(phi).*(Rl.*lim_l./lim_k.^2.*sin(phi.'-phi+alpha)-Rk.*sin(alpha)./lim_k), 2);
H = 1/N*sum(Rl.*lim_l./(Rk.*lim_k).*sin(phi.'-phi+alpha)-sin(alpha), 2);

%calculate rhs
dR = F + K*G;
dp = omega + K*H;

out = [dR; dp];

end
