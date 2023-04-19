function dph = rhs_K2d0_slow(ph, m, alpha)
% compute P^{2,0}, as defined in the paper.
% if N is large this function is slower than "rhs_K2d0", but the result is the same.

N = length(ph);
R01 = 1/(N*m)*sum(-cos(ph'-ph+alpha)+cos(alpha), 2);

grad0H = 1/N*sin(ph-ph'+alpha);
grad0H = grad0H - diag(diag(grad0H));
grad0H = grad0H + diag(-sum(grad0H,1));

dph = grad0H'*R01;



end