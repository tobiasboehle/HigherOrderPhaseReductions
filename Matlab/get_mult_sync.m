function mult = get_mult_sync(rhs, drhs, dim, times, omega)
epsilon = 1e-8;
delta_t = mean(diff(times));

transition_matrix = eye(dim);
for i = 1:length(times)
    xi = times(i)*omega; %position on the periodic orbit (xi,xi,xi,...,xi)

    %next calculate derivative of rhs at position on periodic orbit:
    if ~isempty(drhs)
        A = drhs(xi);
    else
        A = zeros(dim);
        ph = xi*ones(dim,1);
        center = rhs(ph);
        for j = 1:dim
            ej = zeros(dim,1); ej(j)=1;
            A(:,j) = (rhs(ph+epsilon*ej)-center)/epsilon;
        end
    end
    
    %update of transition matrix with the matrix exponential
    transition_matrix = expm(A*delta_t)*transition_matrix;
end
% now transition_matrix approximates the monodromy matrix, whose eigenvalues are 1,mu_1,...,mu_{N-1}, where mu_1,\dots,mu_{N-1} are the multipliers
e = abs(eig(transition_matrix));
logi = abs(e-1)<1e-6; %exclude 1
if ~isempty(logi)
    e(find(logi, 1))=[];
end
mult = max(abs(e));

end