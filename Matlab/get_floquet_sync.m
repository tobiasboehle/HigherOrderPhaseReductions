function f = get_floquet_sync(rhs, drhs, dim, times, omega)
epsilon = 1e-8;
delta_t = mean(diff(times));

transition_matrix = eye(dim);
for i = 1:length(times)
    xi = times(i)*omega;
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
    transition_matrix = expm(A*delta_t)*transition_matrix;
end
e = abs(eig(transition_matrix));
logi = abs(e-1)<1e-6;
if ~isempty(logi)
    e(find(logi, 1))=[];
end
f = max(abs(e));

end