function f = get_floquet_splay(ode_rhs, ini, period, steps)
dim = length(ini);
epsilon = 1e-8;
times = period * linspace(1/2/steps, 1-1/2/steps, steps);
delta_t = mean(diff(times));

opt = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);
sol = ode45(ode_rhs, [0,period], ini, opt);


transition_matrix = eye(dim);
for i = 1:steps
    state = deval(sol, times(i));
    A = zeros(dim);
    center = ode_rhs(times(i), state);
    for j = 1:dim
        ej = zeros(dim,1); ej(j)=1;
        A(:,j) = (ode_rhs(times(i), state+epsilon*ej)-center)/epsilon;
    end
    transition_matrix = expm(A*delta_t)*transition_matrix;
end
e = abs(eig(transition_matrix));
[mi,ind]=min(abs(e-1));
if mi > 1e-4
    warning("no eigenvalue 1")
end
e(ind)=[];
f = max(abs(e));
end