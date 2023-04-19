clear
% set system parameters:

N = 3; %don't change
m = -1; %can be changed
omega = 1; %can be changed
alpha = pi/2+.05; %can be changed

% define variable parameters. the multipliers will be computed for every parameter combination.
K_vec = linspace(-.1, 0.1, 21)';
delta_vec = linspace(0, .4, 5);

% function should not be changed. This script can only compute the multipliers if g=sin.
g = @(phi) sin(phi);
dg = @(phi) cos(phi);

%initialize arrays for both the period of the orbit and its critical multiplier.
period_full  = zeros(length(K_vec), length(delta_vec));
mult_full = zeros(length(K_vec), length(delta_vec));
period_1all  = zeros(length(K_vec), length(delta_vec));
mult_1all = zeros(length(K_vec), length(delta_vec));
period_22    = zeros(length(K_vec), length(delta_vec));
mult_22   = zeros(length(K_vec), length(delta_vec));

optimopt = optimset('TolX', 1e-4);

steps = 20;

parfor kk = 1:length(K_vec) %this can be changed to parfor
    K = K_vec(kk);
    
    %for given value of K, these vectors will contain the periods and critical multipliers for every delta. splitting the arrays in slices like that
    %allows parallel processing
    period_full_K  = zeros(1, length(delta_vec));
    mult_full_K = zeros(1, length(delta_vec));
    period_1all_K  = zeros(1, length(delta_vec));
    mult_1all_K = zeros(1, length(delta_vec));
    period_22_K    = zeros(1, length(delta_vec));
    mult_22_K   = zeros(1, length(delta_vec));

    for dd = 1:length(delta_vec)
        delta = delta_vec(dd);

        period_estimate = 2*pi/(omega-K*sin(alpha));
        ph_splay = 2*pi*linspace(0,1-1/N, N)';

        % full system ===========================================================================================
        ode_rhs_full = @(t,q) rhs_full(q, K, delta, m, omega, g, dg, alpha);
        
        %first search periodic orbit
        %initial conditions for search
        R_ini = ones(N,1);
        ini = [R_ini; ph_splay; period_estimate];
        
        mf = @(x) minfun12_full(ode_rhs_full, x(end), x(N+1:2*N), x(1:N));
        [argmin, fval] = fminsearch(mf, ini, optimopt);
        if fval > 1e-8
            warning("fval too large: %e", fval)
        end
        
        % extract point on periodic orbit:
        R_min = argmin(1:N);
        ph_min = argmin(N+1:2*N);
        
        % get period and multiplier of this periodic orbit.
        period_full_K(dd) = argmin(end);
        mult_full_K(dd) = get_mult_splay(ode_rhs_full, [R_min; ph_min], period_full_K(dd), steps);

        % (1,infty) phase reduction================================================================================================
        % same procedure for (1,infty) phase reduction

        ode_rhs_1all = @(t,ph) omega + K*rhs_K1dALL(ph,alpha,delta);
        ini = [ph_splay; period_estimate];

        mf = @(x) minfun12_Phase(ode_rhs_1all, x(end), x(1:N)) + max(0, abs(period_full_K(dd)-x(end))-1);
        [argmin, fval] = fminsearch(mf , ini, optimopt);
        if fval > 1e-8
            warning("fval too large: %e", fval)
        end

        ph_min = argmin(1:N);
        period_1all_K(dd) = argmin(end);
        mult_1all_K(dd) = get_mult_splay(ode_rhs_1all, ph_min, period_1all_K(dd), steps);

        % (2,2) phase reduction ==================================================================================================
        % same procedure for (2,2) phase reduction
        ode_rhs_22 = @(t,ph) omega + K*rhs_K1d0(ph,alpha) + K*delta*rhs_K1d1(ph,alpha) +K*delta^2*rhs_K1d2(ph,alpha)...
            + K^2*rhs_K2d0(ph,m,alpha) + K^2*delta*rhs_K2d1_slow(ph,m,omega,alpha) + K^2*delta^2*rhs_K2d2_slow(ph,m,omega,alpha);
        ini = [ph_splay; period_estimate];

        mf = @(x) minfun12_Phase(ode_rhs_22, x(end), x(1:N)) + max(0, abs(period_full_K(dd)-x(end))-1);
        [argmin, fval] = fminsearch(mf , ini, optimopt);
        if fval > 1e-8
            warning("fval too large: %e", fval)
        end

        ph_min = argmin(1:N);
        period_22_K(dd) = argmin(end);
        mult_22_K(dd) = get_mult_splay(ode_rhs_22, ph_min, period_22_K(dd), steps);

    end

    period_full(kk,:) = period_full_K;
    mult_full(kk,:) = mult_full_K;
    period_1all(kk,:) = period_1all_K;
    mult_1all(kk,:) = mult_1all_K;
    period_22(kk,:) = period_22_K;
    mult_22(kk,:) = mult_22_K;
end

% omega-K*sin(alpha) is the exact period of the splay orbit when delta=0.
normalize_fct = 2*pi./(omega-K_vec.*sin(alpha));
period_full_normalized = period_full./normalize_fct;
period_1all_normalized = period_1all./normalize_fct;
period_22_normalized = period_22./normalize_fct;

%% plotting period
x_vec = delta_vec;
y_vec = K_vec;

[Y_grid, X_grid] = meshgrid(y_vec, x_vec);

%cel = {period_full, period_1all-period_full, period_22-period_full};
cel = {period_full, period_1all, period_22};
mi = min(cellfun(@(x) min(x, [], 'all'), cel(1:3)));
ma = max(cellfun(@(x) max(x, [], 'all'), cel(1:3)));

fig1 = figure('Name',sprintf("alpha = %g, N=%d, omega = %g", alpha, N, omega));
fig1.Position = [25, 500, 950, 335];

t1 = tiledlayout(1,3);

ax(1) = nexttile;
imagesc(x_vec, y_vec, cel{1})
hold on
contour(X_grid, Y_grid, cel{1}', 'ShowText','on', 'LineColor','black');
%title('Period in Full System')

ax(2) = nexttile;
imagesc(x_vec, y_vec, cel{2})
hold on
contour(X_grid, Y_grid, cel{2}', 'ShowText','on', 'LineColor','black');
%title('Period in 1All System')

ax(3) = nexttile;
imagesc(x_vec, y_vec, cel{3})
hold on
contour(X_grid, Y_grid, cel{3}', 'ShowText','on', 'LineColor','black');
%title('Period in 22 System')

for i = 1:3
    ax(i).YDir = 'normal';
    %xlabel(ax(i), '\delta')
    %ylabel(ax(i), 'K')
    caxis(ax(i), [mi,ma]);
end

col = colorbar;
col.Layout.Tile = 'east';


%===========================================================================================================

%the next paragraph only adjusts the colorbar c, that creates the color gradient from ret to white to green
mi = min(cellfun(@(x) min(x, [], 'all'), {mult_full, mult_1all, mult_22}));
ma = max(cellfun(@(x) max(x, [], 'all'), {mult_full, mult_1all, mult_22}));
space_above = ma-1;
space_below = 1-mi;
spacing = .0001; % ..., 1-spacing/2, 1+spacing/2, 1+3spacing/2, ...
tiles_above = ceil(space_above/spacing-.5);
tiles_below = ceil(space_below/spacing-.5);
tiles = tiles_above + tiles_below + 1;
mi = 1 - spacing*(tiles_below + .5);
ma = 1 + spacing*(tiles_above + .5);
c_shift = @(x) x.^7;
c_step = 1/(max(tiles_below, tiles_above));
c_below = [c_shift(linspace(1-c_step*tiles_below, 1, tiles_below+1)'), ones(tiles_below+1,1), c_shift(linspace(1-c_step*tiles_below, 1, tiles_below+1)')];
c_above = [ones(tiles_above+1,1), c_shift(linspace(1,1-tiles_above*c_step, tiles_above+1)'), c_shift(linspace(1,1-tiles_above*c_step, tiles_above+1)')];
c_below(end,:) = [];
c_above(1,:) = [];
c = [c_below; ones(1,3); c_above];

fig2 = figure('Name',sprintf("alpha = %g, N=%d, omega = %g", alpha, N, omega));
fig2.Position = fig1.Position;
fig2.Position(2) = 50;
t2 = tiledlayout(1,3);

ax(4) = nexttile;
imagesc(x_vec, y_vec, mult_full)
hold on
contour(X_grid, Y_grid, mult_full', 'ShowText','on', 'LineColor','black');
%title('Floquet in Full System')

ax(5) = nexttile;
imagesc(x_vec, y_vec,  mult_1all)
hold on
contour(X_grid, Y_grid, mult_1all', 'ShowText','on', 'LineColor','black');
%title('Floquet in 1All System')

ax(6) = nexttile;
imagesc(x_vec, y_vec,  mult_22)
hold on
contour(X_grid, Y_grid, mult_22', 'ShowText','on', 'LineColor','black');
%title('Floquet in 22 System')


for i = 4:6
    ax(i).YDir = 'normal';
    %xlabel(ax(i), '\delta')
    %ylabel(ax(i), 'K')
    caxis(ax(i), [mi,ma]);
    colormap(ax(i), c);
end

col = colorbar;
col.Layout.Tile = 'east';