% set fixed system parameters
N = 3; %this is the number of oscillators
m = -1;
omega = 1;

% the parameters alpha, K, delta can take arbitrary values in the following vectors. The multipliers are then computed for every combination of
% parameters. The final figure will only be displayed for one value of alpha.
alpha_vec = pi/2+.05;
K_vec = linspace(-.1,0.1, 21);
delta_vec = linspace(0,.4, 21);


% Define the function g and its derivative dg. These definitions should not be changed! This script can only calculate multipliers when g(phi) = sin(phi).
g = @(phi) sin(phi);
dg = @(phi) cos(phi);


sz = [length(alpha_vec), length(K_vec), length(delta_vec)];
% in total there are prod(sz) combinations of parameters for which we compute critical multipliers. We compute them for the full system, the (1,1) phase
% reduction, the (1,\infty) phase reduction, the (2,1) phase reduction and the (2,2) phase reduction. The following 5 lines of code initialize vectors
% that will eventually contain the computation result.
% Remark: one could of course also use 3d-arrays with sizes as in sz. The use of vectors, however, allows parallel computing.
crit_mult_full_ind = zeros(prod(sz),1);
crit_mult_11_ind   = zeros(prod(sz),1);
crit_mult_1all_ind = zeros(prod(sz),1);
crit_mult_21_ind   = zeros(prod(sz),1);
crit_mult_22_ind   = zeros(prod(sz),1);

% the following two functions are analytical results about the multipliers for the (1,\infty) and the (2,2) phase reduction.
mult_fun_1all = @(K, alpha, omega) exp(-2*pi*K/omega*cos(alpha));
mult_fun_22 = @(K,alpha,omega,delta) exp(-2*pi*K/(m*omega*(m^2+omega^2))...
    *(m^3*cos(alpha) + m*omega^2 *cos(alpha) - K*m^2*sin(alpha)^2 - 2*K*m^2*delta^2*sin(alpha)^2-K*omega^2*sin(alpha)^2));

% to numerically calculate PRMMs, we divide the periodic orbit in a few steps. See functions "mult_sync_full_system" and "get_mult_sync" for details.
steps = 20;
sync_period = 2*pi/omega; %period of synchronized orbit.
delta_t = sync_period/steps;
times = linspace(delta_t/2, sync_period-delta_t/2, steps); %discrete times, where we divide the periodic orbit.

%epsilon = 1e-8;

tic

for i = 1:prod(sz) %this can be replaced by parfor for parallel processing
    % for a linear index i, get parameters alpha, K and delta.
    [aa,kk,dd] = ind2sub(sz, i);
    alpha = alpha_vec(aa);
    K = K_vec(kk);
    delta = delta_vec(dd);

    % full system=====================================================
    % numerically calculate the critical PRMM for the full/unreduced system.
    crit_mult_full_ind(i) = mult_sync_full_system(m,alpha,omega,K,delta,g,dg,2*N, times);

    switch "numeric" %"analytic"
                    %there are two ways how to compute multipliers. One is numerically, by approximating the monodromy matrix and then calculating its eigenvalues.
                    % the other way is analytically, by the results from the paper. Choose "numeric" or "analytic" for one of these options. The
                    % result should be independent of this choice. Computation time is much better for "analytic".
        case "numeric"
            % 1st order approximation=====================================================================
            % numerically calculate the critical PRMM for the (1,1) phase reduction
            rhs1order = @(ph) omega + K*rhs_K1d0(ph,alpha) + K*delta*rhs_K1d1(ph,alpha);
            crit_mult_11_ind(i) = get_mult_sync(rhs1order, [], N, times, omega);

            % 1st order approximation, all delta=====================================================================
            % numerically calculate the critical PRMM for the (1,\infty) phase reduction.
            rhs1order = @(ph) omega + K*rhs_K1dALL(ph,alpha, delta);
            crit_mult_1all_ind(i) = get_mult_sync(rhs1order, [], N,times, omega);

            %2nd order approximation, 1st in delta========================================================
            % numerically calculate the critical PRMM for the (2,1) phase reduction.
            rhs2order = @(ph) omega + K*rhs_K1d0(ph,alpha) + K*delta*rhs_K1d1(ph,alpha) + K^2*rhs_K2d0(ph,m,alpha) + K^2*delta*rhs_K2d1(ph,m,omega,alpha);
            crit_mult_21_ind(i) = get_mult_sync(rhs2order, [], N,times,omega);

            %2nd order approximation, 2nd in delta ========================================================
            %numerically calculate the critical PRMM for the (2,2) phase reduction.
            rhs2order = @(ph) omega + K*rhs_K1d0(ph,alpha) + K*delta*rhs_K1d1(ph,alpha) +K*delta^2*rhs_K1d2(ph,alpha)...
                + K^2*rhs_K2d0(ph,m,alpha) + K^2*delta*rhs_K2d1(ph,m,omega,alpha) + K^2*delta^2*rhs_K2d2(ph,m,omega,alpha);
            crit_mult_22_ind(i) = get_mult_sync(rhs2order, [], N,times,omega);

        case "analytic"
            % this uses analytic results to calculate the critical PRMMs
            crit_mult_11_ind(i) = mult_fun_1all(K,alpha,omega);
            crit_mult_1all_ind(i) = mult_fun_1all(K,alpha,omega);
            crit_mult_21_ind(i) = mult_fun_22(K,alpha,omega,0); %agrees with (2,2) phase reduction when delta = 0 (shown in paper).
            crit_mult_22_ind(i) = mult_fun_22(K,alpha,omega,delta);
    end
end
toc

%now reshape the vectors to 3d-arrays that contain the final results.
crit_mult_full = reshape(crit_mult_full_ind, sz);
crit_mult_11   = reshape(crit_mult_11_ind, sz);
crit_mult_1all = reshape(crit_mult_1all_ind, sz);
crit_mult_21   = reshape(crit_mult_21_ind, sz);
crit_mult_22   = reshape(crit_mult_22_ind, sz);


%% Plotting Figure1

aa = 1; %Plotted results only include the aa-th value of alpha in alpha_vec
%create figure 1
fig1 = figure('Name',sprintf("alpha = %g, N=%d, omega = %g", alpha_vec(aa), N, omega));
% adjust position of figure 1. This can be omitted
fig1.Position(1) = 25;
fig1.Position(3) = 900;

delta_cross = 0; % a value of delta, for which a cross section will be plotted in figure 2. This cross section will be indicated by a orange line in figure 1.
[~,delta_cross_ind] = min(abs(delta_cross - delta_vec), [], 2);
delta_cross = delta_vec(delta_cross_ind); %adjusts delta_cross to closest point in delta_vec.
coll = [0.9294    0.6941    0.1255]; %set color.

x_vec = delta_vec; %define vector vor x-axis
y_vec = K_vec; %define vector for y-axis
[Y_grid, X_grid] = meshgrid(y_vec, x_vec); %grids required for surf plots

% find max and min of all data
data_full = squeeze(crit_mult_full(aa,:,:));
data_1all = squeeze(crit_mult_1all(aa,:,:));
data_22 = squeeze(crit_mult_22(aa,:,:));
data_cell = {data_full, data_1all, data_22};
ma = max(cellfun(@(x) max(x, [], 'all'), data_cell));
mi = min(cellfun(@(x) min(x, [], 'all'), data_cell));

%next paragraph is only there to create the colormap c, which defines the color gradient from red to white to green
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


t = tiledlayout(1,3);
axc = 1; %axes counter

% plot PRMMs of full/unreduced  system =======================================================================
ax(axc) = nexttile; axc = axc+1;
imagesc(x_vec, y_vec, data_full)
hold on
contour(X_grid, Y_grid, data_full', 'ShowText','on', 'LineColor','black');
for i = 1:length(delta_cross)
    plot(delta_cross(i)*ones(2,1), [min(K_vec);max(K_vec)], 'Color',  .5*ones(3,1), 'LineWidth',2);
end


% plot PRMM of (1,\infty) phase reduction =======================================================================
ax(axc) = nexttile; axc = axc+1;
data = data_1all;
imagesc(x_vec, y_vec, data)
hold on
contour(X_grid, Y_grid, data', 'ShowText','on', 'LineColor','black');
for i = 1:length(delta_cross)
    plot(delta_cross(i)*ones(2,1), [min(K_vec);max(K_vec)], 'Color',  coll, 'LineStyle', '--', 'LineWidth',2);
end

% plot PRMM of (2,2) phase reduction ==============================================================
ax(axc) = nexttile; axc = axc+1;
data = data_22;
imagesc(x_vec, y_vec, data)
hold on
contour(X_grid, Y_grid, data', 'ShowText','on', 'LineColor','black');
for i = 1:length(delta_cross)
    plot(delta_cross(i)*ones(2,1), [min(K_vec);max(K_vec)], 'Color',  coll, 'LineStyle', '-.', 'LineWidth',2);
end

%general settings
col = colorbar;
col.Layout.Tile = 'south';
for i = 1:3
    ax(i).YDir = 'normal';
    caxis(ax(i), [mi,ma]);
    colormap(ax(i), c);
end
linkaxes(ax(1:3));

%% Plotting Figure 2
% create figure 2, assign name and adjust location.
fig2 = figure('Name',sprintf("alpha = %g, N=%d, omega = %g", alpha_vec(aa), N, omega));
fig2.Position(1) = 950;
fig2.Position(3) = 300;

t2 = tiledlayout(length(delta_cross),1);

%plots the PRMMs for a selected value of delta.
ax(4) = nexttile;
plot(K_vec, data_full(:,delta_cross_ind),'Color', .5*ones(3,1), 'DisplayName', 'Full System', 'LineWidth', 3);
hold on
plot(K_vec, data_1all(:,delta_cross_ind), 'LineStyle', '--', 'Color', coll, 'DisplayName', 'First Order', 'LineWidth', 1.5);
plot(K_vec, data_22(:,delta_cross_ind), 'LineStyle', '-.', 'Color', coll,'DisplayName', 'Second Order', 'LineWidth', 1.5);
grid on


legend(ax(4), 'Location', 'south');
