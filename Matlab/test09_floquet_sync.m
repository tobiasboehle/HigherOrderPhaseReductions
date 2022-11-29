N = 3;

alpha_vec = pi/2+.05;
K_vec = linspace(-.1,0.1, 51);
delta_vec = linspace(0,.4, 51);

% K_vec = .1;
% delta_vec = .1;

m = -1;
omega = 1;

%[K_grid, alpha_grid] = meshgrid(K_vec, alpha_vec);

g = @(phi) sin(phi);
dg = @(phi) cos(phi);
%ddg = @(phi) -sin(phi);

%opt = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);

sz = [length(alpha_vec), length(K_vec), length(delta_vec)];
max_floquet_full_ind = zeros(prod(sz),1);
max_floquet_11_ind   = zeros(prod(sz),1);
max_floquet_1all_ind = zeros(prod(sz),1);
max_floquet_21_ind   = zeros(prod(sz),1);
max_floquet_22_ind   = zeros(prod(sz),1);

floq_fun_1all = @(K, alpha, omega) exp(-2*pi*K/omega*cos(alpha));
floq_fun_22 = @(K,alpha,omega,delta) exp(-2*pi*K/(m*omega*(m^2+omega^2))...
    *(m^3*cos(alpha) + m*omega^2 *cos(alpha) - K*m^2*sin(alpha)^2 - 2*K*m^2*delta^2*sin(alpha)^2-K*omega^2*sin(alpha)^2));

steps = 20;
sync_period = 2*pi/omega;
delta_t = sync_period/steps;
times = linspace(delta_t/2, sync_period-delta_t/2, steps);

epsilon = 1e-8;

tic

parfor i = 1:prod(sz)
    [aa,kk,dd] = ind2sub(sz, i);
    alpha = alpha_vec(aa);
    K = K_vec(kk);
    delta = delta_vec(dd);

    % full system=====================================================
    max_floquet_full_ind(i) = floquet_sync_full_system(m,alpha,omega,K,delta,g,dg,2*N, times);

    switch "numeric"
        case "numeric"
            % 1st order approximation=====================================================================
            rhs1order = @(ph) omega + K*rhs_K1d0(ph,alpha) + K*delta*rhs_K1d1(ph,alpha);
            max_floquet_11_ind(i) = get_floquet_sync(rhs1order, [], N, times, omega);

            % 1st order approximation, all delta=====================================================================
            rhs1order = @(ph) omega + K*rhs_K1dALL(ph,alpha, delta);
            max_floquet_1all_ind(i) = get_floquet_sync(rhs1order, [], N,times, omega);

            %2nd order approximation, 1st in delta========================================================
            rhs2order = @(ph) omega + K*rhs_K1d0(ph,alpha) + K*delta*rhs_K1d1(ph,alpha) + K^2*rhs_K2d0(ph,m,alpha) + K^2*delta*rhs_K2d1(ph,m,omega,alpha);
            max_floquet_21_ind(i) = get_floquet_sync(rhs2order, [], N,times,omega);

            %2nd order approximation, 2nd in delta ========================================================
            rhs2order = @(ph) omega + K*rhs_K1d0(ph,alpha) + K*delta*rhs_K1d1(ph,alpha) +K*delta^2*rhs_K1d2(ph,alpha)...
                + K^2*rhs_K2d0(ph,m,alpha) + K^2*delta*rhs_K2d1(ph,m,omega,alpha) + K^2*delta^2*rhs_K2d2(ph,m,omega,alpha);
            max_floquet_22_ind(i) = get_floquet_sync(rhs2order, [], N,times,omega);

        case "analytic"
            max_floquet_11_ind(i) = floq_fun_1all(K,alpha,omega);
            max_floquet_1all_ind(i) = floq_fun_1all(K,alpha,omega);
            max_floquet_21_ind(i) = floq_fun_22(K,alpha,omega,0);
            max_floquet_22_ind(i) = floq_fun_22(K,alpha,omega,delta);
    end
end
toc

max_floquet_full = reshape(max_floquet_full_ind, sz);
max_floquet_11   = reshape(max_floquet_11_ind, sz);
max_floquet_1all = reshape(max_floquet_1all_ind, sz);
max_floquet_21   = reshape(max_floquet_21_ind, sz);
max_floquet_22   = reshape(max_floquet_22_ind, sz);


%%
aa = 1;
fig1 = figure('Name',sprintf("alpha = %g, N=%d, omega = %g", alpha_vec(aa), N, omega));
fig1.Position(1) = 25;
fig1.Position(3) = 900;

delta_cross = [0];
[~,delta_cross_ind] = min(abs(delta_cross - delta_vec), [], 2);
delta_cross = delta_vec(delta_cross_ind);
coll = [0.9294    0.6941    0.1255];

x_vec = delta_vec;
y_vec = K_vec;

[Y_grid, X_grid] = meshgrid(y_vec, x_vec);

data_full = squeeze(max_floquet_full(aa,:,:));
data_1all = squeeze(max_floquet_1all(aa,:,:));
data_22 = squeeze(max_floquet_22(aa,:,:));
data_cell = {data_full, data_1all, data_22};
ma = max(cellfun(@(x) max(x, [], 'all'), data_cell));
mi = min(cellfun(@(x) min(x, [], 'all'), data_cell));

% ma = max(data_full, [], 'all');
% mi = min(data_full, [], 'all');
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
% c_below = [c_shift(linspace(0,1,tiles_below+1)'), ones(tiles_below+1,1), c_shift(linspace(0,1,tiles_below+1)')];
% c_above = [ones(tiles_above+1,1), c_shift(linspace(1,0, tiles_above+1)'), c_shift(linspace(1,0, tiles_above+1)')];
c_below(end,:) = [];
c_above(1,:) = [];
c = [c_below; ones(1,3); c_above];

% ax = arrayfun(@(i) subplot(1,4,i), 1:3);
% ax(4) = subplot(2,4,4);
% ax(5) = subplot(2,4,8);
t = tiledlayout(1,3);

axc = 1; %axes counter

%=======================================================================
ax(axc) = nexttile; axc = axc+1;
imagesc(x_vec, y_vec, data_full)
hold on
contour(X_grid, Y_grid, data_full', 'ShowText','on', 'LineColor','black');
for i = 1:length(delta_cross)
    plot(delta_cross(i)*ones(2,1), [min(K_vec);max(K_vec)], 'Color',  .5*ones(3,1), 'LineWidth',2);
end
%=======================================================================
ax(axc) = nexttile; axc = axc+1;

data = data_1all;
imagesc(x_vec, y_vec, data)
hold on
contour(X_grid, Y_grid, data', 'ShowText','on', 'LineColor','black');
for i = 1:length(delta_cross)
    plot(delta_cross(i)*ones(2,1), [min(K_vec);max(K_vec)], 'Color',  coll, 'LineStyle', '--', 'LineWidth',2);
end
% incorrect = xor(data_1all > 1, data_full >1);
% K_start = min(K_vec)*2;
% delta_start = min(delta_vec);
% K_step = .01;
% ang = pi/8;
% t = linspace(0,.5,1000);
% for ll = 1:100
%     x = delta_start + cos(ang)*t;
%     y = K_start + K_step*(ll-1) + sin(ang)*t;
%     [~,xmind] = min(abs(x-delta_vec'));
%     [~,ymind] = min(abs(y-K_vec'));
%     linind = sub2ind([length(K_vec), length(delta_vec)], ymind, xmind);
%     logi = incorrect(linind);
%     x(~logi) = NaN;
%     plot(x,y, 'blue', 'LineWidth', 1.5)
% end
%title("Floquet Diff K1dALL")

%==============================================================
ax(axc) = nexttile; axc = axc+1;
data = data_22;
imagesc(x_vec, y_vec, data)
hold on
contour(X_grid, Y_grid, data', 'ShowText','on', 'LineColor','black');
for i = 1:length(delta_cross)
    plot(delta_cross(i)*ones(2,1), [min(K_vec);max(K_vec)], 'Color',  coll, 'LineStyle', '-.', 'LineWidth',2);
end
%title("Floquet Diff K2d2")

col = colorbar;
col.Layout.Tile = 'south';

for i = 1:3
    ax(i).YDir = 'normal';
    caxis(ax(i), [mi,ma]);
    colormap(ax(i), c);
end
linkaxes(ax(1:3));

%% ==============================================================
fig2 = figure('Name',sprintf("alpha = %g, N=%d, omega = %g", alpha_vec(aa), N, omega));
fig2.Position(1) = 950;
fig2.Position(3) = 300;

t2 = tiledlayout(length(delta_cross),1);

for i = 1:length(delta_cross)
    ax(3+i) = nexttile;
    %coll = uisetcolor;
    plot(K_vec, data_full(:,delta_cross_ind(1)),'Color', .5*ones(3,1), 'DisplayName', 'Full System', 'LineWidth', 3);
    hold on
    plot(K_vec, data_1all(:,delta_cross_ind(1)), 'LineStyle', '--', 'Color', coll, 'DisplayName', 'First Order', 'LineWidth', 1.5);
    plot(K_vec, data_22(:,delta_cross_ind(1)), 'LineStyle', '-.', 'Color', coll,'DisplayName', 'Second Order', 'LineWidth', 1.5);
    grid on
end

legend(ax(4), 'Location', 'south');
