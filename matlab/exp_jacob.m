clear; close all; clc;

b_km = Constants.EARTH_R_km + linspace(0, 100, 101)';   % [B x 1]

% Atmosphere model
atmo = generate_atmosphere_model();

% Curve arc-length parameter grid
s_km = (0:10:2800)';    % [S x 1]

% Create radical refractive functions
nr_func = create_nr_func(550, atmo);

% Integrate and get the raypath
fprintf('Computing direct ray paths...');
raypath_store = compute_raypath(nr_func, b_km, s_km); % [B x S x 6]
fprintf('done!\n');

ray_angle = squeeze(atan2(raypath_store(:,end,6), raypath_store(:,end,5)));

% Solve for sun and moon position
sun_sol = solve_raypath(raypath_store, Constants.EARTH_SUN_DIST_km);
sun_sol(:,1) = -sun_sol(:,1);
sun_sol(:,4) = mod(pi - sun_sol(:,4), 2*pi);
moon_sol = solve_raypath(raypath_store, Constants.EARTH_MOON_DIST_km);

% Compute angles
sun_emit = pi - sun_sol(:,4) - ray_angle;
sun_mid = pi - sun_sol(:,4);
sun_exit = pi - sun_sol(:,4) + ray_angle;

% Interpolate and compute derivatives
pp_mid_emit = fnder(pchip(sun_emit, sun_mid), 1);
pp_exit_mid = fnder(pchip(sun_mid, sun_exit), 1);
pp_exit_emit = fnder(pchip(sun_emit, sun_exit), 1);

%%
figure(1); clf;
plot(sun_emit, sun_exit, 'o-', ...
    sun_emit, sun_mid, '^-');

% Results look solid, which is reasonable. It is basically a direct deduction of chain rule.
figure(2); clf;
plot(sun_emit, ppval(pp_exit_emit, sun_emit), 'o-', ...
    sun_emit, ppval(pp_mid_emit, sun_emit), '.-', ...
    sun_emit, ppval(pp_mid_emit, sun_emit) .* ppval(pp_exit_mid, sun_emit), '*-');
