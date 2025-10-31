clear; close all; clc;

EARTH_R_km = 6371;
EARTH_SUN_DIST_km = 149.6e6;
EARTH_MOON_DIST_km = 384.4e3;
SUN_R_km = 696e3;
SUN_R_arcmin = atand(SUN_R_km / EARTH_SUN_DIST_km) * 60;   % ~16 arcmin

% Atmosphere model
z_max = 90;
z_grid_km = (0:.1:z_max)';
atmo = us_standard_atmosphere_1976(z_grid_km);
atmo.o3_number_density = (ozone_vmr_profile_simple(z_grid_km)' * 1e-6) .* ...
    atmo.air_number_density;   % Simplified o3 profile

% Refractive index functions
wavelength_ref_nm = 550;
wavelength_grid_nm = (400:5:780)';
[func_n_r, func_dn_dr] = create_func_n(wavelength_ref_nm, z_grid_km);
sigma_rayleigh_grid = rayleigh_sigma_bodhaine(wavelength_grid_nm, atmo);
sigma_o3_grid = read_chappuis_data(wavelength_grid_nm, atmo.T);

% ASTM Solar spectrum
E_sun = read_astm_e490_data(wavelength_grid_nm);

s_max = 1.0*EARTH_R_km;
b_grid_km = EARTH_R_km + linspace(0, 100, 101)';
s_grid_km = [linspace(0, 500, 11), linspace(600, 1500, 11), linspace(2000, s_max, 7)]';

% Integral to get ray curve
raypath_sol_store = nan(length(b_grid_km), length(s_grid_km), 6);     % [x, y, r, theta, tx, ty]
for i = 1:length(b_grid_km)
    tmp_sol = ode45(@(s,f) ray_curve_ode(s, f, func_n_r, func_dn_dr), [0, s_max], [b_grid_km(i); 0; 0; 1]);
    raypath_sol_store(i, :, [1,2,5,6]) = interp1(tmp_sol.x', tmp_sol.y', s_grid_km, 'pchip');
end
[tmp_theta, tmp_r] = cart2pol(raypath_sol_store(:,:,1), raypath_sol_store(:,:,2));
raypath_sol_store(:,:,3) = tmp_r;
raypath_sol_store(:,:,4) = tmp_theta;
clear tmp_*;

% Find solutions where r = EARTH_SUN_DIST_km
% This requires solving a quadratic equation - a line starting at (x, y) with direction (tx, ty),
% and a circle of radius EARTH_SUN_DIST_km centered at (0,0).
tmp_b = 2 * dot(raypath_sol_store(:, end, 1:2), raypath_sol_store(:, end, 5:6), 3);
tmp_c = squeeze(raypath_sol_store(:, end, 3).^2 - EARTH_SUN_DIST_km.^2);
tmp_s = (-tmp_b + sqrt(tmp_b.^2 - 4*tmp_c)) / 2;
tmp_xy = squeeze(raypath_sol_store(:, end, 1:2) + tmp_s .* raypath_sol_store(:, end, 5:6));
[tmp_theta, tmp_r] = cart2pol(tmp_xy(:,1), tmp_xy(:,2));
sun_sol = [tmp_xy, tmp_r, tmp_theta];   % [x, y, r, theta]_sun
clear tmp_*;

% Emit -> Exit angle from sun
tmp_ray_angle = squeeze(atan2d(raypath_sol_store(:,end,6), raypath_sol_store(:,end,5)));
tmp_sun_theta = atan2d(sun_sol(:,2), sun_sol(:,1));
sun_emit_arcmin = (tmp_ray_angle - tmp_sun_theta)*60;
sun_exit_arcmin = (180 - tmp_ray_angle - tmp_sun_theta)*60;

% Geometric Jacobbian
% TODO: should do some smooth?
emit_fun = pchip(b_grid_km, sun_emit_arcmin);
exit_fun = pchip(b_grid_km, sun_exit_arcmin);
geo_jacob = fnval(fnder(exit_fun, 1), b_grid_km) ./ fnval(fnder(emit_fun, 1), b_grid_km);
clear *_fun;

% Find solutions where r = EARTH_MOON_DIST_km
tmp_b = 2 * dot(raypath_sol_store(:, end, 1:2), raypath_sol_store(:, end, 5:6), 3);
tmp_c = squeeze(raypath_sol_store(:, end, 3).^2 - EARTH_MOON_DIST_km.^2);
tmp_s = (-tmp_b + sqrt(tmp_b.^2 - 4*tmp_c)) / 2;
tmp_xy = squeeze(raypath_sol_store(:, end, 1:2) + tmp_s .* raypath_sol_store(:, end, 5:6));
tmp_xy = [tmp_xy(:,1), -tmp_xy(:,2)];   % Reflect across x-axis to get correct moon position
[tmp_theta, tmp_r] = cart2pol(tmp_xy(:,1), tmp_xy(:,2));
moon_sol = [tmp_xy, tmp_r, tmp_theta];   % [x, y, r, theta]_moon
clear tmp_*;
moon_beta_arcmin = rad2deg(pi + moon_sol(:,4) - sun_sol(:,4)) * 60;

% \kappa_{sct} = N * \sigma_{sct}
raypath_kappa_sct_air = nan(length(b_grid_km), length(s_grid_km), length(wavelength_grid_nm));
raypath_kappa_abs_o3 = nan(length(b_grid_km), length(s_grid_km), length(wavelength_grid_nm));
for i = 1:length(wavelength_grid_nm)
    raypath_kappa_sct_air(:,:,i) = interp1(z_grid_km, atmo.air_number_density .* sigma_rayleigh_grid(i,:), ...
        raypath_sol_store(:,:,3)-EARTH_R_km, 'linear', nan);
    raypath_kappa_abs_o3(:,:,i) = interp1(z_grid_km, atmo.o3_number_density .* sigma_o3_grid(i,:), ...
        raypath_sol_store(:,:,3)-EARTH_R_km, 'linear', nan);
end
raypath_kappa_sct_air(isnan(raypath_kappa_sct_air)) = 0;
raypath_kappa_abs_o3(isnan(raypath_kappa_abs_o3)) = 0;

% optical depth \tau = \int_{s}{\kappa ds}
raypath_tau = squeeze(2 * trapz(s_grid_km*1e3, raypath_kappa_sct_air + raypath_kappa_abs_o3, 2))';

%% Sample rays across sun disk and moon radius
moon_screen_grid_arcmin = linspace(0, 80, 241)';
moon_tan_y = tand(moon_screen_grid_arcmin/60);
sun_radius_grid_arcmin = linspace(0, atand(SUN_R_km / EARTH_SUN_DIST_km) * 60, 91)';
sun_sin_r = sind(sun_radius_grid_arcmin/60);
sun_cos_r = cosd(sun_radius_grid_arcmin/60);

total_sample_num = length(moon_screen_grid_arcmin) * length(sun_radius_grid_arcmin)^2;
moon_data_store = nan(total_sample_num, 4);  % [dist_km, beta_arcmin, idx, w]
sample_count = 0;
for i = 1:length(moon_screen_grid_arcmin)
    for j = 1:length(sun_radius_grid_arcmin)
        tmp_sample_num = 2*j-1;
        sun_phi_grid_rad = linspace(-pi/2, pi/2, tmp_sample_num)';   % From -90 to +90 deg
        sun_sin_phi = sin(sun_phi_grid_rad);
        sun_cos_phi = cos(sun_phi_grid_rad);
        tmp_one = ones(tmp_sample_num, 1);
        tmp_w = tmp_one; tmp_w(1) = 0.5; tmp_w(end) = 0.5;

        dist_moon = EARTH_MOON_DIST_km * ...
            (sun_cos_r(j) - sun_sin_r(j) .* sun_sin_phi .* moon_tan_y(i));  % N x 1

        z_axis = [-sun_cos_phi .* sun_sin_r(j), ...
            -sun_sin_phi .* sun_sin_r(j), ...
            sun_cos_r(j) * tmp_one];   % N x 3
        y_axis = [sun_cos_phi .* sun_sin_r(j) .* (sun_cos_r(j) - sun_sin_phi .* sun_sin_r(j) .* moon_tan_y(i)), ...
            moon_tan_y(i) + sun_sin_phi .* sun_sin_r(j) .* (sun_cos_r(j) - sun_sin_phi .* sun_sin_r(j) .* moon_tan_y(i)), ...
            1 + sun_cos_r(j) .* (-sun_cos_r(j) + sun_sin_phi .* sun_sin_r(j) .* moon_tan_y(i))];  % N x 3
        y_axis = y_axis ./ vecnorm(y_axis, 2, 2);
        x_axis = cross(y_axis, z_axis, 2);   % N x 3

        y_moon = sum([0, moon_tan_y(i) * EARTH_MOON_DIST_km, EARTH_MOON_DIST_km] .* y_axis, 2);
        y_moon_arcmin = atand(y_moon / EARTH_MOON_DIST_km) * 60;    % N x 1
        moon_data_store(sample_count + (1:tmp_sample_num), :) = ...
            [dist_moon, y_moon_arcmin, i*tmp_one, tmp_w];

        sample_count = sample_count + tmp_sample_num;
    end
end
moon_data_store(isnan(moon_data_store(:,2)),:) = [EARTH_MOON_DIST_km, 0, 1, 0.5];

clear tmp_* sun_sin* sun_cos* sun_phi* *_axis y_moon* dist_moon moon_tan_y i j sample_count;

%%
spec_store = nan(length(moon_screen_grid_arcmin), length(wavelength_grid_nm));

geo_disk_val = tand(moon_data_store(:,2)/60) .* moon_data_store(:,1) ./ ...
    (EARTH_SUN_DIST_km+EARTH_MOON_DIST_km) > (EARTH_R_km/EARTH_SUN_DIST_km);
geo_disk_store = nan(length(moon_screen_grid_arcmin), 1);
moon_pos_idx = moon_data_store(:,3);
for i = 1:length(moon_screen_grid_arcmin)
    tmp_idx = moon_pos_idx == i;
    geo_disk_store(i) = sum(geo_disk_val(tmp_idx));
end

% Positive beta
moon_pos_idx = moon_data_store(:,3);
b_interp = interp1(moon_beta_arcmin, b_grid_km, moon_data_store(:,2), 'linear', nan);
jacob_interp = interp1(b_grid_km, geo_jacob, b_interp, 'linear', nan);
transit_factor_interp = exp(-interp1(b_grid_km, raypath_tau', b_interp, 'linear', nan));

jacob_interp(isnan(jacob_interp)) = 1;
transit_factor_interp(isnan(transit_factor_interp)) = 1;

spec_interp = transit_factor_interp ./ jacob_interp .* moon_data_store(:,4);
for i = 1:length(moon_screen_grid_arcmin)
    tmp_idx = moon_pos_idx == i;
    tmp_spec = spec_interp(tmp_idx, :);
    spec_store(i, :) = sum(tmp_spec);
end

[~, tmp_moon_beta_idx] = sort(moon_data_store(:,2));
line_colors = colormap('lines');
figure(1); clf; hold on;
scatter(moon_screen_grid_arcmin(moon_pos_idx(1:80:end)), 1./jacob_interp(1:80:end), 8, 'o', 'filled', ...
    'MarkerEdgeColor', line_colors(1,:), 'MarkerFaceColor', line_colors(1,:), ...
    'MarkerEdgeAlpha', 0.10, 'MarkerFaceAlpha', 0.05);
scatter(moon_screen_grid_arcmin(1:10:end), spec_store(1:10:end, 75), 28, 'o', 'filled', ...
    'MarkerEdgeColor', line_colors(3,:), 'MarkerFaceColor', 'none', ...
    'LineWidth', 1.5);
scatter(moon_screen_grid_arcmin(1:10:end), spec_store(1:10:end, 10), 28, '*', 'filled', ...
    'MarkerEdgeColor', line_colors(3,:), 'MarkerFaceColor', 'none', ...
    'LineWidth', 1.5);
h_line = plot(moon_data_store(tmp_moon_beta_idx,2), 1./jacob_interp(tmp_moon_beta_idx), '-', ...
    moon_screen_grid_arcmin, spec_store(:, 75), '-', ...
    moon_screen_grid_arcmin, spec_store(:, 10), '-', ...
    moon_screen_grid_arcmin, geo_disk_store, '--');
h_line(1).LineWidth = 2.0; h_line(1).Color = line_colors(2,:);
h_line(2).LineWidth = 1.5; h_line(2).MarkerSize = 4; h_line(2).Color = line_colors(3,:);
h_line(3).LineWidth = 1.5; h_line(3).MarkerSize = 4; h_line(3).Color = line_colors(3,:);
h_line(4).LineWidth = 1.5; h_line(4).Color = line_colors(4,:);
box on;
set(gca, 'YScale', 'log','YGrid', 'on', 'XGrid', 'on',...
    'xlim', [0, max(moon_screen_grid_arcmin)]);

% Negative beta
b_interp = interp1(moon_beta_arcmin, b_grid_km, -moon_data_store(:,2), 'linear', nan);
tmp_idx = ~isnan(b_interp);
b_interp = b_interp(tmp_idx);
moon_pos_idx = moon_data_store(tmp_idx,3);
jacob_interp = interp1(b_grid_km, geo_jacob, b_interp, 'linear', nan);
transit_factor_interp = exp(-interp1(b_grid_km, raypath_tau', b_interp, 'linear', nan));

jacob_interp(isnan(jacob_interp)) = 1;
transit_factor_interp(isnan(transit_factor_interp)) = 1;

spec_interp = transit_factor_interp ./ jacob_interp .* moon_data_store(tmp_idx,4);
for i = 1:length(moon_screen_grid_arcmin)
    tmp_idx = moon_pos_idx == i;
    tmp_spec = spec_interp(tmp_idx, :);
    spec_store(i, :) = spec_store(i, :) + sum(tmp_spec);
end

clear tmp_* *_interp i moon_pos_idx;

%%
figure(3); clf;
subplot(2,1,1);
rgb_color = spec_to_rgb([wavelength_grid_nm, spec_store' .* E_sun], 'Mixed', true, 'Y', 2^18);
imagesc(moon_screen_grid_arcmin, 1, reshape(rgb_color, 1, [], 3));
subplot(2,1,2);
rgb_color = spec_to_rgb([wavelength_grid_nm, spec_store' .* E_sun], 'Mixed', true, 'Y', 2^0);
imagesc(moon_screen_grid_arcmin, 1, reshape(rgb_color, 1, [], 3));
