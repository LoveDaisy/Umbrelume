clear; close all; clc;

b_km = Constants.EARTH_R_km + linspace(0, 100, 101)';
spec_wl_nm = (400:5:780)';
moon_y_arcmin = linspace(0, 80, 241)';

SUN_R_arcmin = atand(Constants.SUN_R_km/Constants.EARTH_SUN_DIST_km)*60;

% ASTM Solar spectrum
E_sun = read_astm_e490_data(spec_wl_nm);

% Generate sun disk sample points
[sun_r_theta_arcmin, w] = generate_sun_half_disk_sample(SUN_R_arcmin, 91);

% Compute grid data
grid_data = compute_grid_data(spec_wl_nm, b_km);

% Compute attenuation factors
[transit_factor, geo_factor] = compute_attenuation_factors(sun_r_theta_arcmin, moon_y_arcmin, grid_data);

%%
% Get the spectrum
spec = squeeze(sum(transit_factor .* geo_factor .* w, 1))' .* E_sun;

figure(1); clf;
subplot(2,1,1);
rgb_color = spec_to_rgb([spec_wl_nm, spec], 'Mixed', true, 'Y', 2^18);
imagesc(moon_y_arcmin, 1, reshape(rgb_color, 1, [], 3));
subplot(2,1,2);
rgb_color = spec_to_rgb([spec_wl_nm, spec], 'Mixed', true, 'Y', 2^0);
imagesc(moon_y_arcmin, 1, reshape(rgb_color, 1, [], 3));
