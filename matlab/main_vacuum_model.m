clear; close all; clc;

spec_wl_nm = (400:5:780)';
moon_y_arcmin = linspace(0, 80, 241)';

% ASTM Solar spectrum
E_sun = read_astm_e490_data(spec_wl_nm);

SUN_R_arcmin = atand(Constants.SUN_R_km/Constants.EARTH_SUN_DIST_km)*60;
[sun_r_theta_arcmin, w] = generate_sun_half_disk_sample(SUN_R_arcmin, 91);

N = size(sun_r_theta_arcmin, 1);
M = length(moon_y_arcmin);

cq = cosd(sun_r_theta_arcmin(:,2)/60);
sq = sind(sun_r_theta_arcmin(:,2)/60);
cr = cosd(sun_r_theta_arcmin(:,1)/60);
sr = sind(sun_r_theta_arcmin(:,1)/60);
ty = tand(moon_y_arcmin/60);

sun_xyz = [cq.*sr, sq.*sr, -cr] * Constants.EARTH_SUN_DIST_km;  % [N x 3]
moon_xyz = Constants.EARTH_MOON_DIST_km * [zeros(M, 1), ty, ones(M, 1)]; % [M x 3]

d_vec = cat(3, sun_xyz(:,1) - moon_xyz(:,1)', sun_xyz(:,2) - moon_xyz(:,2)', ...
    sun_xyz(:,3) - moon_xyz(:,3)');     % [N x M x 3]
d_vec2 = sum(d_vec .* d_vec, 3);        % [N x M]
s = -sum(reshape(moon_xyz, [1, M, 3]) .* d_vec, 3) ./ d_vec2;
r = d_vec .* s + reshape(moon_xyz, [1, M, 3]);
r2 = sum(r .* r, 3);
lit_factor = sum(r2 > Constants.EARTH_R_km^2, 1);

spec = lit_factor .* E_sun;

%%
figure(1); clf;
subplot(2,1,1);
rgb_color = spec_to_rgb([spec_wl_nm, spec], 'Mixed', true, 'Y', 2^18);
imagesc(moon_y_arcmin, 1, reshape(rgb_color, 1, [], 3));
subplot(2,1,2);
rgb_color = spec_to_rgb([spec_wl_nm, spec], 'Mixed', true, 'Y', 2^0);
imagesc(moon_y_arcmin, 1, reshape(rgb_color, 1, [], 3));
