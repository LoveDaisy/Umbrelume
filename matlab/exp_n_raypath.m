clear; close all; clc;

b_km = Constants.EARTH_R_km + linspace(0, 100, 101)';
spec_wl_nm = (400:5:780)';

% Atmosphere model
z_km = (0:.1:(max(b_km) - Constants.EARTH_R_km))';
atmo = us_standard_atmosphere_1976(z_km);
atmo.o3_number_density = (ozone_vmr_profile_simple(z_km)' * 1.1 * 1e-6) .* atmo.air_number_density; % Simplified o3 profile

% Curve arc-length parameter grid
s_km = [linspace(0, 500, 11), linspace(600, 1500, 11), linspace(2000, 1.0 * Constants.EARTH_R_km, 7)]';    % [S x 1]

% Integrate and get the raypath
nr_func = create_nr_func(spec_wl_nm(1), atmo);
raypath_store_1 = compute_raypath(b_km, s_km, nr_func); % [B x S x 6]

nr_func = create_nr_func(spec_wl_nm(end), atmo);
raypath_store_2 = compute_raypath(b_km, s_km, nr_func); % [B x S x 6]
