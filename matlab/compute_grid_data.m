function grid_data = compute_grid_data(b_grid_km, wavelength_grid_nm)
% COMPUTE_GRID_DATA Precompute the grid data for attenuation factor computation.
%   Inputs:
%       b_grid_km:                 [B x 1] array of b values in kilometers
%       wavelength_grid_nm:        [K x 1] array of wavelength grid in nanometers
%   Outputs:
%       grid_data:                 Struct containing b grid data with fields:
%          - b_grid_km:                 [B x 1] array of b values in kilometers
%          - wavelength_grid_nm:        [K x 1] array of wavelength grid in nanometers
%          - geo_jacob_grid:            [B x 1] array of geometric Jacobian values
%          - beta_grid_arcmin:          [B x 1] array of beta angles in arcminutes
%          - tau_grid:                  [B x K] array of optical depths

% Make sure inputs are column vectors
b_grid_km = b_grid_km(:);       % [B x 1]
wavelength_grid_nm = wavelength_grid_nm(:);  % [K x 1]

% Atmosphere model
z_max = max(b_grid_km) - Constants.EARTH_R_km;
z_grid_km = (0:.1:z_max)';
atmo = us_standard_atmosphere_1976(z_grid_km);
atmo.o3_number_density = (ozone_vmr_profile_simple(z_grid_km)' * 1.1 * 1e-6) .* atmo.air_number_density; % Simplified o3 profile

% Curve arc-length parameter grid
s_grid_km = [linspace(0, 500, 11), linspace(600, 1500, 11), linspace(2000, 1.0 * Constants.EARTH_R_km, 7)]';    % [S x 1]

% Integrate and get the raypath
nr_func = create_nr_func(mean(wavelength_grid_nm), atmo);
raypath_store = compute_raypath(b_grid_km, s_grid_km, nr_func); % [B x S x 6]

% Solve raypath at r = EARTH_SUN_DIST_km and r = EARTH_MOON_DIST_km
sun_sol = solve_raypath(raypath_store, Constants.EARTH_SUN_DIST_km);   % [x, y, r, theta] at sun position
moon_sol = solve_raypath(raypath_store, Constants.EARTH_MOON_DIST_km); % [x, y, r, theta] at moon position
moon_sol(:,2) = -moon_sol(:,2);
moon_sol(:,4) = -moon_sol(:,4);

% Compute geometric Jacobian
geo_jacob_grid = compute_geometric_jacobian(b_grid_km, raypath_store, sun_sol); % [B x 1]

% Compute moon beta angles
beta_grid_arcmin = rad2deg(pi + moon_sol(:,4) - sun_sol(:,4)) * 60;     % [B x 1]

% Compute optical depth tau
sigma_rayleigh_grid = rayleigh_sigma_bodhaine(wavelength_grid_nm, atmo);
sigma_o3_grid = read_chappuis_data(wavelength_grid_nm, atmo.T);

% \kappa_{sct} = N * \sigma_{sct}
raypath_kappa_sct_air = nan(length(b_grid_km), length(s_grid_km), length(wavelength_grid_nm));
raypath_kappa_abs_o3 = nan(length(b_grid_km), length(s_grid_km), length(wavelength_grid_nm));
for i = 1:length(wavelength_grid_nm)
    raypath_kappa_sct_air(:,:,i) = interp1(z_grid_km, atmo.air_number_density .* sigma_rayleigh_grid(i,:), ...
        raypath_store(:,:,3)-Constants.EARTH_R_km, 'linear', nan);
    raypath_kappa_abs_o3(:,:,i) = interp1(z_grid_km, atmo.o3_number_density .* sigma_o3_grid(i,:), ...
        raypath_store(:,:,3)-Constants.EARTH_R_km, 'linear', nan);
end
raypath_kappa_sct_air(isnan(raypath_kappa_sct_air)) = 0;
raypath_kappa_abs_o3(isnan(raypath_kappa_abs_o3)) = 0;

% optical depth \tau = \int_{s}{\kappa ds}
tau_grid = squeeze(2 * trapz(s_grid_km*1e3, raypath_kappa_sct_air + raypath_kappa_abs_o3, 2));  % [B x K]

% Construct output grid data
grid_data.b_grid_km = b_grid_km;
grid_data.wavelength_grid_nm = wavelength_grid_nm;
grid_data.geo_jacob_grid = geo_jacob_grid;
grid_data.beta_arcmin_grid = beta_grid_arcmin;
grid_data.tau_grid = tau_grid;
end

function sol = solve_raypath(raypath_store, target_r)
% SOLVE_RAYPATH Solve the raypath to find the position where r = target_r
%   Inputs:
%       raypath_store:      [S x 6] array of raypath data along s
%       target_r:           [1 x 1] target r value
%   Outputs:
%       sol:                [S x 4] array of [x, y, r, theta] at target r

% Basically we need to find the intersection of a circle with radius target_r,
% and the tangent line at the end of the raypath.
coef_b = 2 * dot(raypath_store(:, end, 1:2), raypath_store(:, end, 5:6), 3);
coef_c = squeeze(raypath_store(:, end, 3).^2 - target_r.^2);
sol_s = (-coef_b + sqrt(coef_b.^2 - 4*coef_c)) / 2;
sol_xy = squeeze(raypath_store(:, end, 1:2) + sol_s .* raypath_store(:, end, 5:6));
[sol_theta, sol_r] = cart2pol(sol_xy(:,1), sol_xy(:,2));
sol = [sol_xy, sol_r, sol_theta];   % [x, y, r, theta]
end

function geo_jacob = compute_geometric_jacobian(b_grid_km, raypath_store, sun_sol)
% COMPUTE_GEOMETRIC_JACOBIAN Compute the geometric Jacobian for given b values.
%   Inputs:
%       b_grid_km:          [B x 1] array of b values in kilometers
%       raypath_store:      [B x S x 6] array of raypath data
%       sun_sol:            [B x 4] array of [x, y, r, theta] at sun position
%   Outputs:
%       geo_jacob:          [B x 1] array of geometric Jacobian values

% Emit -> Exit angle from sun
ray_angle = squeeze(atan2d(raypath_store(:,end,6), raypath_store(:,end,5)));
sun_theta = atan2d(sun_sol(:,2), sun_sol(:,1));
sun_emit_arcmin = (ray_angle - sun_theta)*60;
sun_exit_arcmin = (180 - ray_angle - sun_theta)*60;

% Geometric Jacobian
% TODO: should do some smooth?
emit_fun = pchip(b_grid_km, sun_emit_arcmin);
exit_fun = pchip(b_grid_km, sun_exit_arcmin);
geo_jacob = fnval(fnder(exit_fun, 1), b_grid_km) ./ fnval(fnder(emit_fun, 1), b_grid_km);  % [B x 1]
end