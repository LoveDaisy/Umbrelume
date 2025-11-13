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
atmo.o3_number_density = (ozone_vmr_profile_simple(z_grid_km)' * 1e-6) .* atmo.air_number_density; % Simplified o3 profile

% Refractive index functions
nr_func = create_nr_func(mean(wavelength_grid_nm), z_grid_km);

% Curve arc-length parameter grid
s_grid_km = [linspace(0, 500, 11), linspace(600, 1500, 11), linspace(2000, 1.0 * Constants.EARTH_R_km, 7)]';

% Integrate and get the raypath
raypath_store = compute_raypath(b_grid_km, s_grid_km, nr_func);

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

function raypath_store = compute_raypath(b_grid_km, s_grid_km, nr_func)
% COMPUTE_RAYPATH Compute the raypath for given b and s values.
%   Inputs:
%       b_grid_km:          [B x 1] array of b values in kilometers
%       s_grid_km:          [S x 1] array of s values in kilometers
%       nr_func:            A struct contains two functions:
%           - n_r:          Function handle of n(r)
%           - dn_dr:        Function handle of dn/dr(r)
%   Outputs:
%       raypath_store:      [B x S x 6] array of raypath data

% Initialize raypath store
raypath_store = nan(length(b_grid_km), length(s_grid_km), 6);

% Compute raypath for each b and s
for i = 1:length(b_grid_km)
    % Solve ray curve ODE
    tmp_sol = ode45(@(s,f) ray_curve_ode(s, f, nr_func.n_r, nr_func.dn_dr), [0, max(s_grid_km)], [b_grid_km(i); 0; 0; 1]);
    raypath_store(i, :, [1,2,5,6]) = interp1(tmp_sol.x', tmp_sol.y', s_grid_km, 'pchip');
end
[tmp_theta, tmp_r] = cart2pol(raypath_store(:,:,1), raypath_store(:,:,2));
raypath_store(:,:,3) = tmp_r;        % r
raypath_store(:,:,4) = tmp_theta;    % theta
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