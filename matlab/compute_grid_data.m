function grid_data = compute_grid_data(wavelength_grid_nm, b_grid_km)
% COMPUTE_GRID_DATA Precompute the grid data for attenuation factor computation.
%   Inputs:
%       wavelength_grid_nm:        [L x 1] array of wavelength grid in nanometers
%       b_grid_km:                 [B x 1] array of b values in kilometers
%   Outputs:
%       grid_data:                 Struct containing b grid data with fields:
%          - b_grid_km:                 [B x 1] array of b values in kilometers
%          - wavelength_grid_nm:        [L x 1] array of wavelength grid in nanometers
%          - geo_jacob_grid:            [B x 1] array of geometric Jacobian values
%          - beta_grid_arcmin:          [B x 1] array of beta angles in arcminutes
%          - tau_grid:                  [B x K] array of optical depths

% Make sure inputs are column vectors
b_grid_km = b_grid_km(:);       % [B x 1]
wavelength_grid_nm = wavelength_grid_nm(:);  % [L x 1]

% Atmosphere model
atmo = generate_atmosphere_model();

% Curve arc-length parameter grid
s_grid_km = [linspace(0, 500, 11), linspace(600, 1500, 11), [1800, 2300, 2800]]';    % [S x 1]

% Integrate and get the raypath
nr_func = create_nr_func(mean(wavelength_grid_nm), atmo);
raypath_store = compute_raypath(nr_func, b_grid_km, s_grid_km); % [B x S x 6]

% Solve raypath at r = EARTH_SUN_DIST_km and r = EARTH_MOON_DIST_km
sun_sol = solve_raypath(raypath_store, Constants.EARTH_SUN_DIST_km);   % [x, y, r, theta] at sun position
moon_sol = solve_raypath(raypath_store, Constants.EARTH_MOON_DIST_km); % [x, y, r, theta] at moon position
sun_sol(:,1) = -sun_sol(:,1);
sun_sol(:,4) = mod(pi - sun_sol(:,4), 2*pi);

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
    raypath_kappa_sct_air(:,:,i) = interp1(atmo.z, atmo.air_number_density .* sigma_rayleigh_grid(i,:), ...
        raypath_store(:,:,3)-Constants.EARTH_R_km, 'linear', nan);
    raypath_kappa_abs_o3(:,:,i) = interp1(atmo.z, atmo.o3_number_density .* sigma_o3_grid(i,:), ...
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

function geo_jacob = compute_geometric_jacobian(b_grid_km, raypath_store, sun_sol)
% COMPUTE_GEOMETRIC_JACOBIAN Compute the geometric Jacobian for given b values.
%   Inputs:
%       b_grid_km:          [B x 1] array of b values in kilometers
%       raypath_store:      [B x S x 6] array of raypath data
%       sun_sol:            [B x 4] array of [x, y, r, theta] at sun position
%   Outputs:
%       geo_jacob:          [B x 1] array of geometric Jacobian values

% Emit -> Exit angle from sun
ray_angle = squeeze(atan2(raypath_store(:,end,6), raypath_store(:,end,5)));
sun_theta = sun_sol(:,4);
sun_emit = pi - ray_angle - sun_theta;
sun_exit = pi - ray_angle + sun_theta;

% Geometric Jacobian
% TODO: should do some smooth?
emit_fun = pchip(b_grid_km, sun_emit);
exit_fun = pchip(b_grid_km, sun_exit);
geo_jacob = abs(fnval(fnder(exit_fun, 1), b_grid_km) ./ fnval(fnder(emit_fun, 1), b_grid_km));  % [B x 1]
end