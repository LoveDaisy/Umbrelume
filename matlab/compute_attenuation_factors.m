function [transit_factor, geo_factor] = compute_attenuation_factors(sun_r_theta_arcmin, moon_y_arcmin, grid_data)
% COMPUTE_ATTENUATION_FACTORS Compute the attenuation factors for given
% sun positions and moon positions at specified wavelengths using precomputed
% grid data.
%   Inputs:
%       sun_r_theta_arcmin:         [N x 2] array of sun positions in arcminutes
%       moon_y_arcmin:              [M x 1] array of moon positions in arcminutes
%       grid_data:                  Struct containing b grid data with fields:
%          - b_grid_km:                 [B x 1] array of b values in kilometers
%          - wavelength_grid_nm:        [L x 1] array of wavelength grid in nanometers
%          - geo_jacob_grid:            [B x 1] array of geometric Jacobian values
%          - beta_arcmin_grid:          [B x 1] array of beta angles in arcminutes
%          - tau_grid:                  [B x L] array of optical depths
%   Outputs:
%       transit_factor:     [N x M x L] array of transit attenuation factors
%       geo_factor:         [N x M] array of geometric attenuation factors

moon_y_arcmin = moon_y_arcmin(:);

% We need to normalize the moon position into an earth-centered local coordinate system.
%   - new z-axis: from sun sample point to earth center;
%   - new y-axis: keep the moon point on y;
%   - new x-axis: cross product of new y- and z- axis;
% But we DO NOT need all coordinates. We need only:
%   1) beta angle between earth-moon line and new z-axis (beta)

% Compute beta using the new y- and z- axis.
cq = cosd(sun_r_theta_arcmin(:,2)/60);
sq = sind(sun_r_theta_arcmin(:,2)/60);
cr = cosd(sun_r_theta_arcmin(:,1)/60);
sr = sind(sun_r_theta_arcmin(:,1)/60);
ty = tand(moon_y_arcmin/60);

z_moon_new = cr - sq .* sr .* ty';  % [N x M]
y_vec_new = cat(3, cq.*sr.*(cr - sq.*sr.*ty'), ...
    cr.*sq.*sr+(cr.^2+cq.^2.*sr.^2).*ty', ...
    sr.*(sr + cr.*sq.*ty'));        % [N x M x 3]
y_moon_new = sqrt(sum(y_vec_new.^2, 3));    % [N x M]
beta_arcmin = atan2d(y_moon_new, z_moon_new)*60;  % [N x M]

% Interpolate grid data to get tau and geo_jacob at required beta and wavelength
[transit_factor_pos, geo_factor_pos] = interpolate_data(beta_arcmin, grid_data);
[transit_factor_neg, geo_factor_neg] = interpolate_data(-beta_arcmin, grid_data);

transit_factor = transit_factor_pos + transit_factor_neg;
geo_factor = geo_factor_pos + geo_factor_neg;
end

function [transit_factor, geo_factor] = interpolate_data(beta_arcmin, grid_data)
% INTERPOLATE_DATA Interpolate the precomputed grid data to get tau and
% geometric Jacobian at specified beta angles and wavelengths.
%   Inputs:
%       beta_arcmin:            [N x M] array of beta angles in arcminutes
%       grid_data:              Struct containing b grid data with fields:
%          - b_grid_km:                 [B x 1] array of b values in kilometers
%          - wavelength_grid_nm:        [L x 1] array of wavelength grid in nanometers
%          - geo_jacob_grid:            [B x 1] array of geometric Jacobian values
%          - beta_arcmin_grid:          [B x 1] array of beta angles in arcminutes
%          - tau_grid:                  [B x L] array of optical depths
%   Outputs:
%       transit_factor:         [N x M x L] array of transit attenuation factors
%       geo_factor:             [N x M] array of geometric attenuation factors

[N, M] = size(beta_arcmin);
L = length(grid_data.wavelength_grid_nm);

beta_arcmin = beta_arcmin(:);   % [N*M x 1]
valid_idx = beta_arcmin >= min(grid_data.beta_arcmin_grid) & ...
    beta_arcmin <= max(grid_data.beta_arcmin_grid);
pos_idx = beta_arcmin > max(grid_data.beta_arcmin_grid);
beta_arcmin_valid = beta_arcmin(valid_idx);  % [valid x 1]

transit_factor = zeros(N*M, L); % [N*M x L]
tau_interp = interp1(grid_data.beta_arcmin_grid, grid_data.tau_grid, beta_arcmin_valid, 'linear', nan);   % [valid x L]
transit_factor(valid_idx, :) = exp(-tau_interp);  % [valid x L]
transit_factor(pos_idx, :) = 1;  % Outside grid, set transit_factor = 1
transit_factor = reshape(transit_factor, [N, M, L]);  % [N x M x L]

geo_factor = zeros(N*M, 1); % [N*M x 1]
geo_jacob_interp = interp1(grid_data.beta_arcmin_grid, grid_data.geo_jacob_grid, beta_arcmin_valid, 'linear', nan);   % [valid x 1]
geo_factor(valid_idx) = 1./geo_jacob_interp;  % [valid x 1]
geo_factor(pos_idx) = 1;  % Outside grid, set geo_factor = 1
geo_factor = reshape(geo_factor, [N, M]);  % [N x M]
end