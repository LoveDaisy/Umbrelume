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
    tmp_sol = ode45(@(s,f) ray_curve_ode(s, f, nr_func.n_r, nr_func.dn_dr), [0, max(s_grid_km)], [0; b_grid_km(i); 1; 0]);
    raypath_store(i, :, [1,2,5,6]) = interp1(tmp_sol.x', tmp_sol.y', s_grid_km, 'pchip');
end
[tmp_theta, tmp_r] = cart2pol(raypath_store(:,:,1), raypath_store(:,:,2));
raypath_store(:,:,3) = tmp_r;        % r
raypath_store(:,:,4) = tmp_theta;    % theta
end