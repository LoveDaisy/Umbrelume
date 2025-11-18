function raypath_store = compute_raypath(nr_func, b_grid_km, s_grid_km, varargin)
% COMPUTE_RAYPATH Compute the raypath for given b and s values.
%   Inputs:
%       nr_func:            A struct contains two functions:
%           - n_r:          Function handle of n(r)
%           - dn_dr:        Function handle of dn/dr(r)
%       b_grid_km:          [B x 1] array of b values in kilometers
%       s_grid_km:          [S x 1] array of s values in kilometers
%   Optional Inputs:
%       t_grid_deg:         [T x 1] array of t values in degrees
%   Outputs:
%       raypath_store:      [B x S x 6] array of raypath data, or
%                           [B x S x T x 6] array if t_grid_deg is provided

B = length(b_grid_km);
S = length(s_grid_km);
if ~isempty(varargin) && ~isempty(varargin{1})
    t_grid_deg = varargin{1}(:);
    T = length(t_grid_deg);
else
    t_grid_deg = 0;
    T = 1;
end

% Initialize raypath store
raypath_store = nan(B, S, T, 6);

% Compute raypath for each b and s
for i = 1:B
    for j = 1:T
        % Solve ray curve ODE
        t = [cosd(t_grid_deg(j)); sind(t_grid_deg(j))];
        sol = ode45(@(s,f) ray_curve_ode(s, f, nr_func.n_r, nr_func.dn_dr), [0, max(s_grid_km)], [0; b_grid_km(i); t]);
        sol_idx = all(~isnan(sol.y));
        sol_x = sol.x(sol_idx)';        % [N x 1]
        sol_y = sol.y(:, sol_idx)';     % [N x 4]
        if sum(sol_idx) > 1
            raypath_store(i, :, j, [1,2,5,6]) = interp1(sol_x, sol_y, s_grid_km, 'pchip', nan);
        end
    end
end
[theta, r] = cart2pol(raypath_store(:,:,:,1), raypath_store(:,:,:,2));
raypath_store(:,:,:,3) = r;        % r
raypath_store(:,:,:,4) = theta;    % theta
raypath_store = squeeze(raypath_store);  % [B x S x 6] or [B x S x T x 6]
end