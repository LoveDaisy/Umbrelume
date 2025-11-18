function sol = solve_raypath(raypath_store, target_r)
% SOLVE_RAYPATH Solve the raypath to find the position where r = target_r
%   Inputs:
%       raypath_store:      [B x S x 6] array of raypath data along s, or
%                           [B x S x T x 6] array if multiple t are provided
%       target_r:           [1 x 1] target r value
%   Outputs:
%       sol:                [B x 4] array of [x, y, r, theta] at target r, or
%                           [B x T x 4] array if multiple t are provided

if ndims(raypath_store) == 3
    raypath_store = reshape(raypath_store, size(raypath_store,1), size(raypath_store,2), 1, size(raypath_store,3));
end

% Basically we need to find the intersection of a circle with radius target_r,
% and the tangent line at the end of the raypath.
coef_b = 2 * dot(raypath_store(:, end, :, 1:2), raypath_store(:, end, :, 5:6), 4);
coef_c = raypath_store(:, end, :, 3).^2 - target_r.^2;
sol_s = (-coef_b + sqrt(coef_b.^2 - 4*coef_c)) / 2;
sol_xy = raypath_store(:, end, :, 1:2) + sol_s .* raypath_store(:, end, :, 5:6);    % [B x 1 x T x 2]
[sol_theta, sol_r] = cart2pol(sol_xy(:,:,:,1), sol_xy(:,:,:,2));
sol = cat(4, sol_xy, sol_r, mod(sol_theta, 2*pi));   % [x, y, r, theta]
sol = squeeze(sol);   % [B x 4] or [B x T x 4]
end