function T = direct_term(wavelength_nm, atmo, r_arcmin)
% length(wavelength_nm) = m
% length(atmo.z) = n
% length(r_arcmin) = k

% Constants
RE_km = 6371.0;
D_SE_km = 1.496e8;
D_EM_km = 38.4e4;

r_rad = deg2rad(r_arcmin(:) / 60)';     % [1xk]
z = atmo.z(:)';                         % [1xn]

% Geometry
u_S = asin(RE_km / D_SE_km);
u_M = asin(RE_km / D_EM_km .* cos(r_rad)) - r_rad;
alpha_r = max(u_S + u_M, 0);                        % bending angle at radius r. [1xk]
alpha_z = bending_angle(wavelength_nm, atmo);       % bending angle along z [mxn]

% Inverse mapping from bending angle to z.
% We use middle wavelength as default.
mid_i = floor(length(wavelength_nm) / 2);
[tmp_alpha, tmp_idx] = unique(alpha_z(mid_i, :));
tmp_z = z(tmp_idx);
z_star = interp1(tmp_alpha, tmp_z, alpha_r, 'linear', 'extrap');    % [mxk]

% Optical depth along z_star.
[tau_air, tau_o3] = optical_depth_along_raypath(wavelength_nm, atmo, z_star);   % [mxk]
% Convert to transmitance
% T_air = exp(-tau_air);      % [mxk]
% T_o3 = exp(-tau_o3);        % [mxk]
T = exp(-(tau_air + tau_o3)); % [mxk]
end