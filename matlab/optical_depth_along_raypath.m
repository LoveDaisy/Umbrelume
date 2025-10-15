function [tau_air, tau_o3] = optical_depth_along_raypath(wavelength_nm, atmo, z0)
% 在每个薄层用折射不变量 n(r) * r * sinψ(r) 计算 cosψ(r),
% 离散近似 ds = Δr / cosψ(r). 两侧对称最后需乘以 2

% Constants
RE_KM = 6371.0;

lam_nm = wavelength_nm(:);      % [mx1]
[z, idx] = unique(atmo.z(:)');  % [1xn]
T = atmo.T(idx);                % [1xn]
p = atmo.p(idx);                % [1xn]
z0 = max(z0(:)', 0.05);         % [1xn1]

r = RE_KM + z;                          % [1xn]
dr = diff(r);
r_mid = (r(1:end-1) + r(2:end)) / 2;    % [1x(n-1)]
T_mid = (T(1:end-1) + T(2:end)) / 2;    % [1x(n-1)]
p_mid = (p(1:end-1) + p(2:end)) / 2;    % [1x(n-1)]

n_mid = ciddor_n_air(lam_nm, T_mid, p_mid); % [mx(n-1)]
nr_mid = n_mid .* r_mid;           % [mx(n-1)]

% interpolate density to midpoints
N_air_mid = interp1(r, atmo.air_number_density(idx), r_mid, 'linear', 'extrap');     % [1x(n-1)]
N_o3_mid = interp1(r, atmo.o3_number_density(idx), r_mid, 'linear', 'extrap');       % [1x(n-1)]

% Cross section
sigma_air = rayleigh_sigma_bodhaine(wavelength_nm, atmo);
sigma_air = sigma_air(:, idx);
sigma_o3 = read_chappuis_data(wavelength_nm, atmo.T);
sigma_o3 = sigma_o3(:, idx);
sigma_air_mid = (sigma_air(:,1:end-1) + sigma_air(:,2:end)) / 2;    % [mx(n-1)]
sigma_o3_mid = (sigma_o3(:,1:end-1) + sigma_o3(:,2:end)) / 2;       % [mx(n-1)]

% ds along the raypath at z_mid
a = nr_mid(:, 1);                     % [mx1], 折射不变量
sin_phi = a ./ nr_mid;                % [mx(n-1)]
cos_phi = sqrt(max(1 - sin_phi.^2, 1e-12));  % [mx(n-1)]
ds_mid = dr ./ cos_phi * 1e3;         % [mx(n-1)], km -> m

tau_air = zeros(length(lam_nm), length(z0));    % [mxn1]
tau_o3 = zeros(length(lam_nm), length(z0));     % [mxn1]
arc_length = zeros(length(lam_nm), length(z0));
for i = 1:length(z0)
    curr_idx = r_mid > (z0(i) + RE_KM);

    % raypath at one side
    curr_ds = ds_mid(:, curr_idx);        % [mxn2]
    arc_length(:,i) = sum(curr_ds, 2);

    % path integration. [1x(n-1)] .* [mx(n-1)] .* [mx(n-1)] -> [mx1]
    tau_air(:,i) = sum(N_air_mid(curr_idx) .* sigma_air_mid(:,curr_idx) .* curr_ds, 2) * 2;     % [mx1]
    tau_o3(:,i) = sum(N_o3_mid(curr_idx) .* sigma_o3_mid(:,curr_idx) .* curr_ds, 2) * 2;        % [mx1]
end
end
