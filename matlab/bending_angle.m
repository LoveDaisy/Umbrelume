function alpha = bending_angle(wavelength_nm, atmo)
% Constants
RE_KM = 6371.0;

z = atmo.z(:)';                                   % should be a row vector [1xn]
T = atmo.T(:)';
p = atmo.p(:)';
lam_nm = wavelength_nm(:);                        % a column vector [mx1]

% 参考弯折角积分式: α(a) = -2 a ∫ (d ln n/dr) / sqrt((n r)^2 - a^2) dr
% 文献: Fjeldbo-Abel框架 (见 Syndergaard 2005 教程) 与 Auer–Standish 非奇异积分方法.
r = RE_KM + z;                          % [1xn]
r_mid = (r(1:end-1) + r(2:end)) / 2;    % [1x(n-1)]
T_mid = (T(1:end-1) + T(2:end)) / 2;
p_mid = (p(1:end-1) + p(2:end)) / 2;
n = ciddor_n_air(lam_nm, T, p);                 % [mxn]
n_mid = ciddor_n_air(lam_nm, T_mid, p_mid);     % [mx(n-1)]
nr = n .* r;                                    % [mxn]
nr_mid = n_mid .* r_mid;                        % [mx(n-1)]

% 计算数值导数: d ln n / dr = (1/n) * dn/dr
dn_dr = pchip(r_mid, n_mid);        % [mx(n-1)]
dn_dr = fnder(dn_dr, 1);
dn_dr = fnval(dn_dr, r);
dlnn_dr = dn_dr ./ n;               % [mx(n-1)]

alpha = zeros(size(nr_mid));        % [mx(n-1)]
for i = 1:length(r_mid)-1
    a = nr(:, i);       % [mx1]
    f = dlnn_dr(:, i) ./ max(sqrt(nr_mid(:, i:end).^2 - a.^2), 1e-12);
    alpha(:, i) = trapz(r_mid(i:end), f')' .* (-2 * a);
end
alpha = interp1(r_mid, alpha', r, 'linear', 'extrap')';

% 单调性修正
for i = 2:length(r)
    alpha(:, i) = min(alpha(:, i-1), alpha(:, i));
end
end
