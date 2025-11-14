function nr_func = create_nr_func(wavelength_nm, atmo)
% CREATE_NR_FUNC Create refractive index functions n(r) and dn/dr(r)
%   Inputs:
%       wavelength_nm:      [L x 1] array of wavelength grid in nanometers
%       atmo:               Struct containing atmospheric profile with fields:
%           - z:                [N x 1] array of altitude in kilometers
%           - T:                [N x 1] array of temperature in Kelvin
%           - p:                [N x 1] array of pressure in Pascal
%   Outputs:
%       nr_func:            A struct contains two functions:
%           - n_r:          Function handle of n(r)
%           - dn_dr:        Function handle of dn/dr(r)

wavelength_nm = wavelength_nm(:);               % [L x 1]
b_grid_km = Constants.EARTH_R_km + atmo.z(:);   % [N x 1]

L = length(wavelength_nm);
r_range = [min(b_grid_km), max(b_grid_km)];

n = ciddor_n_air(wavelength_nm, atmo.T, atmo.p);    % [L x N]
pp_n_r = pchip(b_grid_km, n);   % A function will return an array of shape [L x length(r)]
func_dn_dr = @(r) helper_dn(pp_n_r, r, L, r_range);
func_n_r = @(r) helper_n(pp_n_r, r, L, r_range);

nr_func.n_r = func_n_r;
nr_func.dn_dr = func_dn_dr;
end

function n = helper_n(pp, r, wl_len, r_range)
n = nan(wl_len, length(r));
valid_idx = r >= r_range(1) & r <= r_range(2);
pos_idx = r > r_range(2);
neg_idx = r < r_range(1);
n(:, valid_idx) = ppval(pp, r(valid_idx));
n(:, pos_idx) = 1;
n(:, neg_idx) = nan;
end

function n = helper_dn(pp, r, wl_len, r_range)
der_pp = fnder(pp);
n = nan(wl_len, length(r));
valid_idx = r >= r_range(1) & r <= r_range(2);
pos_idx = r > r_range(2);
neg_idx = r < r_range(1);
n(:, valid_idx) = fnval(der_pp, r(valid_idx));
n(:, pos_idx) = 0;
n(:, neg_idx) = nan;
end