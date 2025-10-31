function [func_n_r, func_dn_dr] = create_func_n(wavelength_nm, z_grid_km)
wavelength_nm = wavelength_nm(:);   % Column [mx1]
z_grid_km = z_grid_km(:)';          % Row    [nx1]

EARTH_R_km = 6371;

atmo = us_standard_atmosphere_1976(z_grid_km);
n = ciddor_n_air(wavelength_nm, atmo.T, atmo.p);
pp_n_r = pchip(z_grid_km + EARTH_R_km, n);
func_dn_dr = @(r) helper_dn(pp_n_r, r, length(wavelength_nm), max(z_grid_km)+EARTH_R_km);
func_n_r = @(r) helper_n(pp_n_r, r, length(wavelength_nm), max(z_grid_km)+EARTH_R_km);
end

function n = helper_n(pp, r, wl_len, r_max)
n = nan(wl_len, length(r));
n(:, r < r_max) = ppval(pp, r(r < r_max));
n(:, r >= r_max) = 1;
end

function n = helper_dn(pp, r, wl_len, r_max)
der_pp = fnder(pp);
n = nan(wl_len, length(r));
n(:, r < r_max) = fnval(der_pp, r(r < r_max));
n(:, r >= r_max) = 0;
end