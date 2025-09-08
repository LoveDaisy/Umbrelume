function n = ciddor_n_air(wavelength_nm, T, p)
% n(λ, T, p) for dry air using Ciddor (1996) standard-air dispersion with density scaling.
% - Step 1: standard dry air (15 °C, 101325 Pa, ~450 ppm CO2):
%     n_std - 1 = 0.05792105/(238.0185 - 1/λ^2) + 0.00167917/(57.362 - 1/λ^2)
%     (λ in µm; result is absolute refractivity, not ×1e-6)
% - Step 2: scale by number density ~ p/T (ideal-gas) to get n(λ,T,p):
%     n - 1 ≈ (n_std - 1) * (p/101325) * (288.15/T)
%   NOTE: This omits humidity/compressibility corrections; for our月食折射用途误差在 ppm 量级。
% References: Ciddor 1996; NIST EMT toolbox notes; refractiveindex.info Ciddor page. 

% --- Ciddor 1996 折射率（干空气；按密度线性缩放），λ 输入 nm ---
lam_nm = wavelength_nm;
lam_um = lam_nm * 1e-3;
s2 = 1.0 ./ (lam_um^2); % [µm^-2]

% Ciddor standard-air dispersion (dry air, 15 °C, 101325 Pa, ~450 ppm CO2)
n_minus_1_std = (0.05792105 ./ (238.0185 - s2)) + (0.00167917 ./ (57.362 - s2));

% density scaling to T,p
scale = (p ./ 101325.0) .* (288.15 ./ T);
n = 1.0 + n_minus_1_std .* scale;
end