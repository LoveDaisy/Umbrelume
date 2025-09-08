function n = ciddor_n_air(wavelength_nm, T, p)
% n(λ, T, p) for dry air using Ciddor (1996) standard-air dispersion with density scaling.
% 
% Matrix convention: Returns n(wavelength, height) as [m×n] where m=length(wavelength_nm), n=length(T)
% 
% - Step 1: standard dry air (15 °C, 101325 Pa, ~450 ppm CO2):
%     n_std - 1 = 0.05792105/(238.0185 - 1/λ^2) + 0.00167917/(57.362 - 1/λ^2)
%     (λ in µm; result is absolute refractivity, not ×1e-6)
% - Step 2: scale by number density ~ p/T (ideal-gas) to get n(λ,T,p):
%     n - 1 ≈ (n_std - 1) * (p/101325) * (288.15/T)
% 
% Inputs:
%   wavelength_nm : wavelength array [nm] (m elements)  
%   T             : temperature array [K] (n elements)
%   p             : pressure array [Pa] (n elements, same size as T)
% 
% Output:
%   n : refractive index array [m×n] where n(i,j) = n(wavelength_nm(i), T(j), p(j))
% 
% References: Ciddor 1996; NIST EMT toolbox notes; refractiveindex.info Ciddor page.

% Ensure proper array shapes
lam_nm = wavelength_nm(:);  % Column vector [m×1]
T = T(:)';  % Row vector [1×n] 
p = p(:)';  % Row vector [1×n]

lam_um = lam_nm * 1e-3;  % [m×1]
s2 = 1.0 ./ (lam_um .^ 2); % [m×1]

% Ciddor standard-air dispersion (dry air, 15 °C, 101325 Pa, ~450 ppm CO2)
n_minus_1_std = (0.05792105 ./ (238.0185 - s2)) + (0.00167917 ./ (57.362 - s2));  % [m×1]

% Density scaling to T,p using broadcasting: [m×1] .* [1×n] -> [m×n]
scale = (p ./ 101325.0) .* (288.15 ./ T);  % [1×n]
n = 1.0 + n_minus_1_std .* scale;  % [m×1] .* [1×n] -> [m×n]
end