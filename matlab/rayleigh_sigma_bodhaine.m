function sigma_rayleigh = rayleigh_sigma_bodhaine(wavelength_nm, atmo)
% Molecular Rayleigh scattering cross-section per molecule σ_R(λ) for *dry air*
% following Bodhaine et al. (1999) recommended approach:
% 
%   σ(λ) = (24 * π^3 / (N^2 * λ^4)) * [ (n(λ)^2 - 1)^2 / (n(λ)^2 + 2)^2 ] * F_air(λ)
% 
% where:
%   - N = molecular number density [m^-3] (height-dependent from atmo struct)
%   - n(λ): refractive index of air (computed using Ciddor 1996 with T,p from atmo)
%   - F_air(λ): King depolarization factor for air (mixture of N2, O2, Ar, CO2)
% 
% Inputs:
%   wavelength_nm : wavelength(s) in nanometers
%   atmo          : atmospheric structure from us_standard_atmosphere_1976
%                   with fields: z, T, p, air_number_density, o3_number_density
% 
% Return:
%   σ_R(λ,z) in m^2 per molecule
%   - Single height: column vector [m×1] matching wavelength_nm
%   - Multiple heights: matrix [m×n] where m=length(wavelength_nm), n=length(atmo.z)
%   - Matrix convention: σ_R(i,j) = σ_R(wavelength_nm(i), atmo.z(j))
% 
% References:
%   - Core σ formula / King factor discussion (Bodhaine 1999).
%   - Ciddor (1996) refractive index with T,p dependence.
%   - Mixture King factor using Bates (1984) components as compiled/used in later ATBDs.

lam_nm = wavelength_nm(:);  % Ensure column vector [mx1]
lam_um = lam_nm * 1e-3;
lam_m = lam_nm * 1e-9;

T = atmo.T(:)';
p = atmo.p(:)';
N_air = atmo.air_number_density(:)';      % [1xn]

% 1) Refractive index using Ciddor (1996) with T,p dependence
% This replaces the standard air approximation with height-dependent values
% Returns [n_lam × n_z] matrix using broadcasting
ns = ciddor_n_air(lam_nm, T, p);    % [mxn]

% 2) King factor for air mixture F_air(λ) (λ in µm)
%    Component King factors (Bates 1984; used e.g. in AerGOM / GOMOS ATBD):
%      F_N2 = 1.034 + 3.17e-4 / λ^2
%      F_O2 = 1.096 + 1.385e-3 / λ^2 + 1.448e-4 / λ^4
%      F_Ar = 1.0
%      F_CO2 = 1.15
ln2 = 1.0 ./ (lam_um .^ 2);     % [mx1]
ln4 = ln2 .^ 2;
F_N2 = 1.034 + 3.17e-4 * ln2;
F_O2 = 1.096 + 1.385e-3 * ln2 + 1.448e-4 * ln4;
F_Ar = 1.0;
F_CO2 = 1.15;

% Volume fractions in % by volume (CO2 at 300 ppm → 0.03%)
C_N2 = 78.084;
C_O2 = 20.946;
C_Ar = 0.934;
C_CO2 = 0.03;

F_air = (C_N2 .* F_N2 + C_O2 .* F_O2 + C_Ar .* F_Ar + C_CO2 .* F_CO2) ./ (C_N2 + C_O2 + C_Ar + C_CO2);  % [mx1]

% 3) Rayleigh σ(λ,z) [m^2] per molecule using height-dependent parameters
% Use broadcasting for efficient computation: [m×1] and [1×n] -> [m×n]

% Bodhaine formula with broadcasting
ns2 = ns .^ 2;                                          % [mxn]
refractive_term = ((ns2 - 1.0) ./ (ns2 + 2.0)) .^ 2;    % [mxn]

% Broadcasting: [1×1] ./ ([1×n].^2 .* [m×1].^4) -> [m×n]
geometric_term = 24.0 * pi^3 ./ (N_air.^2 .* (lam_m.^4));  % [mxn]

% Final result: [mxn] .* [mxn] .* [mx1] -> [mxn]
sigma_rayleigh = geometric_term .* refractive_term .* F_air;
end
