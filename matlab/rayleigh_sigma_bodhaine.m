function sigma_rayleigh = rayleigh_sigma_bodhaine(wavelength_nm)
% Molecular Rayleigh scattering cross-section per molecule σ_R(λ) for *dry air*
% following Bodhaine et al. (1999) recommended approach:
% 
%   σ(λ) = (24 * π^3 / (N_s^2 * λ^4)) * [ (n_s(λ)^2 - 1)^2 / (n_s(λ)^2 + 2)^2 ] * F_air(λ)
% 
% where:
%   - N_s = molecular number density at STP (dry air, 1013.25 hPa, 288.15 K)
%           ≈ 2.5469×10^25 m^-3 (2.5469×10^19 cm^-3).  [Colour/ CODATA]
%   - n_s(λ): refractive index of standard dry air incl. CO2 correction (Bodhaine-style).
%   - F_air(λ): King depolarization factor for air (mixture of N2, O2, Ar, CO2).
% 
% Inputs:
%   lam_nm : wavelength(s) in nanometers
% 
% Return:
%   σ_R(λ) in m^2 per molecule
% 
% References:
%   - Core σ formula / King factor discussion (Bodhaine 1999).
%   - n_s(λ) parameterization with CO2 correction (λ in µm).
%   - Mixture King factor using Bates (1984) components as compiled/used in later ATBDs.
%   - Standard molecular density N_s at STP (288.15 K).

lam_nm = wavelength_nm;
lam_um = lam_nm * 1e-3;
lam_m = lam_nm * 1e-9;

% 1) Refractive index of standard dry air (Bodhaine-style, λ in µm)
%    ns = 1 + 1e-8 * [8060.51 + 2480990/(132.274 - 1/λ^2) + 17455.7/(39.32957 - 1/λ^2)]
%    × [1 + 0.54(CCO2 - 0.0003)] ; here take CCO2=300 ppm → correction factor = 1
inv_l2 = 1.0 ./ (lam_um .^ 2);
ns = 1.0 + 1e-8 * (8060.51 + 2480990.0 ./ (132.274 - inv_l2) + 17455.7 ./ (39.32957 - inv_l2));
% (optional CO2 correction if you later expose it)
% CCO2 = 0.0003  % 300 ppm
% ns *= (1.0 + 0.54 * (CCO2 - 0.0003))

% 2) King factor for air mixture F_air(λ) (λ in µm)
%    Component King factors (Bates 1984; used e.g. in AerGOM / GOMOS ATBD):
%      F_N2 = 1.034 + 3.17e-4 / λ^2
%      F_O2 = 1.096 + 1.385e-3 / λ^2 + 1.448e-4 / λ^4
%      F_Ar = 1.0
%      F_CO2 = 1.15
ln2 = 1.0 ./ (lam_um .^ 2);
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

F_air = (C_N2 .* F_N2 + C_O2 .* F_O2 + C_Ar .* F_Ar + C_CO2 .* F_CO2) ./ (C_N2 + C_O2 + C_Ar + C_CO2);

% 3) Loschmidt number (STP, 288.15 K) in m^-3
N_s = 2.546902e25;  % consistent with 2.5469e19 cm^-3  (Colour docs / CODATA)

% 4) Rayleigh σ(λ) [m^2] per molecule
ns2 = ns .^ 2;
term = ((ns2 - 1.0) ./ (ns2 + 2.0)) .^ 2;
sigma_rayleigh = (24.0 * pi .^ 3 ./ (N_s .^ 2 * lam_m .^ 4)) .* term .* F_air;
end
