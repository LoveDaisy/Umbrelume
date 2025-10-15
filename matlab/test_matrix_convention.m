% Test script to verify the wavelength×height matrix convention
% for ciddor_n_air and rayleigh_sigma_bodhaine functions

clear; clc;

% Test parameters
wavelengths = [400, 500, 600]';  % 3 wavelengths [nm]
heights = [0, 10, 20];           % 3 heights [km]

% Create atmospheric profile
atmo = us_standard_atmosphere_1976(heights);

fprintf('=== Matrix Convention Test ===\n');
fprintf('Wavelengths: [%s] nm\n', sprintf('%.0f ', wavelengths));
fprintf('Heights: [%s] km\n', sprintf('%.0f ', heights));
fprintf('Expected matrix size: %d×%d (wavelength×height)\n\n', length(wavelengths), length(heights));

% Test ciddor_n_air function
fprintf('1. Testing ciddor_n_air function:\n');
n_air = ciddor_n_air(wavelengths, atmo.T, atmo.p);
fprintf('   Result size: %d×%d\n', size(n_air, 1), size(n_air, 2));
fprintf('   n_air(400nm, 0km) = %.8f\n', n_air(1,1));
fprintf('   n_air(600nm, 20km) = %.8f\n', n_air(3,3));

% Test rayleigh_sigma_bodhaine function  
fprintf('\n2. Testing rayleigh_sigma_bodhaine function:\n');
sigma_R = rayleigh_sigma_bodhaine(wavelengths, atmo.T, atmo.p, atmo.air_number_density);
fprintf('   Result size: %d×%d\n', size(sigma_R, 1), size(sigma_R, 2));
fprintf('   σ_R(400nm, 0km) = %.3e m²\n', sigma_R(1,1));
fprintf('   σ_R(600nm, 20km) = %.3e m²\n', sigma_R(3,3));

% Test single height case
fprintf('\n3. Testing single height case:\n');
atmo_single = us_standard_atmosphere_1976(10);  % Single height
sigma_R_single = rayleigh_sigma_bodhaine(wavelengths, atmo_single.T, atmo_single.p, atmo_single.air_number_density);
fprintf('   Result size: %d×%d\n', size(sigma_R_single, 1), size(sigma_R_single, 2));
fprintf('   Should be column vector matching wavelength input\n');

fprintf('\n=== Test completed ===\n');