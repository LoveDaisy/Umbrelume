clear; close all; clc;

wavelength_nm = 400:2:780;
z_grid_km = 0:.1:90;
r_grid_arcmin = 0:.2:70;

% Atmosphere model
atmo = us_standard_atmosphere_1976(z_grid_km);
% Simplified o3 profile
vmr_ppmv = ozone_vmr_profile_simple(z_grid_km);
atmo.o3_number_density = (vmr_ppmv * 1e-6) .* atmo.air_number_density;

E_sun = read_astm_e490_data(wavelength_nm);

% Direct term: transmitance
T_dir = direct_term(wavelength_nm, atmo, r_grid_arcmin);

%%
E_dir_rgb = spec_to_rgb([wavelength_nm', E_sun .* T_dir], 'mixed', true, 'Y', 0.85);

figure(1);
imagesc(r_grid_arcmin, 1, reshape(E_dir_rgb, 1,[],3));
