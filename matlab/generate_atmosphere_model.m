function atmo = generate_atmosphere_model()
z_km = (0:.1:110)';
atmo = us_standard_atmosphere_1976(z_km);
ozone_vmr = ozone_vmr_profile_simple(z_km)' * 1.1 * 1e-6;
atmo.o3_number_density = ozone_vmr .* atmo.air_number_density; % Simplified o3 profile
end