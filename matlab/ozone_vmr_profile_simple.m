function vmr = ozone_vmr_profile_simple(z_grid_km)
base = 0.04;

peak0 = 8.0;
z0 = 25.0;
s0 = 7.0;

peak1 = 0.6;
z1 = 35.0;
s1 = 12.0;

z2 = 8.0;

gauss = peak0 * exp(-0.5 * ((z_grid_km - z0) ./ s0) .^ 2);
upper = peak1 * min(exp(-(z_grid_km - z1) ./ s1), (z_grid_km - z1) ./ s0 + 1);
upper = max(upper, 0);
vmr = base + gauss + upper;

scale = max(min(0.5 / (z1 - z2) .* (z_grid_km - z2), 1.0), 0.5);
vmr = vmr .* scale;
end