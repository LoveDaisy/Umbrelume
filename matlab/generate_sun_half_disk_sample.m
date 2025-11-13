function [r_theta_arcmin, w] = generate_sun_half_disk_sample(sun_r_arcmin, n_bands)
total_num = n_bands^2;
r_theta_arcmin = nan(total_num, 2);  % [r_arcmin, theta_arcmin]
w = nan(total_num, 1);          % weights

r_grid_arcmin = linspace(0, sun_r_arcmin, n_bands);
count = 0;
for j = 1:n_bands
    band_sample_num = 2*j-1;
    theta_grid_rad = linspace(-pi/2, pi/2, band_sample_num)';   % From -90 to +90 deg

    curr_w = ones(size(theta_grid_rad));
    curr_w(1) = 0.5; curr_w(end) = 0.5;

    r_theta_arcmin(count + (1:band_sample_num), :) = ...
        [r_grid_arcmin(j) * ones(band_sample_num,1), rad2deg(theta_grid_rad)*60];
    w(count + (1:band_sample_num)) = curr_w;
    count = count + band_sample_num;
end
end