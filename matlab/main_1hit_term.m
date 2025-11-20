clear; close all; clc;

b_km = Constants.EARTH_R_km + linspace(0, 100, 101)';   % [B x 1]
t_deg = linspace(-12, 12, 41)';                         % [T x 1]
spec_wl_nm = (400:5:780)';                              % [K x 1]
moon_y_arcmin = linspace(0, 80, 241)';

B = length(b_km);
T = length(t_deg);
L = length(spec_wl_nm);
M = length(moon_y_arcmin);

SUN_R_arcmin = atand(Constants.SUN_R_km/Constants.EARTH_SUN_DIST_km)*60;

% Generate sun disk sample points
[sun_r_theta_arcmin, w] = generate_sun_half_disk_sample(SUN_R_arcmin, 91);
N = length(w);

% Atmosphere model
atmo = generate_atmosphere_model();

% Curve arc-length parameter grid
ds = 10;
s_km = (0:ds:2800)';    % [S x 1]
S = length(s_km);
s_full_km = [-wrev(s_km(2:end)); s_km];

% Create radical refractive functions
nr_func = create_nr_func(mean(spec_wl_nm), atmo);

% Integrate and get the raypath
fprintf('Computing direct ray paths...');
raypath_half_store = compute_raypath(nr_func, b_km, s_km); % [B x S x 6]
fprintf('done!\n');

% Solve for sun and moon position
sun_sol = solve_raypath(raypath_half_store, Constants.EARTH_SUN_DIST_km);   % [B x 4]
sun_sol(:,1) = -sun_sol(:,1);
sun_sol(:,4) = mod(pi - sun_sol(:,4), 2*pi);
moon_sol = solve_raypath(raypath_half_store, Constants.EARTH_MOON_DIST_km); % [B x 4]

% Concatenate half paths to full paths
raypath_full_store = [fliplr(raypath_half_store(:,2:end,:)), raypath_half_store];
raypath_full_store(:,1:S,[1,6]) = -raypath_full_store(:,1:S,[1,6]);
raypath_full_store(:,1:S,4) = pi - raypath_full_store(:,1:S,4);

% Compute raypath for all scattering directions
fprintf('Computing scattering ray paths...');
raypath_sct_store = compute_raypath(nr_func, b_km, s_km, t_deg); % [B x S x T x 6]
moon_sct_sol = solve_raypath(raypath_sct_store, Constants.EARTH_MOON_DIST_km);  % [B x T x 4]
fprintf('done!\n');

%%
% Compute Jacobian
fprintf('Computing Jacobian...');
ray_tangent = squeeze(atan2(raypath_full_store(:,:,6), raypath_full_store(:,:,5)));   % [B x 2S]
ray_inc = pi - sun_sol(:,4) + ray_tangent;
sun_emit = pi - sun_sol(:,4) - ray_tangent(:,end);
geo_mid = interp1((sun_emit(2:end)+sun_emit(1:end-1))/2, ...
    (ray_inc(2:end,:)-ray_inc(1:end-1,:)) ./ (sun_emit(2:end)-sun_emit(1:end-1)), sun_emit, ...
    'linear', 'extrap');
geo_mid = min(1./abs(geo_mid), 3.0);   % [B x 2S]

ray_tangent = squeeze(atan2(raypath_sct_store(:,end,:,6), raypath_sct_store(:,end,:,5)));   % [B x T]
geo_sct = nan(size(ray_tangent));     % [B x T]
for i = 1:B
    tmp_idx = ~isnan(ray_tangent(i,:));
    tmp_x = deg2rad(t_deg(tmp_idx));
    tmp_y = ray_tangent(i,tmp_idx);
    tmp_val = ppval(fnder(pchip(tmp_x, tmp_y), 1), tmp_x);
    geo_sct(i,tmp_idx) = 1./abs(tmp_val);
end
clear tmp_* ray_tangent ray_inc sun_emit;
fprintf('done\n');

% Compute optical depth tau
fprintf('Computing optical depth...');
sigma_air = rayleigh_sigma_bodhaine(spec_wl_nm, atmo);      % [L x B]
sigma_o3 = read_chappuis_data(spec_wl_nm, atmo.T);          % [L x B]

kappa_air_full = nan(B, 2*S-1, L);   % [B x 2S x L]
kappa_o3_full = nan(B, 2*S-1, L);    % [B x 2S x L]
for i = 1:L
    kappa_air_full(:,:,i) = interp1(atmo.z, atmo.air_number_density .* sigma_air(i,:), ...
        raypath_full_store(:,:,3) - Constants.EARTH_R_km, 'linear', 0);
    kappa_o3_full(:,:,i) = interp1(atmo.z, atmo.o3_number_density .* sigma_o3(i,:), ...
        raypath_full_store(:,:,3) - Constants.EARTH_R_km, 'linear', 0);
end
tau_mid = squeeze(cumtrapz(s_full_km*1e3, kappa_air_full + kappa_o3_full, 2));    % [B x 2S x L]
transit_factor_mid = exp(-tau_mid);     % [B x 2S x L]

kappa_air_sct = nan(B, S, T, L);
kappa_o3_sct = nan(B, S, T, L);
for i = 1:L
    kappa_air_sct(:,:,:,i) = interp1(atmo.z, atmo.air_number_density .* sigma_air(i,:), ...
        raypath_sct_store(:,:,:,3) - Constants.EARTH_R_km, 'linear', 0);
    kappa_o3_sct(:,:,:,i) = interp1(atmo.z, atmo.o3_number_density .* sigma_o3(i,:), ...
        raypath_sct_store(:,:,:,3) - Constants.EARTH_R_km, 'linear', 0);
end
tau_sct = squeeze(trapz(s_km*1e3, kappa_air_sct + kappa_o3_sct, 2));    % [B x T x L];
transit_factor_sct = exp(-tau_sct);     % [B x T x L]
fprintf('done!\n');

%%
cq = cosd(sun_r_theta_arcmin(:,2)/60);
sq = sind(sun_r_theta_arcmin(:,2)/60);
cr = cosd(sun_r_theta_arcmin(:,1)/60);
sr = sind(sun_r_theta_arcmin(:,1)/60);

cy = cosd(moon_y_arcmin/60);
sy = sind(moon_y_arcmin/60);

ray_terminator_phi = (0:pi/45:(2*pi-1e-8))';
cphi = cos(ray_terminator_phi);
sphi = sin(ray_terminator_phi);

z1_axis = -[cq.*sr, sq.*sr, -cr];
x1_axis = [cr, zeros(size(cq)), cq.*sr];
x1_axis = x1_axis ./ sqrt(dot(x1_axis, x1_axis, 2));
y1_axis = cross(z1_axis, x1_axis, 2);

moon0_dir = [zeros(M,1), sy, cy];   % [M x 3], original moon direction vectors
moon1_dir = cat(3, sum(reshape(moon0_dir, [1,M,3]) .* reshape(x1_axis, [N,1,3]), 3), ...
    sum(reshape(moon0_dir, [1,M,3]) .* reshape(y1_axis, [N,1,3]), 3), ...
    sum(reshape(moon0_dir, [1,M,3]) .* reshape(z1_axis, [N,1,3]), 3));  % [N x M x 3]
moon1_dir = reshape(moon1_dir, [], 3);  % [N*M x 3], in sun-earth frame

spec = zeros(N*M, L);
for i = 1:length(ray_terminator_phi)
    % Rotate by terminator angle (along local z-axis)
    moon2_dir = moon1_dir * [cphi(i), sphi(i), 0; -sphi(i), cphi(i), 0; 0, 0, 1];   % [N*M x 3]

    for j = 1:B
    % for j = 15
        curr_raypath = squeeze(raypath_full_store(j, :, :));    % [2S x 6]
        curr_raypath_idx = curr_raypath(:, 3) < max(b_km);

        curr_moon_sol = squeeze(moon_sct_sol(j,:,:));           % [T x 4]
        curr_moon_sol(:,4) = mod(curr_moon_sol(:,4)+pi, 2*pi)-pi;
        curr_moon_sol_idx = all(~isnan(curr_moon_sol), 2);

        y3_axis = [zeros(2*S-1,1), curr_raypath(:,[2,1])];
        y3_axis = y3_axis ./ curr_raypath(:,3);
        x3_axis = [ones(2*S-1,1), zeros(2*S-1,2)];
        z3_axis = cross(x3_axis, y3_axis, 2);
        fprintf('Computing raypath %d/%d ...', j, B);
        tic;
        for k = 1:2*S-1
            if ~curr_raypath_idx(k)
                continue;
            end
            % Transform to earth-terminator frame
            moon3_dir = moon2_dir * [x3_axis(k,:); y3_axis(k,:); z3_axis(k,:)]';
            moon3_yaw = atan2d(moon3_dir(:,1), moon3_dir(:,3));   % [N*M x 1]
            moon3_beta = asind(moon3_dir(:,2));   % [N*M x 1]

            if all(moon3_beta < min(t_deg(curr_moon_sol_idx))) || ...
                    all(moon3_beta > max(t_deg(curr_moon_sol_idx)))
                continue;
            end

            tmp_t = interp1(curr_moon_sol(curr_moon_sol_idx,4), t_deg(curr_moon_sol_idx), ...
                moon3_beta, 'linear', nan);     % [N*M x 1]
            tmp_idx = ~isnan(tmp_t);
            tmp_cos_theta = cosd(tmp_t(tmp_idx)) .* cosd(moon3_yaw(tmp_idx));   % [valid-N*M x 1]

            tmp_transmit_mid = squeeze(transit_factor_mid(j, k, :));            % [L x 1]
            tmp_geo_mid = geo_mid(j,k);

            tmp_sct_prob = squeeze(kappa_air_full(j, k, :)) * ds*1e3;           % [L x 1]
            tmp_transit_sct = interp1(t_deg, squeeze(transit_factor_sct(j,:,:)), ...
                tmp_t(tmp_idx), 'linear', 0);   % [valid-N*M x L]
            tmp_geo_sct = interp1(t_deg, geo_sct(j,:)', ...
                tmp_t(tmp_idx), 'linear', 0);   % [valid-N*M x 1]

            tmp_spec = 3/16/pi*(tmp_cos_theta.^2 + 1) .* ...
                (tmp_transmit_mid .* tmp_geo_mid .* tmp_sct_prob)' .* ...
                (tmp_transit_sct .* tmp_geo_sct);   % [valid-N*M x L]

            spec(tmp_idx, :) = spec(tmp_idx, :) + tmp_spec;
        end
        curr_raypath_duration = toc;
        fprintf('done. Elapsed: %.2f sec\n', curr_raypath_duration);

        if mod(j, 2) == 1
            figure(1); clf;
            plot(moon_y_arcmin, squeeze(sum(reshape(spec(:,[1,L]), [N,M,2]), 1)), '.-');

            figure(2); clf;
            tmp_rgb=spec_to_rgb([spec_wl_nm, squeeze(sum(reshape(spec,[N,M,L]),1))'], 'Mixed', true, 'Y', 2^0);
            imagesc(moon_y_arcmin, 1, reshape(tmp_rgb, [1,M,3]));
            drawnow;
        end
    end
    break;          % tmp break
end
