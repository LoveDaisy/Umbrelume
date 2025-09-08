function atom = us_standard_atmosphere_1976(z_grid_km)
% Build COESA 'U.S. Standard Atmosphere 1976' (0–86 km) on a geometric altitude grid.
% 
% Implements the layer-wise geopotential formulation:
%   - Convert geometric altitude z -> geopotential h = Re * z / (Re + z)
%   - For each layer b with base (h_b, T_b, p_b, lapse L_b), compute:
%       If L_b != 0:
%         T = T_b + L_b * (h - h_b)
%         p = p_b * (T / T_b) ** (-g0 * M / (R* * L_b))
%       If L_b == 0:
%         T = T_b
%         p = p_b * exp(-g0 * M * (h - h_b) / (R* * T_b))
%   - Density from ideal gas: ρ = p / (R_specific * T);
%     Number density N = p / (k_B * T)
% 
% Units:
%   z, h in km (internally converted to m where needed); L_b in K/m.
% 
% Constants (USSA76 / SI):
%   Re = 6356.766 km (Earth radius for geopotential)  [COESA 1976]
%   g0 = 9.80665 m/s^2
%   R* = 8.31432 J/(mol·K)
%   M  = 28.9644e-3 kg/mol (dry air mean molecular weight)
%   kB = 1.380649e-23 J/K
% 
% Layer bases (geopotential, km), lapse rates (K/km), and base temperatures (K) per USSA76:
%   h_b: [0, 11, 20, 32, 47, 51, 71, 84.852]
%   L_b: [-6.5, 0, +1.0, +2.8, 0, -2.8, -2.0]  (K/km)
%   T_b: [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946]
%   p_b[0] = 101325 Pa; other p_b computed recursively from layer formulas.
% 
% References:
%   - COESA / NOAA “U.S. Standard Atmosphere, 1976” formulas & layer limits.
%   - Quick lookup of layer table / lapse rates (cross-check).

z = z_grid_km;

% Constants
Re_km = 6356.766;
g0 = 9.80665;
Rstar = 8.31432;
M = 28.9644e-3;
kB = 1.380649e-23;

% geometric -> geopotential (km)
h = (Re_km * z) ./ (Re_km + z);
h_m = h * 1000.0;

% Layer bases (geopotential, km)
hb = [0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852];
% Lapse rates [K/km] -> [K/m]
Lb_km = [-6.5, 0.0, +1.0, +2.8, 0.0, -2.8, -2.0];
Lb = Lb_km / 1000.0;
% Base temperatures [K]
Tb = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946];
% Base pressures [Pa] (p_b[0] known; others by recurrence)
pb = zeros(size(Tb));
pb(1) = 101325.0;
for i = 1:length(hb)-1
    if Lb(i) == 0.0
        dh = (hb(i + 1) - hb(i)) * 1000.0;
        pb(i + 1) = pb(i) * exp(-g0 * M * dh / (Rstar * Tb(i)));
    else
        Tb_next = Tb(i + 1);
        pb(i + 1) = pb(i) * (Tb_next / Tb(i)) .^ (-g0 * M / (Rstar * Lb(i)));
    end
end

% For each altitude, pick its layer and compute T, p
T = nan(size(h_m));
p = nan(size(h_m));
for idx = 1:length(h_m)
    hm = h_m(idx);
    % clamp to model top
    if hm/1000.0 >= hb(end)
        ib = length(hb) - 1;  % use top layer
        hm = min(hm, hb(end) * 1000.0 - 1e-6);
    else
        ib = find(hm < hb * 1000.0, 1, 'first') - 1;
        ib = min(max(ib, 1), length(hb) - 1);
    end
    h0 = hb(ib) * 1000.0;
    L = Lb(ib);
    T0 = Tb(ib);
    p0 = pb(ib);
    if L == 0.0
        T(idx) = T0;
        p(idx) = p0 * exp(-g0 * M * (hm - h0) / (Rstar * T0));
    else
        T(idx) = T0 + L * (hm - h0);
        p(idx) = p0 * (T(idx) / T0) .^ (-g0 * M / (Rstar * L));
    end
end

% Number density (dry air)
air_number_density = p ./ (kB * T);

% Placeholder ozone number density (all zeros; plug real profile later)
o3_number_density = nan(size(air_number_density));

atom.z = z;
atom.T = T;
atom.p = p;
atom.air_air_number_density = air_number_density;
atom.o3_number_density = o3_number_density;
end