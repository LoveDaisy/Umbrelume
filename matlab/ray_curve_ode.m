function df = ray_curve_ode(s, f, func_n_r, func_dn_dr)
% RAY_CURVE_ODE ODE function for ray curve computation in 2D.
%   Inputs:
%       s:              Arc-length parameter (not used, but required by ODE solver)
%       f:              [4 x 1] state vector [x; y; tx; ty]
%       func_n_r:      Function handle of refractive index n(r)
%       func_dn_dr:    Function handle of dn/dr(r)
%   Outputs:
%       df:             [4 x 1] derivative of state vector

x = f(1); y = f(2); tx = f(3); ty = f(4);

r  = sqrt(x^2 + y^2);
t2 = tx^2 + ty^2;

% Check for |t| = 1
if abs(t2-1) > 1e-12
    sc = 1 / sqrt(t2);
    tx = tx*sc;
    ty = ty*sc;
end

% grad ln n
nl  = func_n_r(r); 
dln = func_dn_dr(r) ./ nl;
gx = dln .* (x./r);
gy = dln .* (y./r);

% proj_perp = (I - t t^T) g
dot = gx.*tx + gy.*ty;
kx  = gx - dot.*tx;
ky = gy - dot.*ty;
df = [ tx; ty;  kx; ky ];
end