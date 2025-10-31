function df = ray_curve_ode(s, f, func_n_r, func_dn_dr)
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