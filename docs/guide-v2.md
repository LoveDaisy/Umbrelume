# 月食颜色模拟（环积分 & 角半径 r 版本）设计文档 v2.0

> 本版按你提出的“同一像素的不同路径叠加”重构：把直达折射项与散射项明确分离；前者在球对称+点太阳极限下退化为单一掩星平面的一条（或太阳盘平均的一簇）射线；后者对整圈方位角 phi 与高度 z 做体积积分。

---

## 1. 目标与范围

* 目标：在“月球距离处”的屏幕上，以角半径 r（弧度/弧分）为自变量，计算从本影中心向外的光谱辐亮度与颜色，并再现“中心暗红—外缘蓝绿—半影灰白”的过渡。
* 范围：

  * 直达折射：球对称大气 + 可选太阳有限角与肢变暗；
  * 散射：单次散射（瑞利必选；气溶胶可选），整圈 phi 与高度 z 体积分；
  * 月面反照率先取中性常数（或后续接入 ROLO / 光度模型）。

---

## 2. 物理分解（像素处的所有路径叠加）

令像素位于角半径 r；其光谱辐亮度分解为

L\_lambda(r) = L\_dir\_lambda(r) + L\_sca\_lambda(r).

### 2.1 直达折射项（unscattered, refracted）

* 在球对称介质中，总偏折角 alpha 仅由冲量参数/切向高度 h 决定；给定 r（小角近似 r≈alpha），存在唯一 h\*(r)。若考虑有限太阳盘，则对太阳盘面元 u 做平均，求解 h\*(r,u) 使 alpha(h) = r - delta\_alpha(u)。
* 光谱透过：T\_dir(lambda; h) = exp( - (tau\_R + tau\_O3 + tau\_aer) ).
* 直达项：
  L\_dir\_lambda(r) = R\_moon(lambda) \* avg\_u\[ E\_sun(lambda,u) \* T\_dir(lambda; h\*(r,u)) ].
  注：在此极限下不需要对 phi 做环积分；问题严格在单一掩星平面内求解。

### 2.2 散射项（single-scattered）

* 对地球周缘整圈 phi 与高度 z 的体元，累积“入射 → 散射 → 出射”的贡献：

  L\_sca\_lambda(r) = R\_moon(lambda) \* ∫*{phi=0..2pi} ∫*{zmin..zmax} beta\_sca(lambda,z,phi) \* J\_lambda(z,phi) \* P(lambda, theta(z,phi→r)) \* T\_out(lambda; z→r,phi) dz dphi.

  其中 J\_lambda(z,phi) = avg\_u\[ E\_sun(lambda,u) \* T\_in(lambda; u→z,phi) ].
* beta\_sca：体散射系数（瑞利：由气体数密度给出；气溶胶：beta = SSA \* beta\_ext）。
* P(lambda, theta)：相函数（瑞利：3/4 \* (1 + cos^2 theta)；气溶胶可用 Henyey–Greenstein，参数 g≈0.6–0.8）。
* T\_in / T\_out：入射/出射路径消光（瑞利 + 臭氧 + 气溶胶）。

---

## 3. 数据与参数

必需：

1. 太阳 AM0 光谱 E\_sun(lambda)（ASTM E490）。
2. 标准大气 p(z), T(z)（US Standard Atmosphere 1976）。
3. 空气折射率 n(lambda, T, p, RH, x\_CO2)（Ciddor / Edlén）。
4. 瑞利散射：sigma\_R(lambda) 与气体数密度 N(z)；推荐 Bodhaine 参数化。
5. 臭氧吸收（Chappuis 450–750 nm 为关键）：sigma\_O3(lambda, T)（Serdyuchenko / Gorshelev 数据集）。
6. 颜色变换：CIE 1931 2° CMFs；sRGB（D65）。

可选/扩展：

* 气溶胶消光与相函数（Ångström 模型；HG/g、SSA）；
* 太阳盘与肢变暗；
* 月面光度/光谱模型（ROLO）；
* Chapman 函数近似作斜程柱含量的快速估计（亦可全数值积分）。

网格建议：

* lambda：380–780 nm，步长 5 nm（或在 520–650 nm 加密到 2 nm）。
* z：10–40 km 主体（可扩 5–80 km），步长 0.5–1 km。
* phi：N\_phi = 360–1440（1°–0.25°）。
* r：步长 0.5′–1′（覆盖 0 → \~80′）。
* 太阳盘采样 K=1（点太阳）或 10–50（含肢变暗）。

---

## 4. 几何要点与数值射线

* 小角近似：月球距离处角半径 r ≈ alpha（弧度）。
* 采用 Auer–Standish 型射线追迹在分层 n(z) 中计算 alpha(h)，对每个 r（二分/样条）反解唯一 h\*(r)。
* 入/出路径的光学厚度：沿几何路径分片积分 tau = ∫ k(lambda,s) ds。

---

## 5. 伪代码（直达 + 散射，全流程）

```text
INPUTS:
  lambda_grid, E_sun[lambda]
  Atmo: p(z), T(z), n_O3(z), (optional) aerosol params (beta_ext, SSA, g)
  Refractive index model n(lambda,z)
  Grids: r_bins, z_grid, phi_grid, solar_disk_samples u (delta_alpha(u), w_u)
  Options: use_solar_disk, use_aerosol

# 预计算：折射表 alpha(h)（参考 lambda0=550 nm；若要色散可多波长）
for h in h_grid: alpha_tab[h] = ray_bending(h; n(lambda0,z))
assert alpha_tab is monotonic

# 分子常数与截面：sigma_R(lambda), sigma_O3(lambda,T)

# —— 直达折射项 ——
for r_bin in r_bins:
  acc_dir[lambda] = 0; W = 0
  for u in solar_disk_samples (K=1 for point Sun):
    r_eff = r_center(r_bin) - delta_alpha(u)
    h_star = invert_alpha(r_eff, alpha_tab)
    if exists(h_star):
      tau_dir[lambda] = optical_depth_along_refracted_path(h_star, lambda)
      T_dir[lambda] = exp(-tau_dir[lambda])
      acc_dir[lambda] += w_u * E_sun[lambda] * T_dir[lambda];  W += w_u
  Ein_dir[lambda, r_bin] = (W>0) ? acc_dir[lambda]/W : 0

# —— 单次散射项 ——
for r_bin in r_bins:
  acc_sca[lambda] = 0
  for phi in phi_grid:
    for z in z_grid:
      # 入射：太阳盘平均的体元辐照 J_lambda
      J[lambda] = 0; W=0
      for u in solar_disk_samples:
        tau_in[lambda] = optical_depth_in(u -> (z,phi), lambda)
        J[lambda] += w_u * E_sun[lambda] * exp(-tau_in[lambda]);  W += w_u
      J[lambda] = (W>0) ? J[lambda]/W : 0

      # 体散射系数与相函数
      beta = beta_sca(lambda, z)  # 瑞利：由 N(z) 与 sigma_R(lambda)；气溶胶：SSA*beta_ext
      theta = scattering_angle((z,phi) -> r_bin)
      P = phase_function(lambda, theta)  # 瑞利或 HG

      # 出射消光
      tau_out[lambda] = optical_depth_out((z,phi) -> r_bin, lambda)

      acc_sca[lambda] += beta * J[lambda] * P * exp(-tau_out[lambda]) * Delta_z * Delta_phi
  Ein_sca[lambda, r_bin] = acc_sca[lambda]

# —— 合成与色彩 ——
for r_bin in r_bins:
  L_lambda[lambda] = (Ein_dir[lambda, r_bin] + Ein_sca[lambda, r_bin]) * R_moon(lambda)
  (X,Y,Z) = sum_over_lambda( L_lambda[lambda] * CMF1931(lambda) * Delta_lambda )
  sRGB[r_bin] = XYZ_to_sRGB_D65(X,Y,Z)

OUTPUT:
  r -> sRGB 颜色带；可选：r -> L_lambda 光谱、亮度 Y、(x,y) 色度剖面
```

---

## 6. 实现细节与数值建议

* alpha(h) 采样：在 h 或 alpha 空间自适应采样，保证边缘（小 alpha）与中心（大 alpha）都分辨充分。
* 数值稳定：当低层 z 导致 tau 很大时可直接忽略该路径；入/出路径剖分时用弧长均匀步长或按密度加权步长。
* 灵敏度：蓝绿缘强度对 O3 柱含量、20–30 km 的温度/密度与气溶胶 SSA/g 敏感。
* 能量与度量：直达项按太阳盘与掩星平面平均；散射项是体积分，单位保持为辐亮度（不要把整圈 phi 当作重复计数的“增亮因子”）。

---

## 7. 验证用例

1. 无臭氧：关闭 O3 → 外缘蓝绿消失，仅暖红—灰过渡。
2. 提高 O3 柱：蓝绿带更显著、位置略外移（更高 z / 更小 alpha）。
3. 平流层气溶胶增加（g≈0.7, SSA≈0.95）：整体变暗、色相偏红，散射项占比上升。
4. 有限太阳盘：边缘更平滑；半影区亮度/色相随可见太阳盘几何分数单调过渡。

---

## 8. 快速起步参数（建议）

* 波段：380–780 nm / 步长 5 nm；
* z：10–40 km / 步长 1 km；phi：720 段；r：步长 1′；
* O3：总柱 320 DU；剖面峰值 25 km，高斯宽度 8 km；
* 气溶胶：能见度 50 km（AOD≈0.05\@550 nm），SSA=0.95，HG g=0.7；
* 太阳盘：K=1（先关肢变暗）；若 tau > 20 则忽略该路径。

---

## 9. 已知局限与后续

* 本版散射只到一次散射；多次散射与偏振效应留作扩展。
* 大气视作球对称（直达项）+“phi 只用于散射体积分”的近似；若要非均匀臭氧/气溶胶/云，可把 p,T,n\_O3,beta 做成 (z,phi) 场并相应修改入/出路径。
* 未考虑地形遮挡与地球扁率；高精度需 DEM 与再分析场。
* 月面 BRDF 未显式建模，渲染级应用可接入 ROLO/光度模型。

---

## 10. 交付结构

* atmo.py：标准大气、n(z) 计算；sigma\_R / sigma\_O3 数据读入。
* raytrace.py：Auer–Standish 射线偏折 alpha(h) 与反解 h\*(r\[,u])。
* transmit.py：入/出路径光学厚度；Beer–Lambert。
* scatter.py：单次散射体积分（beta, P, J, T\_out）。
* color.py：光谱→XYZ→sRGB（D65）。
* demo.ipynb：r-颜色带、亮度/色度剖面、参数扫描。



