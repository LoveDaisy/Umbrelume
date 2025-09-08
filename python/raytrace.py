# ──────────────────────────────────────────────────────────────────────────────
# raytrace.py —— α(h) 表 + 反解 + 太阳盘采样
# ──────────────────────────────────────────────────────────────────────────────
from __future__ import annotations
from dataclasses import dataclass
from typing import List
import numpy as np

from atmo import Atmosphere1D, ciddor_n_air

# 几何常量 (用于 r = R_E + z)
RE_KM = 6371.0
# 折射几何参考波长 (与 transmit.py 的 LAMBDA0_NM 保持一致)
LAMBDA0_NM = 550.0


@dataclass
class AlphaTable:
    h_km: np.ndarray        # 切向几何高度网格 (km) , 升序
    alpha_rad: np.ndarray   # 对应的总弯折角 α(h) (弧度) , 严格单调下降 (h↑ → α↓)


def _prepare_nr_profiles(atmo: Atmosphere1D, lam0_nm: float):
    """
    预计算: r(z)、n(z;λ0)、nr(z) 与 d ln n / dr (在 r 上的导数, 单位 1/m) . 
    参考弯折角积分式: α(a) = -2 a ∫ (d ln n/dr) / sqrt((n r)^2 - a^2) dr
    文献: Fjeldbo-Abel框架 (见 Syndergaard 2005 教程) 与 Auer–Standish 非奇异积分方法.  
    """
    z = np.asarray(atmo.z, dtype=float)
    r = (RE_KM + z) * 1000.0                         # [m]
    # 标准干空气折射率 (密度缩放) , 同 λ0:
    n = ciddor_n_air(lam0_nm, T=atmo.T, p=atmo.p)    # 与 atmo 网格对齐
    # 数值导数: d ln n / dr = (1/n) * dn/dr
    dn_dr = np.gradient(n, r)                        # [1/m]
    dlnn_dr = dn_dr / n
    nr = n * r
    return z, r, n, nr, dlnn_dr


def _alpha_from_h(atmo: Atmosphere1D, lam0_nm: float,
                  r: np.ndarray, n: np.ndarray, nr: np.ndarray, dlnn_dr: np.ndarray,
                  h_star_km: float, z_top_km: float | None = None) -> float:
    """
    对给定切向高度 h* 计算弯折角 α(h*). 
    实现: 按 Fjeldbo/Abel 的一维积分, 
      α(a) = -2 a ∫_{r_t}^{r_top} (d ln n / dr) / sqrt((n r)^2 - a^2) dr  , a = n(r_t) r_t
    采用分段中点积分避免 r_t 端点奇异; r_top 取大气顶 (atmo.z.max) . 
    """
    r_t = (RE_KM + h_star_km) * 1000.0
    if z_top_km is None:
        z_top_km = float(atmo.z.max())
    r_top = (RE_KM + z_top_km) * 1000.0

    # 需要在 r 上找到切向位置索引 (r 单调增)
    if r_t < r[0] or r_t > r[-1]:
        # h* 超出大气网格 (通常 < z_min 或 > z_max) , 返回 0 弯折
        return 0.0

    # 在 r 网格中定位 r_t 所在的区间 [i, i+1], 并从 r_t 积到 r_top
    i0 = int(np.searchsorted(r, r_t, side="right") - 1)
    # 构造分段边界: 首段以 r_t 为下界, 然后跟随 atmo.r 的层边界到顶
    r_edges = np.concatenate(([r_t], r[i0+1:][r[i0+1:] <= r_top], [r_top]))
    if len(r_edges) < 2:
        return 0.0

    # 在每个薄层 [r_j, r_{j+1}] 的中点上评估 integrand
    r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])

    # 插值所需量到 r_mid
    nr_mid = np.interp(r_mid, r, nr)
    dlnn_dr_mid = np.interp(r_mid, r, dlnn_dr)

    # 计算 a = n(r_t)*r_t (用线性插值得到切向处 n_t)
    n_t = float(np.interp(r_t, r, n))
    a = n_t * r_t

    # 被开方量: ((n r)^2 - a^2); 近切向可能很小, 做个数值地板
    denom = np.sqrt(np.maximum(nr_mid**2 - a**2, 1e-32))

    integrand = dlnn_dr_mid / denom
    dr = np.diff(r_edges)
    alpha = -2.0 * a * np.sum(integrand * dr)   # [rad], 应为正值 (dlnn_dr<0)

    # 数值健壮性: 极小的负零清为 0
    return float(max(alpha, 0.0))


def build_alpha_table(atmo: Atmosphere1D, lam0_nm: float = LAMBDA0_NM,
                      h_min_km: float | None = None, h_max_km: float | None = None,
                      dh_km: float = 0.5) -> AlphaTable:
    """
    生成 α(h) 表: 
      - h 网格默认为 atmo.z 覆盖 (或指定 [h_min, h_max]) , 步长 dh_km; 
      - 对每个 h* 用上面的数值积分计算 α(h*); 
      - 返回单调的 AlphaTable (如遇非严格单调, 做极小平滑确保可反解) . 
    参考: Auer & Standish (2000) “适用于任意天顶角的非奇异数值法”; RO/Abel 经典表达.  
    """
    # 预先准备剖面
    z, r, n, nr, dlnn_dr = _prepare_nr_profiles(atmo, lam0_nm)

    if h_min_km is None:
        h_min_km = float(max(z[0], 0.0))
    if h_max_km is None:
        h_max_km = float(z[-1])
    h_grid = np.arange(h_min_km, h_max_km + 1e-9, dh_km)

    alpha_list = []
    for h in h_grid:
        alpha = _alpha_from_h(atmo, lam0_nm, r, n, nr, dlnn_dr, h_star_km=h, z_top_km=h_max_km)
        alpha_list.append(alpha)
    alpha_arr = np.array(alpha_list, dtype=float)

    # 单调性修正: α(h) 随 h 增加应当单调下降 (数值噪声平滑)
    # 若出现微小上翘, 做一次“退火式”单调回归 (简单后向最小化)
    for i in range(1, len(alpha_arr)):
        if alpha_arr[i] > alpha_arr[i-1]:
            alpha_arr[i] = alpha_arr[i-1]

    return AlphaTable(h_km=h_grid, alpha_rad=alpha_arr)


def invert_h_from_alpha(alpha_table: AlphaTable, r_rad: float) -> float | None:
    """
    反解 h*(r): 在小角近似 r≈α 下, 以 α(h) 单调表反查 h (线性插值) . 
    若 r 超出表范围, 返回 None. 
    """
    a = alpha_table
    # α 应随 h 单调下降; 检查范围
    amin, amax = float(a.alpha_rad[-1]), float(a.alpha_rad[0])
    if not (amin <= r_rad <= amax):
        return None
    # 为保证单调, 反向插值: h(α)
    # 注意: alpha_rad 随 h 增加而下降 → 用 np.interp 需要升序自变量
    alpha_sorted = a.alpha_rad[::-1]   # 升序 (小→大)
    h_sorted = a.h_km[::-1]
    h_star = np.interp(r_rad, alpha_sorted, h_sorted)
    return float(h_star)


@dataclass
class SolarDiskSample:
    delta_alpha_rad: float  # 太阳盘样本在“径向方向上的角偏置” (投影到 r 轴)
    weight: float           # 该样本的面积权重 (含肢变暗可在外层改权)


def default_solar_disk_samples(K: int = 1,
                               solar_radius_rad: float = np.deg2rad(0.5/2.0),
                               limb_darkening: bool = False) -> List[SolarDiskSample]:
    """
    K=1 → 点太阳: δ=0, w=1
    K>1 → 近似“太阳盘在 r 方向的 1D 卷积”采样: 在圆盘上做环×角均匀取样, 
           将每个样本的位移投影到 r 轴: δ = ρ cosθ, 权重按环面积 / 样本数. 
    - solar_radius_rad: 太阳角半径 (默认 ~0.25° ≈ 15′ ≈ 0.00436 rad) 
    - limb_darkening: 若 True, 可在外部用 I(μ) 去修改 weight (这里保持简单) 
    说明: 这与严格 2D 卷积的等效一维采样是近似, 但足以平滑直达项边缘. 
    """
    if K <= 1:
        return [SolarDiskSample(0.0, 1.0)]

    # 选环数与每环样本数 (≈K) , 简单规则: N_ring^2 ≈ K
    n_ring = max(2, int(np.sqrt(K)))
    d_rho = solar_radius_rad / n_ring

    samples: List[SolarDiskSample] = []
    total_w = 0.0
    for j in range(1, n_ring + 1):
        rho = (j - 0.5) * d_rho
        # 每环角向样本数与半径成比例
        n_th = max(6, int(np.round(12 * j)))
        dtheta = 2.0 * np.pi / n_th
        # 该环的面积 ~ 2πρ dρ; 均分到 n_th 个样本
        w_per = (2.0 * np.pi * rho * d_rho) / n_th
        for k in range(n_th):
            theta = (k + 0.5) * dtheta
            delta = float(rho * np.cos(theta))  # 投影到“径向”方向
            w = w_per
            samples.append(SolarDiskSample(delta_alpha_rad=delta, weight=w))
            total_w += w

    # 归一化权重
    for s in samples:
        s.weight /= total_w
    return samples
