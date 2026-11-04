# Lunar Eclipse Color

从太阳光谱、大气折射与沿程消光，计算月球距离处的光谱照度，再显示月食的颜色。
这是从空白重建的第一版；旧实现保留在 `legacy/*` 分支。本版未参考旧实现。

当前模型已跑通：红褐色本影、臭氧造成的浅蓝灰色外缘，以及向半影和正常阳光的过渡。
**目前的冷色外缘尚不能称为已经复现了明显的绿松石带。** 模型采用理想晴空大气，
尚未经过观测校准。首轮公开观测对照显示，本影内部明显偏亮，见下方观测对照。

## 运行

需要 [uv](https://docs.astral.sh/uv/) 和 Python 3.12。依赖由 `uv.lock` 固定。

```bash
uv sync
uv run lunar-eclipse simulate --out outputs/baseline
```

默认计算 380–780 nm（5 nm 间隔）、距影轴 0–9000 km（20 km 间隔），
采用 300 DU 臭氧和 550 nm 参考折射路径。太阳侧的光线路径只预计算一次，
各接收位置使用表格插值与有限太阳圆盘积分。计算光谱通常只需几秒，出图另计。

默认月球中心距影轴 2860 km，处在接近本影边缘的全食位置。用同一组物理参数改变月球位置：

```bash
# 更深入本影的月全食
uv run lunar-eclipse simulate --moon-offset-km 0 --exposure-stops 12 --out outputs/central

# 月球跨越本影边缘的偏食示意
uv run lunar-eclipse simulate --moon-offset-km 4200 --out outputs/partial

# 改变臭氧柱量，或提高数值分辨率
uv run lunar-eclipse simulate --ozone-du 400 --height-step-m 50 --wavelength-step-nm 2.5 --out outputs/refined
```

默认输出：

| 文件 | 含义 |
|---|---|
| `profile.png` / `profile.csv` | 从本影到半影的光谱透过率、相对明视觉亮度与颜色 |
| `spectra.png` / `spectra.npz` | 各位置的光谱照度，含正常阳光参考及臭氧关闭对照 |
| `rays.npz` | 掠过高度、真空冲量参数、偏折角、空气/臭氧柱量及带符号的落点 |
| `moon.png` | 共同曝光的中性月面屏示意 |
| `ozone_comparison.png` | 相同曝光下开启/关闭臭氧的对照 |
| `moon_exposures.png` | 同一光场的三档全局曝光 |
| `moon_chromaticity.png` | **逐位置等亮度的色度诊断，不代表真实亮度** |
| `metadata.json` / `report.json` | 参数、近似、数值网格和显示规则 |

预生成结果见 [figures/baseline](figures/baseline)。原始光谱照度以 W m⁻² nm⁻¹ 保存，
图像变换不会修改它。月球图当前是均匀、中性反射屏对地影截面的采样；
尚未包含真实月壤光谱、月面 BRDF、球面入射角明暗或地面观察者的大气。

## 方法与验证

[模型说明](docs/MODEL.md)记录公式、归一化和近似；[数据来源](eclipse/data/SOURCES.md)
记录从天空颜色项目复用的数据、上游来源和校验值。本项目离线运行，不依赖另一个仓库。

```bash
uv run pytest -q
uv run python -m validation.check_physics --out outputs/validation/physics.json
uv run python -m validation.check_convergence --out outputs/validation/convergence.json
```

独立的笛卡尔光线 ODE 核验路径偏折和沿程柱量；整个接收平面的能量亏损与入射端遮挡/消光
预算核对；无大气极限与几何圆盘遮挡核对。分辨率检查分别提高路径表和光谱采样。
这些检查验证声明模型的数值实现，**不等同于观测验证**。

## 观测对照

已将 2011-06-15 月食的公开 503／677 nm 观测曲线与固定首版模型比较。
在最外侧采样点 0.575°，模型约高 9%–15%；越向本影内部，偏亮增大，部分位置达到 4–7 倍。
当前还不能称为定量上足够接近，下一步优先研究缺少的低层消光。

[对照图](figures/observations2011/comparison.png)、[结论与资料目录](docs/OBSERVATIONAL_VALIDATION.md)、
[来源和提取说明](validation/observations/ugolnikov2011/README.md)。使用的是明确标记的论文曲线
数字化读数，未找到原始数值表；未知滤镜响应和逐点误差限制了判断。这组数据没有覆盖绿松石带本身。

```bash
uv run python -m validation.compare_observations --out outputs/observations2011
```

## 首版范围

包含球面标准大气折射、瑞利散射出射束造成的消光、臭氧吸收、有限距离和大小的太阳。
默认太阳亮度均匀，传播使用一阶小角近似；折射路径暂不随波长改变。
臭氧采用归一到指定 DU 的高斯高度分布，属于理想化廓线。

后续可以依据现有结果，分别研究气溶胶消光、太阳临边昏暗、折射色散、散射进入月球方向的光、
真实臭氧廓线及沿地球边缘的天气变化。当前未包含气溶胶、云或散射入射光，
因此不能据此预测某次真实月食的亮暗。
