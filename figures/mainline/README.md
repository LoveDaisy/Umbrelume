# 主线插图：同一生产链、同一显示口径

本目录是 [MAINLINE](../../docs/MAINLINE.zh.md) 的显示 lane 图组，不是观测拟合或新物理证据。
唯一 provenance 为 [metadata.json](metadata.json)：逐态完整命令、未指定 override 的身份、
本次 CLI 的实际模型/几何/数值状态、报告显示参数、输入源码/数据/锁文件指纹、软件版本及逐图 SHA256 均在其中。
不要用旧 `round2/render` 的单层元数据解释这组当前生产图。

## 怎样读图

- 三张月盘分别是分子 `none`、自由对流层 `opac-ft`、当前双层背景 `opac`。
  它们使用同一几何、同一尺寸和**全局 +12 EV**，没有逐态、逐半径提亮。
  +12 EV 是本组显示选择，**不是 CLI 默认曝光的新定义**，也不是相机曝光或观测亮度标定。
- 白点是食外太阳光谱经 Bradford 适应至 D65；场景的亮度参考是食外太阳光的 CIE Y = 1。
  超出线性显示范围时各 RGB 通道按共同最大通道压缩，负线性 sRGB 裁至零。
  因此显示像素不是照度计，不应从已压缩的高光直接量出物理亮度比。
- 月盘是影平面照明场的**中性月面屏**采样，没有月海贴图、月壤反射谱、BRDF、入射角阴影或观察者大气。
  色度为明视觉 CIE 1931 近似，不预测暗适应观感。
- `shadow_profile.png` 直接复用当前双层态生产剖面：上方为相对食外太阳的明视觉照度，
  前两条色带为同曝光臭氧 on/off；**ozone-off 仅移除吸收，不改变折射或其他介质**。
  最后一条 `equal Y*` 色带在每个半径移除了亮度，只是**色度诊断，不是亮度预测**；图内也有此声明。
- `transmission_ladder.png` 是三态 CLI 原始 `spectra.npz` 的 `E_lambda / F_lambda`，
  非峰值归一、非等 Y；三个分面共用绝对对数纵轴与 `T=1` 食外参考。
  半径 0、4000、4600 km 是显示采样而非事件月面站点，最后一个按本次几何已在半影。
  三态用固定颜色和实线/虚线/点线编码；完整数表位于 metadata 的逐态 `spectral_table`，不依赖辨色读数。

## 图与来源

| 发布图 | 生产来源 | 用途 |
|---|---|---|
| [l1_molecular.png](l1_molecular.png) | `none` 的 `moon.png` | 分子机制基线 |
| [l2_free_troposphere.png](l2_free_troposphere.png) | `opac-ft` 的 `moon.png` | 与 L1 同曝光对照 |
| [l25_stratosphere.png](l25_stratosphere.png) | `opac` 的 `moon.png` | 当前双层先验，不含校准增量 |
| [transmission_ladder.png](transmission_ladder.png) | 三态 `spectra.npz` | 分清谱形和绝对幅度 |
| [shadow_profile.png](shadow_profile.png) | `opac` 的 `profile.png` | 本影到半影、臭氧吸收反事实与 equal-Y 诊断 |

本组没有湿度拟合态或 `calibrated:` 面板。`opac-ft`/`opac` 显式采用已发表 50% RH 参考，
不传 `--aod550`，保留文献先验；它们不代表某次月食的实测天气。

## 与历史观测图的边界

MAINLINE 另外引用既有 [multiple_events](../multiple_events/comparison.png)、
[dec2011](../dec2011/comparison.png)、[filter_response_2011](../filter_response_2011/comparison.png) 快照，**没有重新生成**。
它们的固定 OPAC 是历史 v0.2 自由对流层单层，不是本组当前双层，更不是湿度 A/B 重反演：

- multiple_events：1990 整月面明视觉比较，发表 ± 值取角点范围不是置信区间；2007 外推最低值仅作背景。
- dec2011：名义 Tycho 位置的整段曝光平均与图表连续谱读数；读取裕量不冒充测量误差，位置敏感性不是置信区间。
- filter_response_2011：June 历史单层主变体的假设带型比较（其他变体见其产物）；没有实测滤镜响应，图表读数仍有来源限制。

跨代图不能当统一物理对照，也不能以图片替代拟合/留出与残差账户。standing gap 的裁决仍由原文所链证据承载。

## 离线复现与字节契约

预先备妥 `uv.lock` 对应环境；本入口不联网、不下载、不运行观测或反演。
在仓库根目录执行（输出根须不存在或为空）：

```bash
uv run --offline --no-sync python tools/render_mainline.py --out outputs/mainline-check-a
uv run --offline --no-sync python tools/render_mainline.py --out outputs/mainline-check-b
diff -qr outputs/mainline-check-a outputs/mainline-check-b
uv run --offline --no-sync pytest -q tests/test_mainline_docs.py tests/test_mainline_figures.py
```

入口通过当前解释器调用当前 worktree 的 `python -m eclipse simulate`，完整等价 CLI 命令见 metadata。
所有原始 CLI 产品仅存临时 staging；成功校验后发布上述五图和 metadata，不自动清空已有目录。
README 为手写说明，不由入口覆盖；再生成文件应逐项与本目录比较（不把 README 纳入生成文件集合）。

发布 JSON 递归移除时间/耗时，无绝对路径；PNG 无损重编码并验证尺寸/像素完全相同，
最终只保留 `IHDR`/`IDAT`/`IEND` chunk。metadata 不含自身 hash，避免自指。
**同一源码、输入、锁定依赖及受支持运行环境下逐文件 byte-identical**，包括 metadata；
不承诺字体、渲染器、平台或依赖升级后仍字节相同，软件/字体指纹用于解释此类变化。
