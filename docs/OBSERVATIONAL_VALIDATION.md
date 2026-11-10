# 月食观测参照与首轮对照

本文记录 M2 的首版分子模型对照，保留当时的参数和结论。本轮气溶胶改进及新默认值见
[ROUND2.md](ROUND2.md)，整体路线见 [ROADMAP.md](ROADMAP.md)。下列历史数字未被新结果覆盖。
后续资料扩充与观测口径复核见 [MULTI_EVENT_OBSERVATIONS.md](MULTI_EVENT_OBSERVATIONS.md)。

## 首版分子模型结论

第一版已通过数值闭合检查，但与真实月食相比，**本影内部明显偏亮，目前不能称为已经足够接近**。
首次对照使用 [Ugolnikov 等（2011）Figure 2](https://arxiv.org/pdf/1106.6178v2)
的西侧本影测光曲线。下表是首版分子模型预测／论文观测曲线读数；1 表示一致。

| 距影心角距离 | 503 nm | 677 nm |
|---|---:|---:|
| 0.150° | 没有该处观测曲线 | 7.10 |
| 0.400° | 4.07 | 2.30 |
| 0.500° | 1.56 | 1.56 |
| 0.575° | 1.15 | 1.09 |

这是明确的“越往内部，偏亮越多”的结构。0.575° 是本次读出的最外侧点，
仍在本影内，不能据此说本影边缘或绿松石带已被验证。

[对照图](../figures/observations2011/comparison.png)展示原始亮度比和模型／观测比，
[完整结果](../validation/results/observations2011/comparison.json)及 CSV 保留逐点数值。
该轮没有拟合亮度比例、位置偏移或任何物理参数，生产模型当时未改动。

## 口径与证据强度

观测量为同一月面单元在食中／完全食外的亮度比，可与模型的 E_lambda/F_lambda 对照。
从固定版本 PDF 的实线矢量轮廓提取了十八个曲线读数，保留 SHA-256、坐标、路径编号和提取脚本；
**这些是论文图表数字化数据，不是作者发布的原始数值表，也不是十八次独立曝光**。
曲线选择及若干坐标经另一条渲染图读取路径复核。

论文没有提供滤镜响应、FWHM 和逐点误差。因此当前用有效波长的单色预测作近似；
图表读取误差不能冒充观测误差，不能据此给出统计显著性或严格验收阈值。
细节及重建命令见[数据说明](../validation/observations/ugolnikov2011/README.md)。

采用 JPL Horizons DE441 在 2011-06-15 19:30 UTC 的日地、地月距离，代表可见光观测时段。
使用该事件距离而非首版平均距离，部分预测可改变约 27%，所以几何条件需要先匹配。
19:20–19:40 的距离变化在这些点上造成的预测变化至多约 0.115%。
没有重建逐次曝光、月面深度或站点的完整几何。

## 目前能排除与尚不能归因的部分

- 路径表间隔 100→50 米时，本次采样点的预测最大变化约 0.0122%，远小于内部数倍偏差。
- 臭氧 300→250 DU 时，最大变化约 7.55%，且预测更亮，不能解释当前主要偏差。
  250 DU 来自论文的参考模型，不是从这组数据拟合出来的值。
- 在 0.50° 和 0.55°，两个波段的传输比之比（模型／观测）分别约为 1.00 和 1.01，
  但到 0.40° 变为约 1.77。偏差不仅是全局亮度倍率，内部光谱形状也需要解释。

优先研究的假设是：**目前缺少沿较低掠过路径的额外消光**，气溶胶是一个有物理依据的候选。
已有观测论文也把空间分布与气溶胶联系起来，但这些结果尚不能唯一识别根因。
实际温度／臭氧廓线、云层和太阳临边昏暗均未在首版中匹配；滤镜积分也是未闭合口径
（已由下节带通口径增补量化，影响为次要量级）。

该轮据此提出：采用固定观测基准，先加入有来源依据的简单气溶胶高度分布与光谱消光，
检查能否同时降低内部偏亮并保持较外侧的预测；独立评估太阳临边昏暗和大气廓线敏感性。
散射进入月球方向的光增加接收照明，仅靠增加这项不能消除当前偏亮，
应在消光、几何和光源条件厘清后再判断它在哪些波段/位置重要。
每项增量都保留同一组未经重新归一化的对照，不以想要的颜色作为拟合目标。

## 滤镜带通口径（2026-11 增补）

`ROADMAP.md` 登记的「补齐滤镜与观测口径」中，滤镜响应部分已按下述方式闭合：把 503/677 nm
的单色点预测升级为滤镜带通预测，量化带通口径对模型/观测比残差的贡献。
**不改物理模型与默认参数，不引入任何拟合自由度**（`fitted_parameters` 恒为空表）。

### 滤镜信息来源与近似等级

对固定版本 PDF（SHA-256 绑定）重读裁定：论文只说明 SAO 使用有效波长 503/677 nm 的窄带滤光片，
未给出 FWHM、带宽、带型或响应曲线（24 个参数类关键词全文检索零命中，检索记录与逐字引文落盘于
`source.json` 的 `filter_information`，可用提取工具独立复跑）。近似等级的**唯一权威**是
`source.json` → `filter_information.approximation_grade`（假设矩形与高斯两种带型、FWHM
{10, 20, 40} nm 探针包络；宽度与带型未经实测）。本文与其他文档只引用不复述该字段。
503 nm 滤光片的弱红尾已在论文归约中用 677/867 nm 数据改正，无需建模红漏。
「向作者索取实测响应曲线」仅记录为升级点，未发生任何外部通信。

### 带通可观测与对照结果

观测量是同一月面单元过滤光片的食中/食外流量比，食外参考自带 R(λ)·S(λ) 权重，故主定义取
`T_band(θ) = ∫R·S·T dλ / ∫R·S dλ`（T 为模型 `transmission_profile` 输出，S 为 AM0 太阳谱；
月面反照率在窄带内按灰体处理——与单色对照同一假设等级）。「仅 R 权重」敏感变体单独落盘。

逐点对照（[完整结果](../validation/results/filter_response_2011/analysis.json)、
[CSV](../validation/results/filter_response_2011/comparison.csv)、
[对照图](../figures/filter_response_2011/comparison.png)）：

- 主变体（固定 OPAC 先验）带通与单色的模型/观测比之差：503 nm 全包络最大 0.108
  （在 503 弧段内端 0.400°，该点单色残差 1.99）；677 nm 最大 0.064（在内侧 0.150°，残差 2.76）。
  分子 v1 变体同量级（含两变体全包络 0.147 / 0.076）。
- 最外点（0.575°）：主变体下带通把 503 比值上推 +0.0057 至 +0.0337（主变体全带型×宽度组合，
  字段 `main_variant_signed_range`；含分子 v1 全变体至 +0.036）、677 下推 −0.0015 至 −0.0127
  （全变体至 −0.0145）。方向（主变体口径）：两波段位移均推离 1，不减小最外点残差（主变体口径
  偏差比 1.06–1.07，字段 `outermost_deviation_ratio_best_main_variant`）；分子 v1 在 677 nm
  最外点位移朝向 1（全变体口径最优偏差比 0.84，字段 `outermost_deviation_ratio_best_all_variants`，
  分类仍为「仅边际改变」）。区间由脚本确定性生成于 `analysis.json` →
  `per_band.<band>.outermost_band_minus_mono`，本文只引用不复算。
- 结论（双方向都合法，本组数据落在分层表述）：
  - **最深点**（677 nm 0.150°）：单色残差 2.76，带通效应上界 0.064 且方向推离 1——带通不能解释。
  - **503 nm 全段与 677 nm 最外点（主变体口径）**：带通位移全部推离 1，对残差的减小量为零
    （503 nm 的全变体最大缩减亦为负值，字段 `outer_toward_unity_reduction.all_variants`）；
    最外点带通效应可达该点残差的 41.6%（503）与 48.3%（677）（主变体最大 |位移|/|偏差|，字段
    `outermost_band_minus_mono.main_variant_max_abs_shift_over_deviation`）。分子 v1 在
    677 nm 最外点朝向 1（偏差比 0.84，见最外点条目）。
  - **677 nm 外区中间点**（0.400–0.550°）：存在朝向 1 的位移——主变体最大缩减 0.0131
    （0.500°，最宽探针）、占比最高 21.7%（0.550°，最宽探针；名义 20 nm 宽度下 12.0%），字段
    `outer_toward_unity_reduction.main_variant`。量级不足以改变外区残差结构：该弧段单色残差
    仍为 1.06–1.68 倍的系统性偏亮。
  - 总体：带通口径对这组残差结构的贡献是次要量级（全变体包络最大 0.147，对比残差结构
    0.026–2.76），不改变「越往内部偏亮越多」的结构判断；单色有效波长近似在当前证据精度下
    保留；滤镜带通由「未闭合口径」降级为「已量化、影响次要」。

### 近似充分性（H2）与定义敏感度

- 同宽度带型效应（矩形 vs 高斯）随宽度分解（字段 `shape_vs_width.shape_effect_by_width`）：
  10 nm 档 0.015、20 nm 档 0.010、40 nm 档 0.048（同宽度总最大 0.049），宽度展开（10→40 nm）
  效应最大 0.076，比值 0.64：包络宽端（40 nm）带型不可忽略，窄端效应约为宽度效应的 1/5–1/7。
  「矩形/高斯近似够用」仅在真实 FWHM 落在窄端时成立——这是索取实测响应曲线升级点的量化依据。
- 食外参考权重选择（R·S vs 仅 R）造成的比值差（字段 `definition_sensitivity_R_weight_only`）：
  503 nm 全部 8 个读数（按钉死边界 θ ≥ 0.400° 均属外侧区域）为 +0.0314 至 +0.0929，最大值在
  503 弧段内端 0.400°；677 nm 外侧 ≤0.0051。敏感行已落盘，定义选择的影响量级与带型效应
  相当、仍小于同点残差。
- 积分网格 1.0→0.5 nm 加密的最大相对变化 9.1e-4（< 1e-3 阈；0.5 nm 已低于 AM0 表原生分辨率，
  更细加密在线性插值框架下不再改善，终止原因落盘 `analysis.json`）。

### 判读阈值审计

脚本生成的结论字段按确定性规则分类，以下规划判断阈值钉死在
`validation/filter_response_2011.py` 常量区供独立审计：外侧区域边界 θ ≥ 0.400°；
「实质性移动」阈（带通效应 ≥ 0.5 × 同区单色残差）；最外点分类阈（best-variant 偏差比
< 2/3 为实质性减小、> 1.5 为不减小）；带型/宽度比阈 0.1。单色列与两份已提交结果
逐点交叉核对 rtol = 0（逐位一致）后方产出任何带通产物；复核命令见该文件 docstring。

## 跨事件消光缺口判别（2026-11 增补）

四个独立口径（June 2011 曲线读数、1990 整月面、Dec 2011 七点谱、HST 2019 本影谱）
的「贡献路径高度分布」逐一算出，与已提交残差做跨事件相关检验，并对
「完整背景地球大气」缺失消光成分做文献先验反事实探针。本节只做相容性排序，
不做唯一因果归因；下一个物理增量的裁定属 owner。全部产物在
`validation/results/extinction_gap_diagnosis/`（`height_stats.json`、
`correlation.json`、`discrimination.json`、`sources.json`、`prior_parameters.json`），
由 `uv run python -m validation.extinction_gap_diagnosis --out
validation/results/extinction_gap_diagnosis` 离线再生（byte-identical；
sources / prior_parameters 为一次性文献产物，命令侧仅做 schema 校验）。

### 机器与口径

- 「贡献核」= 生成该口径已提交 prediction 的同一加权混合核按切点高度 b 的分解，
  逐项镜像 `eclipse.transport.transmission_profile` 的 tau 构造；零阶矩闸
  Σw + 真空外枝 = 已提交 prediction 以 1e-12 相对容差在 pytest 全口径全变体通过
  （`height_stats.json:summary.record_count` 条记录，最差偏差
  `height_stats.json:summary.worst_relative_gap`；容差语义 = 结构复现使相等由构造
  成立，预算仅覆盖浮点求和序，见
  `height_stats.json:kernel_identity.relative_tolerance`）。
- 独立锚：June 2011 分子套与 round2 已提交
  `molecular_tangent_height_p10_p50_p90_m` 逐点精确相等（pytest 常驻）。
  HST 形状基线同样锚定：已提交 `slope_model_vs_observed` 的带中位聚合作用在
  **全窗口网格**上（430 nm 带含 [425, 430) 行），探针形状分支复用同一聚合结构，
  基线由 pytest 逐值相等断言与 AOD=0 斜率恒等闸双重钉住——探针不引入第二种
  同名口径。本机器是四口径高度统计的唯一权威实现；round2 的
  `altitude_quantiles` 为历史产物生成器，不再扩展。
- 高度带 <2 km / 2–12 km / >12 km 与模型组件边界一一对应
  （`height_stats.json:kernel_identity.band_semantics`）；真空外枝（切点高于
  大气表顶）实测在四口径全部几何关闭，记录值为浮点噪声量级（≤1e-9 噪声底由
  pytest 钉住），不进入大气枝统计。
- 残差一律读已提交产物，不重算模型；亮度轴与形状轴字段分离
  （`correlation.json:scope`）。

### 跨事件相关（亮度轴）

事件级残差 = 各口径定量读数 ln(model/observed) 中位数（HST 与已提交
`model_to_observed_median` 逐值相等抽核，见
`correlation.json:hst_consistency_check`）。固定 OPAC 变体四口径残差对数中位
0.490 / 1.089 / 0.502 / 0.451（June / HST / Dec / 1990 顺序，下同；各
`correlation.json:event_level.june_2011.fixed_opac.residual_ln_ratio_median`
等四个 `correlation.json:event_level.*.fixed_opac.residual_ln_ratio_median`
分支），对应贡献核 P50 中位高度 10.6 / 8.9 / 12.8 / 11.9 km
（`correlation.json:event_level.*.fixed_opac.p50_height_m_median`）。

深度方向：分子变体 ρ(P50 中位, 残差) = −1.0、六对「更深者残差更大」全部成立
（`correlation.json:cross_event.molecular.spearman_rho_residual_ln_vs_p50_median`
与 `...pairwise_direction.deeper_has_larger_residual_count`）；固定 OPAC 变体
ρ = −0.40、四/六对。「贡献路径越深、模型越偏亮」的深度格局在分子变体上跨事件
完美秩一致，在已含自由对流层气溶胶的固定 OPAC 变体上减弱——气溶胶项吸收了
部分深度信号。N=4，只报排序一致性，不做显著性声明。

### 判别探针与相容性排序

12 个缺失消光候选（OPAC 1998 Tables 3/5 十个边界层类型 @RH80%、mineral
transported 2–3.5 km、平流层背景 12–35 km）以文献 AOD550 构造反事实探针，
经同一 kernel 适配器附加在固定 OPAC 基线上（AOD=0 严格复原基线，pytest）；
先验逐参数钉 SHA+页码（`prior_parameters.json`，OPAC PDF SHA 与 round2 已钉值
逐字节一致）。各候选与基线的叠加语义逐行钉在产物中
（`discrimination.json:candidates[antarctic].stacking_semantics`；
antarctic 行另附双重计入的代理解读字段 `...candidates[antarctic].interpretation`）。
判定带（`discrimination.json:verdict_bands`）：先验量级效应闭合
残差对数缺口的 [0.5, 2] 为相容、[0.1, 0.5)∪(2, 10] 为张力、其余为排除倾向；
机械排序键 = 四口径亮度闭合的最小值（最弱链）降序
（`discrimination.json:ranking`）。

主要结果（相容性排序，不是归因）：

1. **0–2 km 边界层类型（continental/maritime/urban/arctic 各湿度型）在四口径
   事件级全部排除倾向**：贡献核对 <2 km 路径的质量占比在 HST/Dec ≈ 0、June
   中位 = 0、1990 仅 0.55%（`correlation.json:event_level.*.fixed_opac.band_fraction_below_2km_median`），
   先验 AOD（如 continental clean 的 0.046）推不动任何事件级残差。它们只触
   June 最深读数——最深单点闭合 20–25%
   （`discrimination.json:candidates[continental_clean].events.june_2011.brightness.max_readout_gap_closure_fraction`
   = 0.229），事件中位不动。
2. **平流层背景（12–35 km，τ=0.005，OPAC 同源先验）四口径亮度全部张力**
   （闭合 0.19–0.44，`discrimination.json:candidates[stratosphere_background].events.*.brightness.gap_closure_fraction`），
   但是**唯一改善 Dec 蓝/红形状的候选**（450/750 向 unity 收敛分数 +0.42、
   500/750 +0.57，
   `discrimination.json:candidates[stratosphere_background].events.dec_2011.shape.450_to750_transmission_ratio_model_over_observed.disagreement_closed_fraction`
   与 `...500_to750_transmission_ratio_model_over_observed.disagreement_closed_fraction`），
   June 形状亦相容（+0.79，
   `discrimination.json:candidates[stratosphere_background].events.june_2011.shape.disagreement_closed_fraction`）；
   HST 斜率微恶化（−0.076，
   `discrimination.json:candidates[stratosphere_background].events.hst_2019.shape.ody902030.430_over_570.disagreement_closed_fraction`）。
3. **antarctic 型（0–10 km 层、≈5×自由对流层 AOD 的弱谱依赖配置）亮度 3/4
   相容**（June 1.13 过闭合、HST 0.53、1990 0.51、Dec 0.23，
   `discrimination.json:candidates[antarctic].events.*.brightness.gap_closure_fraction`），
   但 **Dec 形状强烈恶化**
   （450/750 收敛分数 −0.68、500/750 −0.93，
   `discrimination.json:candidates[antarctic].events.dec_2011.shape.450_to750_transmission_ratio_model_over_observed.disagreement_closed_fraction`
   与 `...500_to750_transmission_ratio_model_over_observed.disagreement_closed_fraction`）。
   该行实质是「2–12 km 更多同类消光」假设的 OPAC 代理解读——与 677 nm 同事件
   校准方向一致，但形状代价与平流层候选相反。湿度化自由对流层本身因先验不可得
   未探针（检索过程见 `prior_parameters.json` 的 unavailable 行）。
4. desert（0–6 km 层）仅 June 相容（0.53）；mineral transported 事件级排除
   （负 Ångström 谱经 validation 侧参数化评估，
   `discrimination.json:scope`）。

### 机制发现：形状响应由高度选择性决定

「低层消光增量会恶化蓝/红形状指标」这一待核查假设获定量证实，且机制超出
Ångström 指数的直觉：形状变化由**贡献核高度分布 × 成分作用带**决定，而非谱斜率
单独决定。Dec 的 450 nm 核在 2–12 km 带的质量仅 9.0%、750 nm 核达 49.2%
（`height_stats.json:records.41.band_fraction_2_to_12km` 与
`height_stats.json:records.53.band_fraction_2_to_12km`，fixed_opac 变体的
450/750 nm 两点），因此作用于低层的候选（含近中性谱的
antarctic 型）打在红核占优的 2–12 km 质量上，f(450)/f(750) > 1，把本已偏蓝的 Dec 比值推得
更远；作用于高层的平流层候选命中 >12 km 质量占优的蓝色核，反而改善。任何「下一个物理增量」的
形状轴代价都必须按此机制用 kernel 评估，不能按谱斜率手推。

## 其他已核实的公开参照

| 参照 | 公开程度与可比较内容 | 当前限制 |
|---|---|---|
| [Youngblood 等（2020），HST/STIS 2019 月食](https://arxiv.org/html/2008.01837) | MAST 项目 15674；本影、半影和正常月面 FITS 可下载；G430L 约 290–570 nm | 需重新提取扩展源光谱、恢复实际月面采样位置并区分临边昏暗修正；不能直接把表中月球中心角当作采样点角度 |
| [Strassmeier 等（2020）](https://doi.org/10.1051/0004-6361/201936091) | CDS 目录 J/A+A/635/A156；高分辨率光谱 | 后续复核确认发布谱经过连续谱归一化，不能验证宽带照度／颜色；波段交叠本身不构成可比性 |
| [Pallé 等（2009），2008 月食](https://doi.org/10.1038/nature08050) | 发表了可见/近红外透射光谱 | 本轮未找到可直接使用的处理后数值表；该次受火山气溶胶影响，更适合后续复合大气测试 |

HST 归档下载入口已通过实际 HTTP 请求确认（本轮没有据此宣称完成 HST 光谱比较）：

- [本影 ody902030_flt.fits](https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/ody902030_flt.fits)
- [半影 ody9a2010_flt.fits](https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/ody9a2010_flt.fits)
- [正常月面 ody904080_flt.fits](https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/ody904080_flt.fits)

HST 位于大气外，能避开地面观察站的大气传输；但月面指向漂移和反照率仍需处理。
该数据不覆盖完整约 600 nm 吸收谷。验证绿松石带还需要覆盖相关波长、接近本影边缘，
并有位置和仪器响应记录的观测。普通经过白平衡、伽马和局部增强的照片只适合作视觉参考。

## 复现

```bash
uv sync --locked
uv run python -m validation.compare_observations --out outputs/observations2011
uv run python -m validation.filter_response_2011 --out outputs/filter_response_2011
```

此命令离线运行，用已提交的图表读数、事件距离和现有首版物理模型生成图、CSV 与 JSON。
滤镜带通对照同样离线运行；其单色列先与上面已提交结果逐点交叉核对后才产出带通产物。
