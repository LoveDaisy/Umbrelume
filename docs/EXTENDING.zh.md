# 如何添加一个物理增量（removable component）

本文是人类视角的上手清单，沉淀自既有任务反复重复的同一模式。AI 代理视角的
硬约束见 [../AGENTS.md](../AGENTS.md)，物理契约见 [MODEL.md](MODEL.md)，
方法论叙事见 [MAINLINE.zh.md](MAINLINE.zh.md)。

English version: [EXTENDING.en.md](EXTENDING.en.md)（译本非权威，以本文件为准）

## 进入前：判别力优先

- **探针先于物理**：先用同口径反事实探针确认候选碰得到当前残差，再实施组件
  （见 MAINLINE「探针先于物理」一节）。
- **预注册**：先钉预测对象、准则与容差，再让另一条执行路径核对；生产实现
  不许追着预期结果调参。
- **证据身份三分**：独立先验 / 同事件校准 / 独立预测。新参数必须标明身份，
  校准值永不进默认（AGENTS.md 既有裁定）。

## 实施清单

1. **独立模块或独立函数**：不与既有物理项合并谱律或作用路径；可拆卸是硬要求。
2. **零值复原**：关闭新组件时逐位恢复前一状态（气溶胶模式：`aerosol=False`
   在同一采样光线上做对照；湿度 f(50%)≡1.0 位精确）。
3. **反事实测试**：组件 on/off 对照 + fail-loud 校验（参数域、来源锚定、
   结构配对不变量）。
4. **L0 独立验证闸**：`check_physics` / `check_aerosol` / 受影响的
   `check_convergence` 与 committed 快照比对——纯等价改动应 BYTE-IDENTICAL，
   物理改动则显式更新快照并说明。
5. **CLI 接线**：参数校验 fail-loud；metadata 写 provenance（来源、身份、哈希）。
6. **数值口径单源**：默认值提为模块常量，签名与 CLI metadata 引用同一常量
   （QUADRATURE_ORDER 模式，2026-11 语义收编确立）。
7. **文档联动**：MODEL.md（契约）、physics/AEROSOL.md（公式与来源）、
   ROADMAP 行（里程碑或进入条件）、SOURCES.md（若引入数据）。
8. **快照纪律**：活性比对快照随输入指纹刷新（`figures/mainline` 由
   `tests/test_mainline_figures` 机械强制）；冻结档案保持产出时原貌。
9. **双语同步**：README 与 ROADMAP 双语对需同步更新；文档改动跑 pytest
   即触发链接闸与字段闸。

## 显示与绘图规则

- matplotlib 基础风格单源：`eclipse/style.py` 的 `BASE_RC`（全库共享 7 键）与
  `HOUSE_STYLE`（report 专属扩展）。新绘图脚本
  `from eclipse.style import BASE_RC`（或 `HOUSE_STYLE`），额外语义键在
  `BASE_RC` 之上局部叠加，**不得重定义共享键**。
- `tools/render_mainline.py` 的 `PLOT["font"]` 与 `BASE_RC["font.family"]`
  同值双写（provenance 元数据复制）；改字体需两处同步。
- 图内关键警告就地重复（如 equal-Y 色条「只是色度诊断，不是亮度预测」）——
  语义单一权威、多处提醒，属允许的重复。

## 验证闸速查

```bash
uv run pytest -q                      # 60s 预算；独立数值交叉检查最后砍
uv run python -m validation.check_physics --out outputs/validation/physics.json
uv run python -m validation.check_aerosol --out outputs/validation/aerosol.json
uv run python -m validation.check_convergence --out outputs/validation/convergence.json
```

观测对比脚本是 validation 产物，不入 pytest 门；重跑后按快照纪律比对。
