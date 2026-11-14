# docs 索引 / Docs index

文档分两层：面向读者的说明放本层，内部开发记录下沉 [`records/`](records/)。
语言约定：**文件名后缀只标语言（`.zh.md` / `.en.md`），权威版本逐对在下方声明**；
契约与记录类文档保持单语（豁免理由见下）。双语内容同步靠定期人工检查；
pairing 机械闸（`tests/test_mainline_docs.py` 的 `test_bilingual_pairs_exist_and_cross_link`）
保证集合内每对文件存在、互链且全部本地链接可解析。

## 面向读者 / Reader-facing

| 文件 | 语言 | 权威 | 职责 |
|---|---|---|---|
| [../README.md](../README.md) | en | 英文 canonical | 项目入口：是什么、怎么跑、结果与边界 |
| [../README.zh.md](../README.zh.md) | zh | 与英文版互为镜像 | 同上（中文） |
| [MODEL.md](MODEL.md) | en | 单语（物理契约不设译本：一词两译必分叉） | 物理与数值契约、限制 |
| [MAINLINE.zh.md](MAINLINE.zh.md) | zh | 单语权威 | 主线叙事：阶梯抽象、方法论五支柱、模块地图、叙事钩子（英文译本在后续 chore 队列） |
| [ROADMAP.zh.md](ROADMAP.zh.md)（权威）/ [ROADMAP.en.md](ROADMAP.en.md)（译本） | zh / en | zh | 里程碑、上一轮物理选择、后续方向与进入条件 |
| [physics/AEROSOL.md](physics/AEROSOL.md) | en | 单语 | 气溶胶组件公式、来源与 standing gap |

## 内部开发记录 / Internal development records

`records/` 是过程证据与历史裁定（append-only 证据簿），不面向普通读者，
**保持产出时语言、不做翻译**——翻译会制造持续付息的同步负债，而读者价值趋近于零：

| 文件 | 职责 |
|---|---|
| [records/OBSERVATIONAL_VALIDATION.md](records/OBSERVATIONAL_VALIDATION.md) | 首轮观测对照、口径、判别证据与机制发现 |
| [records/ROUND2.md](records/ROUND2.md) | 气溶胶选择与证据 |
| [records/DECEMBER2011.md](records/DECEMBER2011.md) | 2011-12 事件完整对照记录 |
| [records/MULTI_EVENT_OBSERVATIONS.md](records/MULTI_EVENT_OBSERVATIONS.md) | 多事件资料状态、使用条件与比较 |
| [records/PRIOR_SOURCE_SURVEY.md](records/PRIOR_SOURCE_SURVEY.md) | 独立先验来源普查 |

数据来源见 [../eclipse/data/SOURCES.md](../eclipse/data/SOURCES.md)（en，单语）。

建议阅读顺序：README → MODEL → MAINLINE → ROADMAP；开发记录按需从 MAINLINE
的钩子进入，不要求通读。
