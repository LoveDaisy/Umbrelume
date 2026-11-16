# 浏览器月食 demo

[English translation](WEB_DEMO.en.md) · **本页为权威版本**

浏览器是 Python 生产模型的一个**解读出口**，不是第二套大气或折射求解器。
它展示中性月面接收屏上的照明场，不是事件预测、精度认证或暗适应视觉模型。

## 怎样操作

- 拖动大月盘或小型影区导航器改变月心位置；聚焦后方向键移动，Shift 加速。
  导航器虚线表示几何本影与半影边界，中央十字为影轴，单位 km。
- 臭氧滑杆只改变臭氧吸收；AOD 是 50% RH 的 FT 垂直参考值。
  RH 只缩放 FT 幅度，双层先验的平流层背景独立保持 0.005。
- AOD 精确零关闭**全部**气溶胶，并恢复分子-only 网格族；此时 RH 不可调。
  回到非零态保留先前 RH。`none` 同时禁用 AOD/RH，臭氧仍可调。
- `calibrated:` 条目逐项来自生产 registry，明确标为已拟合、非物理先验；
  不接受 AOD/RH 覆盖，也不会成为默认状态。
- 曝光只改变显示；读回光谱里的 SI 值不受曝光影响。
  全场使用一个太阳白点，经 Bradford 适应至 D65；先在线性 RGB 域径向插值，
  再统一曝光、裁负、共同通道高光压缩及 sRGB 编码。

本站仅开放固定几何/标准分子大气、380–780nm/5nm、臭氧 0–600DU、
FT 参考 AOD 0–0.2（正数最小 0.001）、RH 50–99%、曝光 −4 至 +20EV。
整个月盘必须位于导出的 0–9000km 域内；不把越界值钳到最后一行颜色。
月心允许二维移动，最大距影轴为 9000km 减月球半径。

## 数据与实现边界

[导出器](../tools/export_web_demo.py) 调用生产 `build_ray_table`、
`receiver_quadrature`、光谱/色彩 API 和校准 loader，保存每个网格族的：

1. 精确去重的求积 `b` 节点、CSR 径向索引和守恒非负面积权重；全部 landing 分支保留。
2. 在 `b` 上插值后的分子柱量、单位 DU 臭氧柱量及逐层单位 AOD tau，
   每层配自己的谱因子；另含解析真空直射项，出食不会错误变黑。
3. 食外 SI 太阳谱、CIE 积分权重、统一太阳白点的线性颜色矩阵。

浏览器按 **tau 求和 → exp → 面积收缩 → SI 光谱/CIE → 统一显示** 执行。
它只插值生产导出的 RH 增长节点，不复制湿度物理、截面公式或折射积分器。
同族 AOD/DU 不改变 W；不同活跃层支持可能有不同节点，不冒称跨 preset 同网格。
精确零切换到 `none`，臭氧 off 不改变折射。

`manifest.json` 给出 schema/软件版本、输入文件 SHA256、规范生成命令、
单位/维序、dtype/小端序、shape/offset/byteLength、gzip 文件校验和、
网格族/preset/层配对、数值参数、几何、默认值与支持域。
只有整组数组字节完全相同才共享网格族。数据用显式 gzip，解码不依赖服务器压缩设置；
缺失资产、未知 schema、损坏 hash、错误 shape 均显示错误并停止，不回退到旧场。

初始加载显示状态，所有 preset 预载、解码和验证完成后才开放操作。
实际规模允许所有精确共享网格族一次解码，随后释放 gzip 副本；避免切换时反复解压。
Worker 队列最多一项运行、一项最新等待，旧结果不能覆盖新参数。
画布固定 560×560 backing store，不按手机 DPR 无界扩张。

## 本地复现与机械闸

从对应 worktree 根目录执行，输出目录必须不存在或为空：

```bash
uv sync --locked --group web
uv run --group web python -m playwright install chromium
uv run python tools/export_web_demo.py --out outputs/web-demo-a
uv run python tools/export_web_demo.py --out outputs/web-demo-b
diff -qr outputs/web-demo-a outputs/web-demo-b
uv run --group web python -m validation.check_web_demo --site outputs/web-demo-a --out outputs/validation/web-demo.json
uv run pytest -q
```

本地浏览可运行 `uv run python -m http.server 8000 --directory outputs`，
打开 `http://localhost:8000/web-demo-a/`。不能用 `file://` 加载 Worker/模块。
没有 npm、后端、Pyodide 或 GPU 必需路径。Linux CI 安装浏览器时加 `--with-deps`。

[独立交叉闸](../validation/check_web_demo.py) 从当前生产 API 新建参考，
启动项目子路径 HTTP，实际运行 app→Worker→Canvas，而非只测一个未被 UI 调用的 helper。
它覆盖参数边界/中间值/固定种子随机值、所有径向节点与波长、真实 CLI CSV/NPZ、
各径向 cell 中点和几何边界附近、像素 mask、快速交错输入、无效参数及故障注入。
几何半影边界超出合法月心域的像素用月盘外缘触及边界来覆盖，不越界钳值。

固定出口容差（逐点，不以平均值替代）：

| 量 | 判定 |
|---|---|
| 同节点 T、relativeY、XYZ/太阳参考Y | atol=2e-12，rtol=2e-8 |
| 同节点未裁负线性 RGB | atol=2e-11，rtol=2e-8 |
| SI 食外参考谱 | atol=1e-14 W m⁻² nm⁻¹，rtol=2e-12 |
| 离节点 relativeY | atol=5e-9，rtol=0.005 |
| 离节点显示 RGB，EV=−4/0/8/12/20 | 每通道绝对误差≤3/255 |
| 同口径 Canvas | 每通道≤1个 8-bit 码值，另验 mask |

报告保存逐量最大误差/最差坐标、参数、指纹、全部计时样本、缓冲清单与桌面/窄屏截图。
故意篡改 W、气溶胶谱因子或颜色矩阵，即使更新 hash，也必须使数值闸失败。
缺浏览器、缺资产、零样本等非零退出，不能 skip 变绿。

性能口径：预载后30次热身、200次预定参数/preset 变更，
从应用接收参数到完整场+CIE+560×560绘制完成，p95和所有样本均≤100ms；
拖拽/曝光 p95≤33ms。全部支持 preset 的发布数据≤32MiB；
明确列出的数值/图像/通信缓冲保守峰值≤192MiB，不把它当整个浏览器 RSS。
手机 viewport 截图只证明布局，不是真机性能证明。
这些是出口复现/响应门槛，不是物理精度认证。

## CI 与发布

[Web workflow](../.github/workflows/web-demo.yml) 对所有 PR、main push 构建，
不使用可能遗漏核心改动的 path filter；重新生成资产与 fresh 参考，再执行 Python、
物理/气溶胶/收敛和完整浏览器闸。失败不得发布。
普通 pytest 不要求安装浏览器。

`v*` tag（必须在 main 历史内）或 main 的显式 workflow_dispatch 才进入部署。
同一 run 保存已验证原始站点、验证报告、源码 SHA 和 manifest hash；部署 job
下载并复核同一 artifact，不重建。权限限于构建 `contents:read`，
部署另加 `pages:write`/`id-token:write` 与 `github-pages` environment/concurrency。
首次部署需本任务合入 main 后由授权 owner 触发；本实现不会自行 merge、push 或打 tag。
Pages 设置/环境保护需要额外授权时上抛，不绕过。

**运行验收与上线状态以实际 Actions/部署记录为准，不能由本文推定已上线。**
真机操作/性能和人工审美须独立验收，本地绿色不能替代。

## 来源与限制

[逐源来源与条款](../eclipse/data/SOURCES.md) 随站点提供副本，包含部分核实状态。
代码 MIT 不等于所有第三方数据 MIT。只发布生产数值衍生表，
不打包原论文、观测下载或 wheel；不隐式下载任何数据。

页面保留三条边界：浏览器为解读出口之一，核心求解以仓 CI 为准；
现象呈现，非精度认证/预测；散射入射光等现象模型本身不含。
中性月面屏不含月壤纹理/BRDF、观察者大气或暗视觉；高光压缩后的像素不是照度计。
权威物理与出口登记见 [MODEL](MODEL.md)。
