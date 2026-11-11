# Search log — task-prior-source-survey

统一条目格式（每条检索动作一行记录；`unavailable` 判定必带本文件条目指针）：

```
### [W<n>.<k>] <日期> <关键词/URL>
- 动作：<查询 / 下载 / 访问>
- 结果：<命中 / 死链 / 不含所需字段 / 超预算放弃>
- 证据：<sha256 + 缓存文件名，或 HTTP 状态码 / 页面要点>
- 放弃理由（如放弃）：<为何不再追>
```

预算闸（plan 默认假设 4）：每个维度候选来源 ≤3 即判定，放弃理由一律记录在案。

---

<!-- 检索记录自下方追加 -->

### [W1.1] 2026-11-11 "Koepke Hess 1997 Global Aerosol Data Set GADS MPI report 243 PDF"
- 动作：WebSearch 找 GADS 报告（MPI Report No. 243）原始 PDF
- 结果：命中引用与镜像线索（GitHub mrke/Global-Aerosol-Data-Set 网格化 NetCDF 重建、多篇引用论文），无报告原始 PDF 直链
- 证据：MPI 官网 URL 线索（mpimet.mpg.de 报告库）；GADS 数据文件（winter.dat/summer.dat）同样在已死 FTP 上，wayback 未存档（见 W1.3）
- 放弃理由：GADS 报告对湿度轴无增量（其组件光学数据与 OPAC 同源，BAMS98 + 组件文件已覆盖）；全球分布图对四事件纬度带类比价值边际；预算闸下先攻更高优先维度。留 1 次后续检索额度（M1 收口后视 Step 5/6 需要决定是否使用）

### [W1.2] 2026-11-11 OPAC 官方在线文档（活链）
- 动作：抓取 https://opac.userweb.mwn.de/radaer/opac.html 与 opac-des.html（官方镜像，LMU）
- 结果：页面存活；opac-des.html 明确数据文件命名（4 字母组件码 + RH 档）、8 个 RH 类（0/50/70/80/90/95/98/99%）、optdat/ 子目录、61 波长
- 证据：页面文本已核（本条为在线文档，未入 sources.json——无可钉字节产物；其内容被 W1.6 的数据文件本体印证）
- 备注：官方页描述的 OPAC31 tar 包指向 ftp.lrz-muenchen.de（见 W1.3）

### [W1.3] 2026-11-11 死链与存档缺口核验（LMU 原始 URL + LRZ FTP）
- 动作：wayback CDX 查 www.meteo.physik.uni-muenchen.de/strahlung/aerosol/*、ftp.lrz-muenchen.de/pub/science/meteorology/aerosol/*、www.lrz-muenchen.de/~uh234an/*
- 结果：1999 年 LMU 主页快照存在（重定向到 lrz-muenchen.de/~uh234an/www/radaer/，gads.html/opac.html 页面均已存档）；**FTP 数据文件本体零存档**（gads22a.tar.gz、opac31 包、全部 optdat 文件在 LRZ FTP 路径下均无 capture）
- 证据：CDX 逐路径查询输出（404/301 或空集）；~uh234an 个人主页路径下仅 HTML 页 + opac.inp/opac.out 被存档
- 放弃理由：该路径下数据不可恢复，转 W1.4 GEISA 分发路径

### [W1.4] 2026-11-11 GEISA/AERIS 门户
- 动作：抓取 https://geisa.aeris-data.fr/opac（WebFetch，JS 渲染页）
- 结果：门户明示承载经典 OPAC（BAMS98）组件数据文件直链（aeris-geisa.ipsl.fr/geisa_files/2011/... 与 cds-espri.ipsl.upmc.fr/espri/pubipsl/...）；另有独立 OPAC 2015 条目（未使用）
- 证据：完整 URL 清单已提取；门户页不可 curl 直取（JS 渲染、0 字节）
- 备注：与 prior_parameters.json search_notes 中「GEISA hosts OPAC 2015, a different dataset version」的前提冲突——门户同时承载经典版与 2015 版；本 survey 只用经版本核验的经典版文件（版本核验见 sources.json `opac1998_component_data_files.version_verification`）

### [W1.5] 2026-11-11 ESPRI 直接下载尝试
- 动作：curl 直连 aeris-geisa.ipsl.fr 与 cds-espri.ipsl.upmc.fr（多 TLS 参数、http/https）
- 结果：全部失败——TCP 连接成功、TLS 握手超时（网络路径不可达，非服务端死链：门户页可从外部正常访问）
- 证据：curl -v 输出（Connected ... TLS handshake, Client hello ... Connection timed out）
- 放弃理由：本地网络路径限制，转 wayback 存档恢复（W1.6）

### [W1.6] 2026-11-11 wayback 恢复 ESPRI 组件文件（成功路径）
- 动作：CDX 定位 ESPRI optdat 存档（2022-2025 多次抓存）；按 `web/<ts>id_/<original_url>` 原始字节模式批量取 38 个组件文件；availability API 逐文件钉 capture 时间戳
- 结果：**25/38 取得**（4 个湿度组件：SUSO 6/8、WASO 4/8、SSAM 8/8、SSCM 5/8；另 INSO/MICM/MITR 干态 3 件）；12 个未存档（waso80/90/95/98、suso50/90、sscm70/98/99、soot00、minm00、miam00）
- 证据：25 文件 SHA256 已入 sources.json（每文件带 original_url + retrieved_url + capture timestamp）；结构核验：61 波长（0.25–40 µm）+ 干态参数与 BAMS98 Table 1c 逐值一致 + 原生 OPAC ASCII 格式
- 备注：缺失 12 文件在两个 ESPRI 主机的 wayback 上均 404（跨 2022/2023/2024/2025 时间戳试过）

### [W1.7] 2026-11-11 FlexAOD 镜像（bundled data/gads/optdat）
- 动作：curl http://pumpkin.aquila.infn.it/flexaod/flexaod.tgz
- 结果：返回 Maintenance 页面（服务器维护中）
- 证据：返回 2500 字节 HTML（title: Maintenance）
- 放弃理由：预算闸（该维度候选源第 3 个：ESPRI 直连 / wayback / FlexAOD）；wayback 已提供 25 文件，缺失文件不再追

### [W2.1] 2026-11-11 Dubovik et al. 2002（AERONET 类型气候学）
- 动作：curl https://data.giss.nasa.gov/gacp/publications/special/dubovik3.pdf
- 结果：取得 968 KB PDF（51 页，JAS 排版稿），Table 1 逐站点 τ 范围/α 范围提取成功
- 证据：sha256 6a8875b1969bbf6486f865d03895e5bd77498163999ea7d56f73812617ab97f5（sources.json `dubovik2002_aeronet_types`）

### [W2.2] 2026-11-11 Holben et al. 2001（AERONET 气候学奠基论文）
- 动作：WebSearch 找开放 PDF
- 结果：无开放直链（AGU/Wiley 付费墙；ResearchGate 页存在但不可靠直下）
- 放弃理由：Dubovik 2002 已覆盖同维度（站点级 τ/α 范围），Holben 2001 增量边际；预算闸（第 2 个候选）

### [W2.3] 2026-11-11 AERONET v3 web service 站点数据（GSFC 2011）
- 动作：旧 cgi 端点（data_plot_aod_v3）404 → 经 /new_web/data.html 定位现行端点 print_web_data_v3 → 下载 GSFC 2011 全年 AOD Level 2.0
- 结果：7,205,149 字节 CSV，含 AOD_500/551/560 nm（近 550 nm 直接可比池值域，无需波长换算）+ 440-870/440-675/500-870 α 对 + 站址元数据
- 证据：sha256 f44fba2d714f9274238225828356744806a3f010f6e83dd01103b3ba6e1a236a（sources.json `aeronet_v3_gsfc_2011_lev20`）
- 类比判据落定：四事件季节跨 Jun/Dec/Jan/Feb、观测站纬度带 ~29°S–50°N（ESO La Silla 至俄罗斯中纬度）→ 选取全年 39.0°N 中纬度站点；GSFC 为连续多年 L2.0 数据的代表站

### [W2.4] 2026-11-11 MODIS 全球 AOD 分布产品
- 动作：未发起检索
- 放弃理由：预算闸（AOD 分布维度第 4 个候选；Dubovik 2002 + AERONET 站点文件已构成站点级分布证据，MODIS 全球栅格产品对 pool 值域判定增量边际）。如实记录「预算内未穷尽」

### [W3.1] 2026-11-11 McClatchey et al. 1972（AFCRL-72-0497，经典标准大气气溶胶廓线）
- 动作：WebSearch 定位 → curl https://apps.dtic.mil/sti/pdfs/AD0753075.pdf（两次 UA 尝试）→ wayback 无存档
- 结果：DTIC 返回 "The request is blocked"（主动拦截）；wayback 对该 PDF 无 capture
- 放弃理由：预算闸（垂直几何维度经典报告路径 2 次尝试失败）；转开放获取卫星气候学路径（W3.3）

### [W3.2] 2026-11-11 Elterman 1968（AFCRL-68-0153）
- 动作：WebSearch 找直链
- 结果：无直接 PDF（NTRS/DTIC 仅有引用文献）
- 放弃理由：同 W3.1 预算闸

### [W3.3] 2026-11-11 Winker et al. 2013（CALIPSO 全球 3-D 垂直分布，ACP 开放获取）
- 动作：Yu et al. 2010/Koffi 线索纠偏（JGR 付费墙）→ 改 Winker 2013 → curl acp.copernicus.org 直下
- 结果：取得 5.46 MB PDF（17 页）；Sec. 5 定义 H63（63% AOD 所在高度，与指数层标高构造等价）；Fig. 9 渲染读图：全球 H63 典型 1–3 km（海洋 1–2 km）、极值 ~4.5–5 km、极地高值文中明示 spurious；Fig. 10：海洋背景 AOD532 ~0.05–0.1
- 证据：sha256 4b150733884da4d0d76eca3b63cd98225c481577a891c161431c5302e1f06233（sources.json `winker2013_calipso_3d`）；读图值已在该条目 figure_read_caveat 标注近似性质

### [W1.8] 2026-11-11 OPAC BAMS98 论文独立重下载
- 动作：curl https://pure.mpg.de/rest/items/item_3365752/component/file_3365754/content
- 结果：取得 14.9 MB 扫描版 PDF（14 页，MPG 仓库版）
- 证据：sha256 207aa8f7da4d2d2974879288b58fc7752a1a09485d2ccf58589462468fdd36af（sources.json `opac1998_bams_paper_mpg_copy`）；正文/表格提取经 pymupdf 只读解析
- 备注：canonical 钉定仍在 prior_parameters.json（68a20ba9...，lweb.cfa.harvard.edu 版，121916 字节）；本件为独立第二 artifact，用于本 survey 的表格提取
