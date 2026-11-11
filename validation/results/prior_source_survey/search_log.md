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

### [W1.8] 2026-11-11 OPAC BAMS98 论文独立重下载
- 动作：curl https://pure.mpg.de/rest/items/item_3365752/component/file_3365754/content
- 结果：取得 14.9 MB 扫描版 PDF（14 页，MPG 仓库版）
- 证据：sha256 207aa8f7da4d2d2974879288b58fc7752a1a09485d2ccf58589462468fdd36af（sources.json `opac1998_bams_paper_mpg_copy`）；正文/表格提取经 pymupdf 只读解析
- 备注：canonical 钉定仍在 prior_parameters.json（68a20ba9...，lweb.cfa.harvard.edu 版，121916 字节）；本件为独立第二 artifact，用于本 survey 的表格提取
