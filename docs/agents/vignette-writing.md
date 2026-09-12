# bulkMAE 可执行 vignette 写作经验

本文总结真实分析教程形成过程中已经验证过的写作、渲染和 CI 经验。新增一条分析主线时，先按
[`AGENTS.md`](../../AGENTS.md) 的 Sequence an analysis mainline 节排序。先定主线，再起草
vignette 以列出需要的分析和图。缺分析才加薄封装。最后按本文和
[`plotting-functions.md`](plotting-functions.md) 画图并写完教程。本文只管 vignette 怎么写。

新增分析型 vignette 前，应结合本文阅读
[`vignettes/airway-qc-de.Rmd`](../../vignettes/airway-qc-de.Rmd)、
[`vignettes/enrichment-analysis.Rmd`](../../vignettes/enrichment-analysis.Rmd)、
[`DESCRIPTION`](../../DESCRIPTION) 和相关工作流配置。

## 先区分 API 总览与真实分析教程

两类文档解决的问题不同：

- API 总览帮助用户快速看见包能做什么，可以使用小型模拟数据和有选择地执行代码。
- 真实分析教程必须围绕一个明确研究问题展开，使用真实数据，解释设计选择，并让正文中的结果、表格和图全部由代码产生。

不要把真实教程写成按函数名排列的调用清单。更有效的叙事顺序是：研究设计 → 数据构建 → QC → 预处理 → 模型 → 结果展示 → 限制与复现信息。每个代码块前说明为什么做，代码块后说明能得出什么、不能得出什么。

## 数据应真实、官方且离线可复现

首选官方 ExperimentData 包或稳定的随包数据，不提交容易漂移的派生副本。airway 教程的做法是：

- 数据来自 `airway` ExperimentData 包；
- 安装依赖时可以联网，但安装完成后的渲染不访问网络；
- 原始 counts、样本元数据和 feature 元数据在文中显式提取；
- ID、因子水平和分析所需注释在建模前固定；
- 数据维度和过滤后数量通过 inline R 动态报告，不硬编码预期结果。

正文中的网页链接只用于来源说明，不能成为渲染数据的运行时依赖。BioMart、KEGG、STRING、OmniPath、LINCS 等在线服务不适合作为普通 vignette 构建的必经步骤。

续篇分析应复用上一篇已经定义的数据、过滤、设计和 contrast。airway 富集教程虽然为保证独立执行而
重复最短 DE 代码，但 counts、`~ cell + dex`、`trt - untrt` 和过滤规则必须与 QC/DE 教程一致；不能
为了得到更好看的富集图换成模拟结果或另一套数据。GO 教学使用已安装的 `org.Hs.eg.db`，安装完成后的
渲染不访问网络。

## 先解释实验设计，再调用模型

真实数据不等于真实分析。教程必须把设计结构写清楚，并让代码与文字保持一致。例如 airway 的四个 cell line 各有 untreated/treated 一对样本，因此：

- 先固定 `dex` 的参考水平和比较方向；
- 用 `~ cell + dex` 吸收 cell line 的基线差异；
- 显式请求 `c("dex", "trt", "untrt")`，不依赖因子当前顺序猜测 contrast；
- PCA/相关性使用合适的 VST assay，DE 仍在原始 counts 上拟合；
- 离群标记只作为回查信号，不在教程中自动删除样本。

模型结果也要用同一组对象和阈值驱动正文。显著基因总数、up/down 数和表格应来自 `de_selected()` 与 `de_table()`，不能手抄一次运行结果。

## 防止展示阶段制造过度结论

分析教程很容易在作图时产生“看完结果再定义问题”的偏差。当前教程采用了两条可复用规则：

- 需要标注的 CRISPLD2 来自原研究的预先关注，而不是按当前数据表现临时挑选。
- 热图的 top 30 feature 在教程代码中显式选择并传给绘图函数，同时明确它是同一数据上的描述性展示，不是独立验证。

同样，筛查阈值不是生物学真理，聚类图不是验证队列，显著性也不等于因果。教程结尾应说明样本量、设计范围、注释版本、模型假设和独立验证需求。

富集教程还要避免先看 p-value 再挑图：明确区分 ORA 的阈值化基因列表与完整受检背景、GSEA 的完整
有方向排序，并在查看本次富集显著性前固定需要展示的 term。绘图可按显式 effect 选择每个 term 的
代表 gene，但正文必须说明该选择只改变显示 membership，不改变 ORA 检验或完整 membership 统计量。

## 所有结果都由可执行代码生成

除安装示例外，分析 chunk 应实际执行。不要用 `eval = FALSE` 隐藏昂贵但关键的步骤，也不要把本地运行得到的数字粘进正文。

package vignette 还需要完整的 front matter。中文文档尤其不能漏掉 UTF-8 声明；标题或文件名改变时，也要有意识地同步检查 index entry 与 `_pkgdown.yml`：

```yaml
---
title: "教程标题"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{教程标题}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```

当前 airway 教程的 setup 可作为结构示例；具体 seed 和画布大小应按教程内容选择，而不是复制成所有文章的固定值：

```r
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.8,
  dpi = 120,
  fig.retina = 2
)
set.seed(20260824)
library(bulkMAE)
library(ggplot2)
```

固定 seed、图宽、图高和 DPI 让本地、R CMD build 与 pkgdown 的输出更可比较。需要注意：knitr 会先创建图形设备，所以不会读取 bulkMAE ggplot 上附着的推荐物理尺寸。vignette 内通过 `fig.width`/`fig.height` 控制嵌入图；最终单图文件才通过 `plot_save()` 控制厘米尺寸。

inline R 适合报告：

- 数据维度与样本数；
- 过滤后保留数量和比例；
- 当前 FDR/effect 阈值；
- 入选、up、down 数量；
- 包版本。

表格使用由当前结果排序得到的对象，再交给 `knitr::kable()`。排序时要显式处理 `NA`，并用 `min()` 防止数据不足时越界。

## 图形既要可访问，也要统计上准确

每幅图使用独立 chunk，并同时给出：

- 与正文语言一致的 `fig.cap`，说明图展示什么（当前 airway 教程使用中文）；
- 与正文语言一致且有意义的 `fig.alt`，说明读者看不到图时仍需要知道的编码和内容；
- 必要时单独设置合适的 `fig.height`；
- 正文中的解释和限制。

包绘图函数提供英文轴名与图例，教程通过 `+ labs(title = ...)` 添加上下文标题。不要调用 `DESeq2::plotMA()`、基础 `plot(prcomp)`、pheatmap 等另一个绘图后端，否则教程无法同时验证统一绘图层。

图注必须描述真实编码。例如 Volcano 的颜色由 adjusted p-value 与 effect 阈值决定，而 y 轴是原始 p-value 的 `-log10`；这两者不能写成同一个阈值。修改函数语义后，应同步检查图注、alt text 和紧邻正文。

## 依赖声明与构建工作流必须成套更新

能在开发机渲染不代表包构建会成功。新增 vignette 依赖时至少检查四处：

1. 数据包和可选分析后端列入 `DESCRIPTION` 的 `Suggests`；渲染引擎列入 `Suggests`，并保持 `VignetteBuilder: knitr`。
2. R CMD check 工作流显式安装教程真正需要的 `knitr`、`rmarkdown`、数据包和后端。
3. pkgdown 工作流安装同一组教学依赖，并真实构建文章。
4. `_pkgdown.yml` 把文章加入预期顺序，README 提供稳定入口。

本仓库有大量不在标准 CRAN/Bioconductor 仓库中的可选 Suggests。pkgdown 工作流不能为了省事使用 `dependencies: "all"`，否则与教程无关的 GitHub-only 后端也可能阻断网站构建。当前做法是 `dependencies: "hard"`，再在 `extra-packages` 中明确列出教程依赖。完整后端和 coverage 工作流可以使用 `"all"`，但必须同时提供那些非标准后端的固定 GitHub 来源。

普通教程构建继续设置 `BULKMAE_RUN_ONLINE_TESTS=false`。该变量只影响主动读取它的测试，不是网络沙箱，也不能保证 vignette 离线；离线性必须由教程自己的数据来源和函数调用保证。

## 分层验证比只点 Knit 更可靠

建议按以下顺序验收：

1. 在干净 R 会话中完整执行 Rmd，确认只有安装示例不运行。
2. 检查每个结果数字是否来自 inline R，每幅图是否有独立 `fig.cap` 和 `fig.alt`。
3. 从源码构建包，使 vignette 在真正的 build 环境中执行：

   ```sh
   R CMD build --no-manual --compact-vignettes=gs+qpdf .
   ```

4. 对生成的 tarball 运行 check：

   ```sh
   R CMD check --no-manual --as-cran bulkMAE_VERSION.tar.gz
   ```

   其中 `VERSION` 替换为上一步实际生成的版本号。本地必须确认最终为 `Status: OK`；CI 通过 `error-on: "warning"` 强制 warning 也使任务失败。

5. 用 pkgdown 的实际配置构建站点，确认文章顺序、图片嵌入和链接。
6. 打开最终 HTML，人工确认结果表、动态摘要、正文排版、alt text 和预期图数；airway QC/DE 教程
   应有 7 幅图，富集教程应有 5 幅图。

直接 Knit 主要验证作者当前环境；`R CMD build` 会在打包路径中实际执行 vignette；只有在干净依赖环境中结合 `R CMD check`，才能可靠发现漏声明依赖。pkgdown 还会验证网站配置。几条路径覆盖的失败面不同。

构建得到的 tarball、check 目录和 pkgdown `docs/` 是派生产物，不应作为教程源码提交。当前 `.gitignore` 已忽略这些路径；预览完成后应清理构建产物，但不要误删仍在使用的 R、Quarto 或编辑器临时目录。

## 常见失败及原因

| 现象 | 根因 | 处理方式 |
|---|---|---|
| 本地能 Knit，CI 缺包 | 依赖只安装在开发机，未进入 Suggests/工作流 | 同步更新 DESCRIPTION 与构建工作流 |
| pkgdown 因无关后端失败 | 对大量非标准 Suggests 使用 `dependencies: "all"` | 使用 hard dependencies，并显式列教程依赖 |
| 正文数字过期 | 把一次运行结果写死 | 用 inline R 从当前对象计算 |
| 图能显示但网站布局失衡 | 误以为 ggplot 推荐尺寸会控制 knitr 设备 | 在 chunk 中设置 `fig.width`/`fig.height` |
| 图注与统计量矛盾 | 把 raw p-value、adjusted p-value 或 effect 阈值混为一谈 | 从绘图实现反查轴与颜色的真实语义 |
| 教程看似完整但没有执行关键分析 | 大量 chunk 使用 `eval = FALSE` 或缓存旧结果 | 只让安装示例跳过，并从干净会话构建 |
| 教程暗示独立验证 | 在同一数据中筛选并展示 top features | 明确称为描述性展示并陈述验证需求 |
| 渲染依赖网络状态 | 在 chunk 中请求远程数据库 | 改用已安装数据包或稳定的随包资源 |

## 新增真实分析 vignette 检查表

- [ ] 有明确研究问题和实验设计说明，而不是函数清单。
- [ ] 使用官方、可追溯且安装后离线可用的数据。
- [ ] 因子水平、配对/批次结构和 contrast 显式固定。
- [ ] 若文章接续上一篇，数据、过滤、模型和 contrast 与上游教程一致。
- [ ] 除安装示例外，所有分析 chunk 实际执行。
- [ ] seed、图形设备参数和关键阈值明确。
- [ ] 数字、比例、结果表和包版本由代码动态生成。
- [ ] 每幅图有与正文语言一致的独立 caption 和有意义的 alt text。
- [ ] 标签基因有事前依据，数据驱动热图明确标为描述性。
- [ ] ORA 背景和 GSEA 完整排序明确；展示 term 不在查看本次 p-value 后临时挑选。
- [ ] 结尾包含研究限制、数据/论文来源和 `sessionInfo()`。
- [ ] DESCRIPTION、R CMD check、pkgdown 与文章索引同步更新。
- [ ] 从源码完成 build、check 和 pkgdown 构建，无 warning。
- [ ] 不提交 HTML 图片目录、tarball、check 目录或临时预览图。
