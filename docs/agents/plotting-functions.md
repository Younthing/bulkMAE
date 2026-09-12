# bulkMAE 绘图函数重构经验

本文记录统一 ggplot2 绘图层形成过程中已经验证过的接口边界、实现方式和验收方法。新增或修改
绘图相关接口前，应先阅读本文，并以当前的
[`R/plotting.R`](../../R/plotting.R)、
[`R/plotting-gsea.R`](../../R/plotting-gsea.R)、
[`R/plotting-ora.R`](../../R/plotting-ora.R) 及其对应的
[`tests/testthat/test-plotting.R`](../../tests/testthat/test-plotting.R)、
[`tests/testthat/test-plotting-gsea.R`](../../tests/testthat/test-plotting-gsea.R)、
[`tests/testthat/test-plotting-ora.R`](../../tests/testthat/test-plotting-ora.R) 为可执行事实来源。

## 先确定绘图层的职责

bulkMAE 的分析函数是对原生后端的薄封装，绘图层也应保持同样的边界：它负责把已经计算好的结果转换为一致、可继续组合的图，不负责偷偷重新分析数据。

稳定契约如下：

- 每个图形构造函数返回一个未打印、未保存的标准 `ggplot` 对象；这里不包括负责落盘的 `plot_save()`。
- 不创建自定义结果类，不实现 `autoplot()`，也不改变分析函数的返回类型。
- 图形构造函数不调用 `theme_set()`、不开图形设备，也不写文件。
- 标题和上下文解释由调用者通过 `+ labs()` 添加；包内提供英文轴名、图例名和通用 alt text。
- 用户可以用普通的 ggplot2 `+` 语法替换主题、比例尺、标签或局部样式。

这条边界让绘图层保持无状态，也避免为了画图而迫使所有 DE 后端转换成新的 bulkMAE 结果类。例如，MA 图需要后端特有的丰度列，但不应因此扩充 `de_table()` 的六列公共契约。

### 参考图脚本是视觉规格，不只是图型名称

当维护者提供现有绘图脚本时，应先把它转换成逐项验收表，再设计 package Interface。至少记录：

- 每一个统计层和装饰层及其绘制顺序；
- 固定颜色、渐变端点、透明度、线宽、点大小与形状；
- 坐标范围、面板比例、标签换行、图例位置和 plot margin；
- 最终物理宽高、DPI，以及这个尺寸下的实际字号；
- 数据选择、排序和节点归属规则。

不能只保留“山脊图”“网络图”等抽象类别，然后用包内通用主题和色板重新设计。参考脚本已经包含
成图经验；除非维护者要求重新设计，默认输出应尽量复现这些具体决定。允许的偏离必须逐项说明，
通常只包括统计纠错、严格 ID 连接、统一 6 pt，以及为返回单个标准 `ggplot` 所需的结构替换。
结构测试应直接断言关键 layer 参数；最终仍生成代表性 PNG 交由维护者判断视觉效果。

### 五个富集图的逐图复核结论

本轮按“一个函数完成源码比对、修改、测试和单图 PNG 后，再进入下一个函数”的顺序复核。后续
维护时应保留以下边界：

| 函数 | 应与参考源码一致 | 已纠正的漂移 | 可以且应保留的差异 |
|---|---|---|---|
| `plot_gsea_classic()` | 0.50/0.20/0.30 三段比例、ES 渐变、hit barcode、200 格 rank 色带、灰色 rank polygon、细边框、无图例和 `5.9 × 5.3 cm` | 逐层复核未发现额外漂移 | 单个 free-space facet 代替 patchwork；自行重建并校验 ES，不调用 `enrichplot:::`；term 显式指定且函数不保存文件 |
| `plot_gsea_ridge()` | `geom_ridgeline()`、0.36 高度、0.18 填充透明度、命中刻线、七色顺序、无图例和 `12 × 5 cm` | GSEA 表同时识别 `qvalue`、`qvalues`、`q_value` | term 和 membership 显式指定；不自动选“正向 top term”；常量 rank 分布报错，不添加假 jitter |
| `plot_ora_bubble()` | 40 字符换行、0.8–2.8 点面积范围、蓝到暗红证据渐变、0.1% 百分比标签、右侧图例、细边框和 `10 × 8 cm` | 固定 Gene count 在上、证据色条在下；证据标题缩写为 `-log10(adj p)`；主刻度约束为 4 个以避免最终尺寸下重叠；显式白底 | 保留正确的 rich factor、gene ratio、fold enrichment 定义和轴名；term 由调用者显式选择并排序；不复制参考脚本中错误的 denominator 或统计量 fallback |
| `plot_ora_network()` | term-specific feature 视觉节点、社区形状、Jaccard term 边、15 层 halo、membership 线、节点范围、默认不标 gene、无图例、白底和 `8 × 8 cm` | 恢复共享 feature 在每个所属社区各画一次；term-named `features` 可逐条目选择显示成员；补齐白底与文字背景参数 | 完整 membership 仍只保留一份并计算 Jaccard；feature 选择显式；确定性布局不读写全局 RNG；gene 标签仅作为显式扩展 |
| `plot_ora_radial()` | 唯一外圈 feature 节点、shared-count term 边、6.3/4.1/6.65 半径、0.075 扇区间隔、节点与 membership 线范围、20 字符换行、径向标签旋转、白底和 `8 × 8 cm` | 恢复 term edge 的 `shared_n^0.85` 强度；使用径向源码独有的短色板规则；补齐白底、字符串换行和 term 文字参数 | 所有几何文字统一为维护者要求的 6 pt，而不复制源码的 5.5/4 pt；gene 标签与显示集合仍由 ID 显式指定 |

网络图与径向图的短色板规则看似相近但参考脚本确实不同，因此实现中不能再次合并为同一策略。
外圈背景弧、中心节点和中心连线在径向参考脚本中默认关闭；bulkMAE 复现默认成图，不为这些关闭的
装饰层增加浅公共参数。

## 主题、字号与物理尺寸是三件事

### 主题只控制样式

`theme_bulkmae()` 是包内定制主题，默认 `base_size = 6` pt。它定义基础字体、网格、图例位置等视觉规则，但 ggplot2 主题无法决定图形设备的宽高。

主题文字本来就使用 pt；数据层文字必须显式写单位：

```r
ggplot2::geom_text(
  size = 6,
  size.unit = "pt"
)
```

否则 `theme()` 与 `geom_text()` 的数字看似相同，实际单位却可能不同。发表图的实用原则是：**先锁定最终物理尺寸，再锁定字号；主题文字由 `theme()` 管，数据标签由 `geom_text(size.unit = "pt")` 管。**

### 推荐尺寸附着在普通 ggplot 上

每个图形构造函数用 `bulkmae_dimensions` 属性附加以厘米为单位的推荐宽高。对象仍然继承 `ggplot`，经过普通的 `+` 操作后属性仍应存在。固定构图使用固定尺寸；样本数或特征数会影响可读性的图，则使用有上下限的动态尺寸，避免无限扩张。

属性是输出建议，不是新的绘图类，也不会被 RStudio Plot 面板或 knitr 自动读取，因为这些工具会在绘图前创建图形设备。

### 最终文件由 `plot_save()` 锁定尺寸

`plot_save()` 读取推荐尺寸并调用 `ggplot2::ggsave()`；默认 `scale = 1`，可用时 PDF 自动选择 Cairo。期刊规定尺寸时，调用者可显式覆盖宽高：

```r
p <- plot_de_volcano(result) +
  ggplot2::labs(title = "Treatment effect")

plot_save("volcano.pdf", p)

journal_width_cm <- 9
journal_height_cm <- 7
plot_save(
  "volcano.pdf", p,
  width = journal_width_cm,
  height = journal_height_cm,
  units = "cm"
)
```

对于不是 bulkMAE 生成的普通 ggplot，必须同时提供 `width` 和 `height`。文档预览应另行设置 chunk 的 `fig.width` 和 `fig.height`；不要误以为主题或对象属性能够改变已经打开的设备。

## 标识符连接必须显式

生物信息图最危险的错误通常不是语法错误，而是元数据与数值错位。绘图函数不得按当前位置猜测对应关系：

- 样本名和 feature ID 必须非空、非缺失且唯一。
- 元数据按 row name、feature name 或显式 ID 连接，并验证接口所要求的覆盖范围。
- 允许输入顺序不同，但连接后必须按主数据的顺序重排。
- embedding 的外部 `sample_data` 和外部命名向量采用精确集合契约，缺失或多出 ID 都报错。
- MAE 容器内的 metadata 列可以包含未入图的行，但必须覆盖全部已选 ID，再按已选 ID 取值。
- 重复 ID、缺失 ID、空 ID 或需要命名却未命名的向量应尽早报错。
- 只连接用户实际映射的元数据列，避免无关列覆盖绘图内部的 `sample`、`x`、`y` 等保留名称。
- `shape` 只接受离散值；连续变量应要求用户先转成因子，而不是让 ggplot2 在构建阶段才报难懂的错误。

这个规则同样适用于热图的 `column_split`、feature label 和 MA 图的命名 abundance 向量。

## 各图族需要保留的统计语义

### Embedding 与 QC

- `plot_embedding()` 只做输入适配和展示。`prcomp` 的轴名动态包含解释率；limma MDS、Rtsne 和命名矩阵保留各自语义。
- 相关性矩阵必须是行列样本集合相同的有限方阵。聚类只重排显示顺序，不改变数值；较长横轴标签默认旋转 45 度。
- 离群图中的阈值与 `flagged` 是审查信号，不应在绘图函数中触发样本删除。
- 离散色使用包内固定的色盲友好颜色，连续色使用 viridis，相关性使用固定发散色阶。

### Volcano 与 MA

两类 DE 图必须共用一套逐行显著性分类语义；当 adjusted p-value 可用时，它与 `de_selected()` 保持一致：

```text
adjusted p-value is not missing AND adjusted p-value <= fdr
AND effect is finite AND effect != 0
AND abs(effect) >= min_abs_effect
```

绘图时，缺失 adjusted p-value 视为不显著，因此即使 adjusted p-value 全部缺失，仍可按 Not significant 展示其它可绘制信息；`de_selected()` 面对完全不可用的 adjusted p-value 则会报错。effect 恰好为零不能归到 up 或 down。颜色和因子顺序固定为蓝色 Down、灰色 Not significant、橙红色 Up，图例放在上方并横排。

还要保留下列边界：

- Volcano 的 y 轴是原始 p-value 的 `-log10`，FDR 只用于点的分类。因此不能在 `-log10(fdr)` 处画所谓“FDR 阈值线”；BH FDR 与单个原始 p-value 没有这种等价关系。
- 原始 p-value 为零时，应使用有限的正数下限：优先取最小正 p-value 的十分之一，没有正值时回退到机器下限；不能得到无限坐标，全为 `NA` 时应报错。
- Volcano 默认以零为中心使用对称 x 轴，便于比较正负效应；保留显式关闭选项。
- 只有用户明确传入 `label_features` 时才加标签，不在看完结果后自动挑前 N 个。
- DESeq2 LRT 没有单一方向，应拒绝 up/down 着色并提示改用 Wald 结果。
- 通用 data frame 的 `F`/`LR` 列因无法安全推断自由度，一律保守地视为无方向；edgeR 单 coefficient 结果可利用 `logFC` 方向绘制，multi-coefficient edgeR 则拒绝。
- MA 图的 y 轴统一叫 `Effect estimate`，不能假设所有模型的效应都是 log2 fold change。
- DESeq2 使用 `baseMean` 并显示 log10 丰度轴；edgeR 使用原尺度 `logCPM`；limma/voom/dream 使用原尺度 `Amean`。
- 普通表只有在丰度列唯一且无歧义时才自动识别，否则要求完整命名的 `abundance` 或明确的 `abundance_column`。
- log10 丰度轴可以警告并忽略零值，但必须拒绝负值。因子型丰度必须报错，不能把 level 编码误当数值。

新增后端适配时，应增加后端测试，同时断言 `de_table()` 仍精确返回：

```text
feature_id, effect, standard_error, statistic, p_value, adjusted_p_value
```

### Assay heatmap

热图只负责展示调用者明确给出的 feature 集合：

- `features` 必填且不能为空，不在函数内部选择“最显著基因”。
- 行标准化遇到常量行时先警告再删除；若无可绘制行则报错。
- 聚类只决定顺序，不绘制树状图。
- 分组和 feature label 仍按名称连接；重复的显示标签用 `make.unique()` 区分，避免多个 feature 合并到同一离散轴位置。
- 行标准化值使用固定发散色阶，未标准化 assay 值使用连续 viridis 色阶。

### GSEA 细节图

GSEA 细节图需要完整 ranked list、命中位置和 running enrichment score，因此
`plot_gsea_classic()` 与 `plot_gsea_ridge()` 使用专门的 ranked-list Adapter。两者共享 ranks、
gene-set membership、leading edge 和证据解析，但保留两个语义明确的公共入口；不要用一个
`view = "classic"/"ridge"` 的大函数把条件参数暴露给调用者。

- 原生 clusterProfiler `gseaResult` 已保存 `geneList`、`geneSets` 和 exponent，可以准确重建图形。
- fgsea 结果表不保存完整 ranks、pathways 或分析时的 `gseaParam`。classic 必须要求调用者显式
  补齐 `ranks`、`gene_sets` 和 `exponent`；不能静默假定 exponent 为 1。
- ranks 必须有限、具名、ID 唯一且已经按非递增顺序排列。绘图层验证但不暗中重排。
- running score 在包内按 weighted hit/miss walk 计算，并与结果中的 ES 做容差校验；不要调用
  `enrichplot:::` 私有函数。
- classic 的 ES 曲线、hit barcode/排序热带和 ranked metric 是同一个统计图型的三个不可拆层，
  不是报告级总拼图。实现用一个标准 ggplot 的 free-space facet，并固定约 0.50/0.20/0.30 层高。
- classic 保留参考图的绿到红 ES 曲线、蓝白红 rank 热带、灰色 rank polygon、细黑边框和
  `5.9 × 5.3 cm` 推荐尺寸。用 facet 代替 patchwork 是为维持标准 `ggplot` 返回契约，不是重新设计图。
- classic 内部为保持三套真实 y 刻度而使用归一化显示带；主题、标题和颜色仍可用 `+` 替换，
  但调用者不应替换内部 y scale。
- ridge 的 x 是 gene-level rank statistic，不是表达量。每条 ridge 在 term 内归一化，只比较分布
  位置与形状；不能用 ridge 高度比较 set size 或富集强度。
- ridge 保留参考脚本的七色顺序、`0.36` 标准化高度、`0.18` 填充透明度、命中刻线和
  `12 × 5 cm` 推荐尺寸；默认关闭原脚本也默认关闭的统计文字框。
- leading edge 与完整 gene set 是两种不同 membership，必须由参数明确选择。少于两个可用且
  不同的 rank 值时直接报错，不通过人为 jitter 伪造核密度。

### ORA bubble 与 membership graph

ORA 专用图需要 `GeneRatio`、`BgRatio`、`Count` 和 enriched-feature membership，因此使用独立
的 ORA Adapter，不向通用 enrichment view 追加大量只对某个后端有意义的可空列。

三个常见横轴必须严格区分：

```text
Gene ratio      = Count / input gene count
Rich factor     = Count / background term size
Fold enrichment = Gene ratio / background ratio
```

clusterProfiler 的 ratio 字符串为 `numerator/denominator`；background term size 是 `BgRatio`
的分子，不是分母。解析时应验证正分母、ratio 范围以及 `GeneRatio` 分子与 `Count` 一致，不能在
某列缺失时把另一种统计量作为 fallback 却沿用原轴名。

`plot_ora_network()` 与 `plot_ora_radial()` 共享同一份经过验证的完整 membership，但两种布局保留
参考源码各自的视觉节点与 term-edge 语义：

- `geneID` 表示富集输入中命中的 feature，不是完整 pathway membership；帮助页和 alt text 应称为
  enriched-feature overlap。
- 内部 feature ID 保持唯一并与 term ID 分开命名；community 布局再把每条显示 membership 展开成
  term-specific 视觉节点，因此共享 feature 会在每个所属社区各出现一次。radial 布局则保留一个
  外圈视觉节点，并连接全部相邻 term。
- network 的 term-term edge 表示完整 enriched-feature membership 的 Jaccard coefficient；radial
  内圈 edge 的宽度和透明度按共享 feature 数的 `shared_n^0.85` 缩放。二者都不是 ontology semantic
  similarity，也不能混写成同一统计量。
- `features` 可以是全局 character vector，也可以是逐 term 命名的 list；它只限制显示 membership，
  不改变用于 network Jaccard 的完整 membership。
- 外部 `feature_values` 按名称覆盖实际显示 feature，允许携带未显示的全基因组额外值；为复现参考
  网络风格，节点大小和透明度使用绝对值，节点颜色表示其构图 community，不伪装成带方向的 effect 图例。
- community layout 的 halo、radial layout 的外圈 ownership 都只是降低遮挡的构图信息，没有
  inferential 含义；radial 中共享 feature 仍保留通向所有 term 的 membership edge。
- term 和 feature 标签仍由显式 ID 连接。`label_features = NULL` 不自动标注“top gene”；参考 network
  默认不标 gene，而参考 radial 可由调用者显式标注全部显示 gene。
- bubble 保留蓝到暗红证据渐变、细黑边框、浅灰读数网格、`adj p` 紧凑色条、四个主刻度和
  `10 × 8 cm` 尺寸；network/radial 保留低饱和 community 色、无图例的 `theme_void()`、显式白底、
  halo/内外圈几何和 `8 × 8 cm` 尺寸。

总拼图脚本只用于吸收单图比例、边框、留白和推荐物理尺寸经验。bulkMAE 不提供 GSEA/ORA
overview 或报告拼图函数；调用者可在包外使用 patchwork/cowplot 组合返回的独立 ggplot。

## 实现时优先复用内部语义

公共接口应保持小而清楚，把适配和校验复杂度封装在私有 helper 中，例如：

- 输入列和 ID 校验；
- 推荐尺寸附加与单位换算；
- embedding 输入适配；
- DE 表提取、显著性分类和颜色比例尺；
- feature 标签解析；
- abundance 后端适配；
- MAE 样本/feature 注释对齐。

复用的目标是让同一个统计概念只有一套语义。优先共享 helper；若现阶段仍有独立实现，就必须用 parity test 防止它们漂移。不要为了减少几行代码而把互不相干的图族塞进一个带大量分支的公共函数。

## 验收分为数值、结构和视觉三层

### 结构与行为测试

绘图测试组合应按适用图族覆盖：

- 每个图形构造函数都返回标准 `ggplot`，且 `ggplot2::ggplot_build()` 成功；
- 轴名、图例名、alt text、因子顺序、固定颜色和阈值分类正确；
- 推荐尺寸存在、单位为 cm，且经过 `+` 修改后仍保留；
- 为 `plot_save()` 做一个代表性集成测试：以固定 DPI 写出 PNG，并从文件头核对像素宽高确实等于推荐厘米数换算后的结果；
- 乱序元数据按 ID 正确连接；重复、缺失，以及精确集合输入中的多余 ID 会报错；
- 零 p-value、全 `NA`、常量热图行、空 feature 集以及不支持结果类型具有明确行为；
- DESeq2、edgeR、limma 等真实后端的 abundance 适配正确，且未改变公共 DE 表。

不要通过截图证明数值语义；直接检查 plot data、scale、layer 或构建后的对象更稳定。

### 视觉回归与人工确认

本仓库不把 vdiffr 或 SVG 快照作为视觉验收基线。需要人工确认时，将 PNG 生成到被忽略的
`docs/plot-previews/`，使用最终推荐尺寸和 600 dpi；不要写入临时目录，也不要提交这些预览图。
DE 文献图的预览应使用 airway 真实结果，不用 `mae_simulate()`。入门教程的合成图除外。
由维护者确认文字是否拥挤、图例是否侵占面板、横轴标签是否重叠以及正负尺度是否均衡。自动化
测试只负责可计算的数值、比例尺和图层结构契约，最终视觉判断由人完成。

agent 不得自行接受视觉变化；应先生成对应的单幅 PNG，并交由维护者确认。

## 新增图形构造函数检查表

- [ ] 返回标准、未打印的单个 `ggplot`。
- [ ] 不重新执行分析，不改变上游结果类型。
- [ ] 报告级标题留给调用者；若 term 名是图形本体的一部分，可作为内置标题；轴、图例和 alt text 完整。
- [ ] 所有元数据按名称连接，并验证该接口规定的精确集合或覆盖关系。
- [ ] 颜色、因子顺序和统计分类可预测。
- [ ] 数据标签显式使用 `size.unit = "pt"`。
- [ ] 推荐物理尺寸合理，必要时随 item 数量有界增长。
- [ ] `plot_save()` 能按 `scale = 1` 导出最终尺寸。
- [ ] 正常输入能 `ggplot_build()`，关键错误输入有测试。
- [ ] 生成忽略的 PNG 供人工视觉确认。

相关事实来源还包括
[`R/result-inputs.R`](../../R/result-inputs.R)、
[`README.md`](../../README.md) 的 publication-size 示例，以及上文链接的结构测试和视觉快照测试。
