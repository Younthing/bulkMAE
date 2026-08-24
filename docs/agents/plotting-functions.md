# bulkMAE 绘图函数重构经验

本文记录统一 ggplot2 绘图层形成过程中已经验证过的接口边界、实现方式和验收方法。新增或修改
绘图相关接口前，应先阅读本文，并以当前的
[`R/plotting.R`](../../R/plotting.R) 和
[`tests/testthat/test-plotting.R`](../../tests/testthat/test-plotting.R) 为可执行事实来源。

## 先确定绘图层的职责

bulkMAE 的分析函数是对原生后端的薄封装，绘图层也应保持同样的边界：它负责把已经计算好的结果转换为一致、可继续组合的图，不负责偷偷重新分析数据。

稳定契约如下：

- 每个图形构造函数返回一个未打印、未保存的标准 `ggplot` 对象；这里不包括负责落盘的 `plot_save()`。
- 不创建自定义结果类，不实现 `autoplot()`，也不改变分析函数的返回类型。
- 图形构造函数不调用 `theme_set()`、不开图形设备，也不写文件。
- 标题和上下文解释由调用者通过 `+ labs()` 添加；包内提供英文轴名、图例名和通用 alt text。
- 用户可以用普通的 ggplot2 `+` 语法替换主题、比例尺、标签或局部样式。

这条边界让绘图层保持无状态，也避免为了画图而迫使所有 DE 后端转换成新的 bulkMAE 结果类。例如，MA 图需要后端特有的丰度列，但不应因此扩充 `de_table()` 的六列公共契约。

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

少量 vdiffr 快照可以作为布局意外漂移的报警器，但不能代替视觉验收。当前代表性快照见
[`tests/testthat/test-plotting-vdiffr.R`](../../tests/testthat/test-plotting-vdiffr.R)。

需要人工确认时，将 PNG 生成到被忽略的 `docs/plot-previews/`，使用最终推荐尺寸和固定 DPI；不要写入临时目录，也不要提交这些预览图。由维护者确认文字是否拥挤、图例是否侵占面板、横轴标签是否重叠以及正负尺度是否均衡。自动化测试负责可计算契约，最终视觉判断由人完成。

agent 不得自行接受视觉变化或更新 vdiffr baseline；应先生成 PNG，并取得维护者确认。

## 新增图形构造函数检查表

- [ ] 返回标准、未打印的单个 `ggplot`。
- [ ] 不重新执行分析，不改变上游结果类型。
- [ ] 标题留给调用者；轴、图例和 alt text 完整。
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
