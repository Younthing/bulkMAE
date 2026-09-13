# Save a bulkMAE plot at its recommended physical size

Plot themes control text and visual styling, but cannot control the
physical size of a graphics device. Each bulkMAE `plot_*()` helper
therefore attaches a recommended width and height in centimetres to its
otherwise standard ggplot object. `plot_save()` reads that
recommendation and delegates to
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
with `scale = 1` by default. Explicit `width` or `height` values
override either recommendation.

## Usage

``` r
plot_save(
  filename,
  plot,
  width = NULL,
  height = NULL,
  units = c("cm", "in", "mm"),
  dpi = 300,
  scale = 1,
  device = NULL,
  bg = NULL,
  ...
)
```

## Arguments

- filename:

  Output file name, including an extension understood by
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

- plot:

  A ggplot object. Plots not produced by bulkMAE must supply both
  `width` and `height`.

- width, height:

  Optional positive physical dimensions. When omitted, use the
  recommendation attached by a bulkMAE plotting helper.

- units:

  Physical units for `width` and `height`.

- dpi:

  Raster resolution passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

- scale:

  Multiplicative scale passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).
  Keep the default of one when exact final dimensions matter.

- device:

  Optional graphics device. When `NULL`, infer it from `filename`, with
  Cairo PDF preferred when available.

- bg:

  Background colour passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

- ...:

  Additional arguments passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Value

`filename`, invisibly.

## Details

PDF output uses
[`grDevices::cairo_pdf()`](https://rdrr.io/r/grDevices/cairo.html)
automatically when Cairo graphics are available. Dimensions remain
attached after ordinary ggplot2 `+` operations, so themes, labels, and
scales can be changed before saving.

The RStudio plot pane and knitr graphics chunks create their graphics
devices before a plot is drawn, so they do not consume the attached
recommendation. Use `plot_save()` for an exact final output size; set
`fig.width` and `fig.height` separately when a rendered document must
use a particular preview size.

## See also

Other plotting:
[`plot_assay_heatmap()`](https://younthing.github.io/bulkMAE/reference/plot_assay_heatmap.md),
[`plot_de_ma()`](https://younthing.github.io/bulkMAE/reference/plot_de_ma.md),
[`plot_de_volcano()`](https://younthing.github.io/bulkMAE/reference/plot_de_volcano.md),
[`plot_embedding()`](https://younthing.github.io/bulkMAE/reference/plot_embedding.md),
[`plot_gsea_classic()`](https://younthing.github.io/bulkMAE/reference/plot_gsea_classic.md),
[`plot_gsea_ridge()`](https://younthing.github.io/bulkMAE/reference/plot_gsea_ridge.md),
[`plot_lincs_heatmap()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_heatmap.md),
[`plot_lincs_overlap()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_overlap.md),
[`plot_lincs_rank()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_rank.md),
[`plot_lincs_waterfall()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_waterfall.md),
[`plot_ora_bubble()`](https://younthing.github.io/bulkMAE/reference/plot_ora_bubble.md),
[`plot_ora_network()`](https://younthing.github.io/bulkMAE/reference/plot_ora_network.md),
[`plot_ora_radial()`](https://younthing.github.io/bulkMAE/reference/plot_ora_radial.md),
[`plot_qc_correlation()`](https://younthing.github.io/bulkMAE/reference/plot_qc_correlation.md),
[`plot_qc_library()`](https://younthing.github.io/bulkMAE/reference/plot_qc_library.md),
[`plot_qc_outliers()`](https://younthing.github.io/bulkMAE/reference/plot_qc_outliers.md),
[`theme_bulkmae()`](https://younthing.github.io/bulkMAE/reference/theme_bulkmae.md)
