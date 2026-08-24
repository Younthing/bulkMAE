# bulkMAE 依赖安装指南（pak）

本指南使用 [`pak`](https://pak.r-lib.org/) 作为唯一的 R 包依赖安装器。
`bulkMAE` 的 MAE/SE 基础设施是必需依赖，分析后端是按函数调用加载的可选
依赖。只安装项目实际使用的后端即可；发布验证环境应安装全部后端。

以下官方仓库地址核验于 2026-08-24。GitHub 默认分支会移动，正式分析应把
实际使用的 commit 记录到项目锁文件或容器定义中。

## 1. 准备 pak 和兼容的 R/Bioconductor 环境

若尚未安装 `pak`，先执行一次：

```r
install.packages("pak")
```

之后的包安装全部使用 `pak::pkg_install()` 或 `pak::local_install()`。pak 会为
当前 R 版本选择相应的 Bioconductor 仓库；不要把不同 Bioconductor 发布线的
包混装到同一库中。安装前先确认当前 R/Bioconductor 组合能够满足
`DESCRIPTION` 中的最低版本，尤其是 edgeR、GSVA 和 decoupleR。

从仓库根目录安装仅含必需依赖的 `bulkMAE`：

```r
pak::local_install(".", dependencies = NA)
```

## 2. 标准 CRAN/Bioconductor 后端

以下包均可由 pak 从当前 R 对应的标准仓库解析。可以按项目选择子集，也可以
一次安装完整的标准后端集合。

```r
cran_backends <- c(
  "cran::glmnet",
  "cran::lme4",
  "cran::metafor",
  "cran::msigdbr",
  "cran::NMF",
  "cran::Rtsne",
  "cran::survival",
  "cran::timeROC",
  "cran::uwot",
  "cran::WGCNA"
)

bioc_backends <- c(
  "bioc::AnnotationDbi",
  "bioc::Biobase",
  "bioc::biomaRt",
  "bioc::clusterProfiler",
  "bioc::ConsensusClusterPlus",
  "bioc::decoupleR",
  "bioc::DESeq2",
  "bioc::diffcoexp",
  "bioc::edgeR",
  "bioc::fgsea",
  "bioc::GENIE3",
  "bioc::goseq",
  "bioc::GSVA",
  "bioc::limma",
  "bioc::maSigPro",
  "bioc::OmnipathR",
  "bioc::ReactomePA",
  "bioc::RUVSeq",
  "bioc::signatureSearch",
  "bioc::singscore",
  "bioc::STRINGdb",
  "bioc::sva",
  "bioc::tximport",
  "bioc::variancePartition"
)

pak::pkg_install(c(cran_backends, bioc_backends))
```

`annotation_orgdb()`、`annotate_ids()` 和 GO 富集还需要与物种匹配的本地
`OrgDb` 数据包。例如人和小鼠项目可分别安装：

```r
pak::pkg_install(c("bioc::org.Hs.eg.db", "bioc::org.Mm.eg.db"))
```

不要为了让示例运行而默认安装错误物种的注释库。

## 3. 官方 GitHub 后端

MuSiC、immunedeconv 和 BayesPrism 不保证出现在每一条标准
CRAN/Bioconductor 发布线中。使用 pak 的规范 GitHub 包引用安装：

```r
github_backends <- c(
  "MuSiC=github::xuranw/MuSiC",
  "immunedeconv=github::omnideconv/immunedeconv",
  "BayesPrism=github::Danko-Lab/BayesPrism/BayesPrism"
)

pak::pkg_install(github_backends)
```

这些规格分别对应维护者提供的安装位置：

- [MuSiC 官方仓库](https://github.com/xuranw/MuSiC)：R 包位于仓库根目录；
- [immunedeconv 官方仓库](https://github.com/omnideconv/immunedeconv)：R 包位于
  仓库根目录，仓库自己的 `Remotes` 依赖由 pak 继续解析；
- [BayesPrism 官方仓库](https://github.com/Danko-Lab/BayesPrism)：R 包位于
  `BayesPrism/` 子目录，因此 pak 规格必须保留最后一个 `/BayesPrism`。

为可复现分析固定 commit 时，在规格末尾添加 `@<commit>`，例如
`MuSiC=github::xuranw/MuSiC@<commit>`。三个包应分别记录 commit，不要假定
它们共享版本号或发布时间。

immunedeconv 本身安装成功不代表所有方法都没有额外条件。其官方文档说明
CIBERSORT/CIBERSORT-ABS 需要用户取得许可并提供 CIBERSORT 脚本和签名矩阵；
这部分资源不会由 pak 或 `bulkMAE` 自动下载。

## 4. 开发和发布工具

构建 vignette、运行测试和执行 Bioconductor 检查所需的工具也通过 pak 安装：

```r
pak::pkg_install(c(
  "cran::devtools",
  "cran::knitr",
  "cran::rmarkdown",
  "cran::testthat",
  "bioc::BiocCheck"
))
```

安装后可先查看依赖树和实际版本，再运行发布检查：

```r
pak::pkg_deps_tree(".", dependencies = TRUE)
pak::pkg_status()
```

## 5. 安装之外仍需准备的资源

以下功能即使依赖包已安装，也可能需要本地数据、凭据或运行时联网：

- `import_tximport()` 需要 Salmon、kallisto 或 RSEM 等工具已经生成的定量文件；
- `annotate_biomart()` 使用 Ensembl 服务；
- `annotate_gene_lengths()` 使用 Ensembl，并且返回的是注释转录本长度汇总，
  不是定量工具给出的样本特异 effective length；
- `gene_sets_msigdb()` 首次取得指定 MSigDB 集合时可能下载并缓存资源；
- KEGG、STRINGdb 和通过 OmniPathR 获取的 decoupleR 资源会访问远程服务或
  下载缓存；
- MuSiC 和 BayesPrism 需要与 bulk 数据匹配的单细胞参考；
- `drug_lincs()` 可用 `cmap`、`lincs`、`lincs2` 等标识让 signatureSearch
  通过 ExperimentHub 下载并缓存参考库，也可使用显式本地数据库；预构建库
  要求 human Entrez ID，可用 `drug_lincs_databases()` 查询；
- CIBERSORT 系列方法需要其许可范围内的外部文件。

默认离线测试不应以这些服务当时可用为前提。在线或授权资源测试应独立运行，
并在发布记录中注明日期、资源版本和是否真实执行。

仓库测试默认跳过会访问远程服务的用例。仅在确认允许联网且相关许可、凭据或
资源已经就绪时，才在当前 R 会话中显式启用：

```r
Sys.setenv(BULKMAE_RUN_ONLINE_TESTS = "true")
devtools::test()
Sys.unsetenv("BULKMAE_RUN_ONLINE_TESTS")
```

该开关只表示主动允许在线测试，不会替用户准备授权文件或保证第三方服务可用。
