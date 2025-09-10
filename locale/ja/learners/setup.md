---
title: Setup
---

RとRStudioを最新バージョンでコンピュータにインストールしていることを確認してください。
インストール手順の詳細については、[R入門](https://carpentries-incubator.github.io/bioc-intro/#r-and-rstudio)
コースの[「すでにRとRStudioをインストール済みの場合」](https://carpentries-incubator.github.io/bioc-intro/#r-and-rstudio)
セクションを参照してください。

さらに、本レッスンで使用する以下のパッケージもインストールする必要があります。

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("tidyverse", "SummarizedExperiment",
                       "ExploreModelMatrix", "AnnotationDbi", "org.Hs.eg.db", 
                       "org.Mm.eg.db", "csoneson/ConfoundingExplorer",
                       "DESeq2", "vsn", "ComplexHeatmap", "hgu95av2.db",
                       "RColorBrewer", "hexbin", "cowplot", "iSEE",
                       "clusterProfiler", "enrichplot", "kableExtra",
                       "msigdbr", "gplots", "ggplot2", "simplifyEnrichment",
                       "apeglm", "microbenchmark", "Biostrings",
                       "SingleCellExperiment"))

```

_ワークショップに参加される方は、開始前までに上記の項目をすべて完了してください。ご不明な点がございましたら、ワークショップ開始30分前から講師が対応いたします。_







