---
source: Rmd
title: Exploratory analysis and quality control
teaching: 120
exercises: 60
editor_options:
  chunk_output_type: console
---



::::::::::::::::::::::::::::::::::::::: objectives

- 遺伝子発現マトリックスの解析方法と、一般的な品質管理手順を習得します。
- 探索的解析用のインタラクティブなアプリケーション環境の構築方法を学びます。

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- RNA-seq解析において、探索的解析がなぜ重要なステップなのでしょうか？
- 探索的解析を行う際、生のカウント行列をどのように前処理すべきですか？
- データを表現するために2次元で十分なのでしょうか？

::::::::::::::::::::::::::::::::::::::::::::::::::

## パッケージの読み込み

RStudioを再起動した直後の状態を想定して、このレッスンで使用するパッケージと、前回のレッスンで作成した `SummarizedExperiment` オブジェクトを読み込みます。


``` r
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(DESeq2)
    library(vsn)
    library(ggplot2)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(hexbin)
    library(iSEE)
})
```


``` r
se <- readRDS("data/GSE96870_se.rds")
```

## 発現していない遺伝子の除去

探索的解析は、データの品質管理と理解において非常に重要なプロセスです。
これにより、データ品質の問題、サンプルの取り違え、コンタミネーションなどを検出できるほか、データ中に存在する顕著なパターンを把握することも可能になります。
本エピソードでは、RNA-seqデータに対する探索的解析の代表的な手法であるクラスタリングと主成分分析（PCA）の2つについて解説します。
これらのツールはRNA-seqデータの解析に限定されたものではなく、他の種類のデータ解析にも応用可能です。
ただし、カウントベースのアッセイにはこの種のデータに適用する際に考慮すべき特有の特徴があります。まず第一に、ゲノム上のすべてのマウス遺伝子が小脳サンプルで発現しているわけではありません。遺伝子の発現が検出可能かどうかを判断するための閾値は複数存在しますが、ここでは非常に厳格な基準を採用します。具体的には、全サンプルを通じて遺伝子の総カウント数が5未満の場合、そもそもデータ量が不足しており、有効な解析を行うことができないものとします。


``` r
nrow(se)
```

``` output
[1] 41786
```

``` r
# Remove genes/rows that do not have > 5 total counts 
se <- se[rowSums(assay(se, "counts")) > 5, ]
nrow(se)
```

``` output
[1] 27430
```

:::::::::::::::::::::::::::::::::::::::  challenge

## 課題：このフィルタリングを通過した遺伝子にはどのような特徴があるのか？

前回のエピソードでは、mRNA遺伝子のみに絞り込むサブセット処理について議論しました。今回はさらに、最低限の発現レベルを基準にサブセット処理を行いました。

1. 各種類の遺伝子のうち、フィルタリングを通過した遺伝子の数はそれぞれどれくらいでしょうか？
2. 異なる閾値を用いてフィルタリングを通過した遺伝子数を比較してください。
3. より厳格なフィルタリングを行う場合の利点と欠点は何でしょうか？また、考慮すべき重要な点にはどのようなものがあるでしょうか？

::::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::: solution

1.

``` r
table(rowData(se)$gbkey)
```

``` output

     C_region          exon     J_segment      misc_RNA          mRNA 
           14          1765            14          1539         16859 
        ncRNA precursor_RNA          rRNA          tRNA     V_segment 
         6789           362             2            64            22 
```

2.

``` r
nrow(se)  # represents the number of genes using 5 as filtering threshold
```

``` output
[1] 27430
```

``` r
length(which(rowSums(assay(se, "counts")) > 10))
```

``` output
[1] 25736
```

``` r
length(which(rowSums(assay(se, "counts")) > 20))
```

``` output
[1] 23860
```

3.
Cons: Risk of removing interesting information
Pros: 
 - Not or lowly expressed genes are unlikely to be biological meaningful.
 - Reduces number of statistical tests (multiple testing).
 - More reliable estimation of mean-variance relationship
 
Potential considerations:
 - Is a gene expressed in both groups?
 - How many samples of each group express a gene?

:::::::::::::::::::::::::::::::::::

## ライブラリサイズの差異について

サンプル間で遺伝子に割り当てられた総リード数に差異が生じる場合、これは主に技術的な要因によるものです。実際には、単に遺伝子の生のリードカウントをサンプル間で直接比較し、リード数が多いサンプルほどその遺伝子の発現量が多いと結論づけることはできません。高いリードカウントは、そのサンプルにおける総リード数が全体的に多いことに起因している可能性があるためです。
本節の残りの部分では、「ライブラリサイズ」という用語を、サンプルごとに遺伝子に割り当てられた総リード数を指すものとして使用します。まずすべてのサンプルのライブラリサイズを比較する必要があります。


``` r
# Add in the sum of all counts

se$libSize <-  colSums(assay(se))

# Plot the libSize by using R's native pipe |>
# to extract the colData, turn it into a regular
# data frame then send to ggplot:

colData(se) |>
  as.data.frame() |>
  ggplot(aes(x = Label, y = libSize / 1e6, fill = Group)) + 
         geom_bar(stat = "identity") + theme_bw() + 
         labs(x = "Sample", y = "Total count in millions") + 
         theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
```

<img src="fig/04-exploratory-qc-rendered-lib-size-1.png" alt="Barplot with total count on the y-axis and sample name on the x-axis, with bars colored by the group annotation. The total count varies between approximately 32 and 43 million." style="display: block; margin: auto;" />

サンプル間のライブラリサイズの差異を適切に補正しなければ、誤った結論を導く危険性があります。RNA-seqデータにおいてこの処理を行う標準的な方法は、2段階の手順で説明できます。
まず、サンプルごとに固有の補正係数である「サイズ因子」を推定します。これらの因子を用いて生のカウント値を除算すれば、サンプル間でより比較可能な値が得られるようになります。
次に、これらのサイズ因子をデータの統計解析プロセスに組み込みます。
重要なのは、特定の解析手法においてこの処理がどのように実装されているかを詳細に確認することです。
場合によっては、アナリスト自身がカウント値をサイズ因子で明示的に除算する必要があることもあります。
また別のケースでは（差異発現解析で後述するように）、これらを解析ツールに別途提供することが重要です。ツールはこれらの因子を適切に統計モデルに適用します。

`DESeq2`では、`estimateSizeFactors()`関数を用いてサイズ因子を計算します。
この関数で推定されるサイズ因子は、ライブラリサイズの差異補正と、サンプル間のRNA組成の差異補正を組み合わせたものです。
後者は特に重要です。RNA-seqデータは組成的な性質を持つため、遺伝子間で分配されるリード数には一定の上限があります。もし特定の遺伝子（あるいは少数の非常に高発現遺伝子）がリードの大部分を占有すると、他のすべての遺伝子は必然的に非常に低いカウント値しか得られなくなります。ここで私たちは、内部構造としてこれらのサイズ因子を格納できる`DESeqDataSet`オブジェクトに`SummarizedExperiment`オブジェクトを変換します。また、主要な実験デザイン（性別と時間）を指定する必要もあります。


``` r
dds <- DESeq2::DESeqDataSet(se, design = ~ sex + time)
```

``` warning
Warning in DESeq2::DESeqDataSet(se, design = ~sex + time): some variables in
design formula are characters, converting to factors
```

``` r
dds <- estimateSizeFactors(dds)

# Plot the size factors against library size
# and look for any patterns by group:

ggplot(data.frame(libSize = colSums(assay(dds)),
                  sizeFactor = sizeFactors(dds),
                  Group = dds$Group),
       aes(x = libSize, y = sizeFactor, col = Group)) + 
    geom_point(size = 5) + theme_bw() + 
    labs(x = "Library size", y = "Size factor")
```

<img src="fig/04-exploratory-qc-rendered-est-size-factors-1.png" alt="Scatterplot with library size on the x-axis and size factor on the y-axis, showing a high correlation between the two variables." style="display: block; margin: auto;" />

## データの前処理

探索的解析のための手法に関する研究文献は数多く存在します。
これらの手法の多くは、入力データ（ここでは各遺伝子）の分散が平均値から比較的独立している場合に最も効果的に機能します。
RNA-seqのようなリードカウントデータの場合、この前提は成立しません。
実際には、分散は平均リードカウントの増加に伴って増大する傾向があります。


``` r
meanSdPlot(assay(dds), ranks = FALSE)
```

<img src="fig/04-exploratory-qc-rendered-mean-sd-plot-raw-1.png" alt="Hexagonal heatmap with the mean count on the x-axis and the standard deviation of the count on the y-axis, showing a generally increasing standard deviation with increasing mean. The density of points is highest for low count values." style="display: block; margin: auto;" />

この問題に対処する方法は2つあります：1つ目は、カウントデータに特化した分析手法を開発する方法、2つ目は既存の手法を適用できるようにカウントデータを変換する方法です。
どちらの方法も研究されてきましたが、現時点では後者のアプローチの方が実際に広く採用されています。具体的には、DESeq2の分散安定化変換を用いてデータを変換した後、平均リードカウントと分散の相関関係が適切に除去されていることを確認できます。


``` r
vsd <- DESeq2::vst(dds, blind = TRUE)
meanSdPlot(assay(vsd), ranks = FALSE)
```

<img src="fig/04-exploratory-qc-rendered-mean-sd-plot-vst-1.png" alt="Hexagonal heatmap with the mean variance-stabilized values on the x-axis and the standard deviation of these on the y-axis. The trend is generally flat, with no clear association between the mean and standard deviation." style="display: block; margin: auto;" />

## ヒートマップとクラスタリング解析

発現パターンの類似性に基づいてサンプルをクラスタリングする方法は複数存在します。最も単純な手法の一つは、すべてのサンプルペア間のユークリッド距離を計算することです（距離が長いほど類似性が低いことを示します）。その後、分岐型デンドログラムとヒートマップの両方を用いて結果を可視化し、距離を色で表現します。この解析から、Day 8のサンプルは他のサンプル群と比較して互いにより類似していることが明らかになりました。ただし、Day 4とDay 0のサンプルは明確には分離していません。代わりに、雄個体と雌個体は確実に分離することが確認されました。


``` r
dst <- dist(t(assay(vsd)))
colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
ComplexHeatmap::Heatmap(
    as.matrix(dst), 
    col = colors,
    name = "Euclidean\ndistance",
    cluster_rows = hclust(dst),
    cluster_columns = hclust(dst),
    bottom_annotation = columnAnnotation(
        sex = vsd$sex,
        time = vsd$time,
        col = list(sex = c(Female = "red", Male = "blue"),
                   time = c(Day0 = "yellow", Day4 = "forestgreen", Day8 = "purple")))
)
```

<img src="fig/04-exploratory-qc-rendered-heatmap-1.png" alt="Heatmap of Euclidean distances between all pairs of samples, with hierarchical cluster dendrogram for both rows and columns. Samples from day 8 cluster separately from samples from days 0 and 4. Within days 0 and 4, the main clustering is instead by sex." style="display: block; margin: auto;" />

## PCA

Principal component analysis is a dimensionality reduction method, which projects the samples into a lower-dimensional space.
This lower-dimensional representation can be used for visualization, or as the input for other analysis methods.
The principal components are defined in such a way that they are orthogonal, and that the projection of the samples into the space they span contains as much variance as possible.
It is an _unsupervised_ method in the sense that no external information about the samples (e.g., the treatment condition) is taken into account.
In the plot below we represent the samples in a two-dimensional principal component space.
For each of the two dimensions, we indicate the fraction of the total variance that is represented by that component.
By definition, the first principal component will always represent more of the variance than the subsequent ones.
The fraction of explained variance is a measure of how much of the 'signal' in the data that is retained when we project the samples from the original, high-dimensional space to the low-dimensional space for visualization.


``` r
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("sex", "time"),
                           returnData = TRUE)
```

``` output
using ntop=500 top features by variance
```

``` r
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sex, shape = time), size = 5) +
    theme_minimal() +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() + 
    scale_color_manual(values = c(Male = "blue", Female = "red"))
```

<img src="fig/04-exploratory-qc-rendered-pca-1.png" alt="Scatterplot of samples projected onto the first two principal components, colored by sex and shaped according to the experimental day. The main separation along PC1 is between male and female samples. The main separation along PC2 is between samples from day 8 and samples from days 0 and 4." style="display: block; margin: auto;" />

:::::::::::::::::::::::::::::::::::::::  challenge

## Challenge: Discuss the following points with your neighbour

1. Assume you are mainly interested in expression changes associated with the time after infection (Reminder Day0 -> before infection). What do you need to consider in downstream analysis?

2. Consider an experimental design where you have multiple samples from the same donor. You are still interested in differences by time and observe the following PCA plot. What does this PCA plot suggest?


``` output
using ntop=500 top features by variance
```

<img src="fig/04-exploratory-qc-rendered-pca-exercise-1.png" alt="Scatterplot of samples projected onto the first two principal components, colored by a hypothetical sample ID annotation and shaped according to a hypothetical experimental day annotation. In the plot, samples with the same sample ID tend to cluster together." style="display: block; margin: auto;" />

::::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::: solution

1. The major signal in this data (37% variance) is associated with sex. As we are not interested in sex-specific changes over time, we need to adjust for this in downstream analysis (see [next episodes](../episodes/05-differential-expression.Rmd)) and keep it in mind for further exploratory downstream analysis. A possible way to do so is to remove genes on sex chromosomes.

2.

- A strong donor effect, that needs to be accounted for.
- What does PC1 (37% variance) represent? Looks like 2 donor groups?
- No association of PC1 and PC2 with time --> no or weak transcriptional effect of time
  \--> Check association with higher PCs (e.g., PC3,PC4, ..)

:::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Challenge: Plot the PCA colored by library sizes.

Compare before and after variance stabilizing transformation.

_Hint: The `DESeq2::plotPCA` expect an object of the class `DESeqTransform` as input. You can transform a `SummarizedExperiment` object using `plotPCA(DESeqTransform(se))`_

::::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::: solution


``` r
pcaDataVst <- DESeq2::plotPCA(vsd, intgroup = c("libSize"),
                              returnData = TRUE)
```

``` output
using ntop=500 top features by variance
```

``` r
percentVar <- round(100 * attr(pcaDataVst, "percentVar"))
ggplot(pcaDataVst, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = libSize / 1e6), size = 5) +
    theme_minimal() +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() + 
    scale_color_continuous("Total count in millions", type = "viridis")
```

<img src="fig/04-exploratory-qc-rendered-pca-lib-1.png" alt="Scatterplot of samples projected onto the first two principal components of the variance-stabilized data, colored by library size. The library sizes are between approximately 32.5 and 42.5 million. There is no strong association between the library sizes and the principal components." style="display: block; margin: auto;" />


``` r
pcaDataCts <- DESeq2::plotPCA(DESeqTransform(se), intgroup = c("libSize"),
                              returnData = TRUE)
```

``` output
using ntop=500 top features by variance
```

``` r
percentVar <- round(100 * attr(pcaDataCts, "percentVar"))
ggplot(pcaDataCts, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = libSize / 1e6), size = 5) +
    theme_minimal() +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() + 
    scale_color_continuous("Total count in millions", type = "viridis")
```

<img src="fig/04-exploratory-qc-rendered-pca-lib-vst-1.png" alt="Scatterplot of samples projected onto the first two principal components of the count matrix, colored by library size. The library sizes are between approximately 32.5 and 42.5 million. The first principal component is strongly correlated with the library size." style="display: block; margin: auto;" />

:::::::::::::::::::::::::::::::::::

## Interactive exploratory data analysis

Often it is useful to look at QC plots in an interactive way to directly explore different experimental factors or get insides from someone without coding experience.
Useful tools for interactive exploratory data analysis for RNA-seq are [Glimma](https://bioconductor.org/packages/release/bioc/html/Glimma.html) and [iSEE](https://bioconductor.org/packages/release/bioc/html/iSEE.html)

:::::::::::::::::::::::::::::::::::::::  challenge

## Challenge: Interactively explore our data using iSEE


``` r
## Convert DESeqDataSet object to a SingleCellExperiment object, in order to 
## be able to store the PCA representation
sce <- as(dds, "SingleCellExperiment")

## Add PCA to the 'reducedDim' slot
stopifnot(rownames(pcaData) == colnames(sce))
reducedDim(sce, "PCA") <- as.matrix(pcaData[, c("PC1", "PC2")])

## Add variance-stabilized data as a new assay
stopifnot(colnames(vsd) == colnames(sce))
assay(sce, "vsd") <- assay(vsd)

app <- iSEE(sce)
shiny::runApp(app)
```

::::::::::::::::::::::::::::::::::::::::::::::::::

## Session info


``` r
sessionInfo()
```

``` output
R version 4.5.1 (2025-06-13)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] iSEE_2.20.0                 SingleCellExperiment_1.30.1
 [3] hexbin_1.28.5               RColorBrewer_1.1-3         
 [5] ComplexHeatmap_2.24.1       ggplot2_3.5.2              
 [7] vsn_3.76.0                  DESeq2_1.48.1              
 [9] SummarizedExperiment_1.38.1 Biobase_2.68.0             
[11] MatrixGenerics_1.20.0       matrixStats_1.5.0          
[13] GenomicRanges_1.60.0        GenomeInfoDb_1.44.0        
[15] IRanges_2.42.0              S4Vectors_0.46.0           
[17] BiocGenerics_0.54.0         generics_0.1.4             

loaded via a namespace (and not attached):
 [1] rlang_1.1.6             magrittr_2.0.3          shinydashboard_0.7.3   
 [4] clue_0.3-66             GetoptLong_1.0.5        compiler_4.5.1         
 [7] mgcv_1.9-3              png_0.1-8               vctrs_0.6.5            
[10] pkgconfig_2.0.3         shape_1.4.6.1           crayon_1.5.3           
[13] fastmap_1.2.0           XVector_0.48.0          labeling_0.4.3         
[16] promises_1.3.3          shinyAce_0.4.4          UCSC.utils_1.4.0       
[19] preprocessCore_1.70.0   xfun_0.52               cachem_1.1.0           
[22] jsonlite_2.0.0          listviewer_4.0.0        later_1.4.2            
[25] DelayedArray_0.34.1     BiocParallel_1.42.1     parallel_4.5.1         
[28] cluster_2.1.8.1         R6_2.6.1                bslib_0.9.0            
[31] limma_3.64.1            jquerylib_0.1.4         Rcpp_1.0.14            
[34] iterators_1.0.14        knitr_1.50              httpuv_1.6.16          
[37] Matrix_1.7-3            splines_4.5.1           igraph_2.1.4           
[40] tidyselect_1.2.1        abind_1.4-8             yaml_2.3.10            
[43] doParallel_1.0.17       codetools_0.2-20        affy_1.86.0            
[46] miniUI_0.1.2            lattice_0.22-7          tibble_3.3.0           
[49] shiny_1.11.0            withr_3.0.2             evaluate_1.0.4         
[52] circlize_0.4.16         pillar_1.10.2           affyio_1.78.0          
[55] BiocManager_1.30.26     renv_1.1.5              DT_0.33                
[58] foreach_1.5.2           shinyjs_2.1.0           scales_1.4.0           
[61] xtable_1.8-4            glue_1.8.0              tools_4.5.1            
[64] colourpicker_1.3.0      locfit_1.5-9.12         colorspace_2.1-1       
[67] nlme_3.1-168            GenomeInfoDbData_1.2.14 vipor_0.4.7            
[70] cli_3.6.5               viridisLite_0.4.2       S4Arrays_1.8.1         
[73] dplyr_1.1.4             gtable_0.3.6            rintrojs_0.3.4         
[76] sass_0.4.10             digest_0.6.37           SparseArray_1.8.0      
[79] ggrepel_0.9.6           rjson_0.2.23            htmlwidgets_1.6.4      
[82] farver_2.1.2            htmltools_0.5.8.1       lifecycle_1.0.4        
[85] shinyWidgets_0.9.0      httr_1.4.7              GlobalOptions_0.1.2    
[88] statmod_1.5.0           mime_0.13              
```

:::::::::::::::::::::::::::::::::::::::: keypoints

- Exploratory analysis is essential for quality control and to detect potential problems with a data set.
- Different classes of exploratory analysis methods expect differently preprocessed data. The most commonly used methods expect counts to be normalized and log-transformed (or similar- more sensitive/sophisticated), to be closer to homoskedastic. Other methods work directly on the raw counts.

::::::::::::::::::::::::::::::::::::::::::::::::::


