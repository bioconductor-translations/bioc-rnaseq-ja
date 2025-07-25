---
title: "Rに量的データをインポートして注釈を付ける"
source: Rmd
teaching: 80
output:
  html_document:
    df_print: ページ付
exercises: 40
---





::::::::::::::::::::::::::::::::::::::: objectives

- 量的データをSummarizedExperimentオブジェクトにインポートする方法を学びます。
- オブジェクトに追加の遺伝子注釈を追加する方法を学びます。

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- 量的遺伝子発現データをRで下流の統計分析に適したオブジェクトにインポートするにはどうすればよいですか？
- 通常使用される遺伝子識別子の種類は何ですか？ それらのマッピングはどのように行われますか？

::::::::::::::::::::::::::::::::::::::::::::::::::

## パッケージを読み込む

このエピソードでは、アドオンRパッケージの関数をいくつか使用します。 それらを使用するためには、`library`から読み込む必要があります。


``` r
suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    library(hgu95av2.db)
    library(SummarizedExperiment)
})
```

`there is no package called 'XXXX'`というエラーメッセージが表示された場合、これはこのバージョンのRに対して、パッケージがまだインストールされていないことを意味します。このワークショップに必要なすべてのパッケージをインストールするには、[Summary and Setup](https://carpentries-incubator.github.io/bioc-rnaseq/index.html)の下部を参照してください。 インストールする必要がある場合は、上記の`library`コマンドを再実行してそれらを読み込むことを忘れないでください。

## データを読み込む

前回のエピソードでは、Rを使用してインターネットから4つのファイルをダウンロードし、それらをコンピュータに保存しました。 しかし、これらのファイルはまだRに読み込まれていないため、作業することができません。 [Blackmore et al. 2017](https://pubmed.ncbi.nlm.nih.gov/28696309/)の元の実験デザインはかなり複雑でした。8週齢のオスとメスのC57BL/6マウスは、インフルエンザ感染前の0日目、感染後の4日目および8日目に収集されました。 各マウスからは、小脳と脊髄の組織がRNA-seq用に取り出されました。 元々は「性別 x 時間 x 組織」群御に4匹のマウスがいましたが、途中でいくつかが失われ、合計で45のサンプルに至りました。 このワークショップでは、分析を簡素化するために22の小脳サンプルのみを使用します。 発現の定量化は、STARを使用してマウスのゲノムにアライメントを行い、その後、遺伝子にマッピングされるリードの数を数えることを通じて行われました。 各遺伝子ごとのサンプルあたりのカウントに加えて、どのサンプルがどの性別/時間点/複製に属するかの情報も必要です。 遺伝子に関しては、注釈と呼ばれる追加の情報があると便利です。
前回ダウンロードしたデータファイルを読み込み、探索し始めましょう：

### カウント


``` r
counts <- read.csv("data/GSE96870_counts_cerebellum.csv", 
                   row.names = 1)
dim(counts)
```

``` output
[1] 41786    22
```

``` r
# View(counts)
```

遺伝子は行に、サンプルは列に含まれています。したがって、41,786の遺伝子と22のサンプルのカウントがあります。 `View()`コマンドはウェブサイト用にコメントアウトされていますが、実行するとRStudioでデータを確認したり、特定の列でテーブルを並べ替えたりできます。 ただし、ビューワは`counts`オブジェクト内部のデータを変更できないため、見るだけで、永久に並べ替えたり編集したりすることはできません。 終了したら、タブのXを使ってビューワを閉じます。 行名は遺伝子シンボルであり、列名はGEOのサンプルIDであるようです。これは、私たちがどのサンプルが何かを教えてくれないため、あまり有益ではありません。

### サンプルの注釈

次に、サンプルの注釈を読み込みます。 カウント行列の列にはサンプルが含まれているため、オブジェクトに`coldata`という名前を付けます：


``` r
coldata <- read.csv("data/GSE96870_coldata_cerebellum.csv",
                    row.names = 1)
dim(coldata)
```

``` output
[1] 22 10
```

``` r
# View(coldata)
```

今、サンプルが行にあり、GEOサンプルIDが行名として付いています。そして、私たちには10列の情報があります。 このワークショップで最も便利な列は、`geo_accession`（再度、GEOサンプルID）、`sex`、および`time`です。

### 遺伝子の注釈

カウントには遺伝子シンボルしかありませんが、それは短くて人間の脳にとってはある程度認識可能ですが、実際に測定された遺伝子の正確な同定子としてはあまり利便性がありません。 そのため、著者によって提供された追加の遺伝子注釈が必要です。 `count`および`coldata`ファイルはカンマ区切り値（.csv）形式でしたが、遺伝子注釈ファイルにはそれが使用できません。なぜなら、説明には、カンマを含む可能性があるため、.csvファイルを正しく読み込むのを妨げるからです。 その代わりに、遺伝子注釈ファイルはタブ区切り値（.tsv）形式です。 同様に、説明には単一引用符`'`（例：5'）が含まれる可能性があり、Rはデフォルトでこれを文字列として扱うためです。 そのため、データがタブ区切りであることを指定するために、より一般的な関数`read.delim()`に追加の引数を使用する必要があります（`sep = "\t"`）で、引用符は使用しない（`quote = ""`）。 さらに、最初の行に列名が含まれていることを指定するためにその他の引数を追加し（`header = TRUE`）、行名として指定される遺伝子シンボルは5列目（`row.names = 5`）であること、NCBIの種特異的遺伝子ID（すなわちENTREZID）は、数字のように見えるが文字列として読み込む（`colClasses`引数）必要があることを指定します。 使用可能な引数に関する詳細については、関数名の先頭にクエスチョンマークを入力することで確認できます。 (例：`?read.delim`)


``` r
rowranges <- read.delim("data/GSE96870_rowranges.tsv", 
                        sep = "\t", 
                        colClasses = c(ENTREZID = "character"),
                        header = TRUE, 
                        quote = "", 
                        row.names = 5)
dim(rowranges)
```

``` output
[1] 41786     7
```

``` r
# View(rowranges)
```

41,786の遺伝子ごとに、`seqnames`（例えば、染色体数）、`start`および`end`位置、`strand`、`ENTREZID`、遺伝子産物説明（`product`）および特徴タイプ（`gbkey`）があります。 これらの遺伝子レベルのメタデータは、下流分析に役立ちます。 たとえば、`gbkey`列から、どのような種類の遺伝子があり、それらがデータセットにどのくらい含まれているかを確認できます：


``` r
table(rowranges$gbkey)
```

``` output

     C_region     D_segment          exon     J_segment      misc_RNA 
           20            23          4008            94          1988 
         mRNA         ncRNA precursor_RNA          rRNA          tRNA 
        21198         12285          1187            35           413 
    V_segment 
          535 
```

:::::::::::::::::::::::::::::::::::::::  challenge

## チャレンジ: 以下のポイントを隣の人と話し合ってください。

1. `counts`、`coldata`、`rowranges`の3つのオブジェクトは、行および列に関してどのように関連していますか？
2. mRNA遺伝子のみを分析したい場合、一般的にはどのようにしてそれらだけを保持しますか？（正確なコードではない）
3. 最初の2つのサンプルが外れ値であると印象付ける場合、それらを削除するにはどうすればよいですか？（一般的には、正確なコードではない）

::::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::: solution

1. `counts`では、行は遺伝子であり、`rowranges`の行と同じです。 `counts`の列はサンプルですが、これは`coldata`の行に対応します。
2. `counts`の行をmRNA遺伝子に限定し、`rowranges`の行もそのようにしなければなりません。
3. `counts`の列と`coldata`の行で、最初の2つのサンプルを除外するために両方の行をサブセットする必要があります。

:::::::::::::::::::::::::::::::::::

関連情報を別のオブジェクトに保持することで、カウント、遺伝子の注釈、サンプルの注釈の間で不一致が発生する可能性があります。 これが、Bioconductorが`SummarizedExperiment`という特殊なS4クラスを作成した理由です。 `SummarizedExperiment`の詳細は、[RNA解析とBioconductorの導入](https://carpentries-incubator.github.io/bioc-intro/60-next-steps.html#next-steps)ワークショップの最後で詳しく説明されています。
リマインダーとして、`SummarizedExperiment`クラスの構造を表す図を見てみましょう：

<img src="https://uclouvain-cbio.github.io/WSBIM1322/figs/SE.svg" alt="Schematic showing the composition of a SummarizedExperiment object, with three assay matrices of equal dimension, rowData with feature annotations, colData with sample annotations, and a metadata list." width="80%" style="display: block; margin: auto;" />

これは、任意の種類の定量的なオミクスデータ（`assays`）と、それにリンクされたサンプル注釈（`colData`）、および（遺伝子）特徴注釈（`rowRanges`）または染色体、開始および終了位置を持たない（`rowData`）形式で保持されるように設計されています。 これらの3つのテーブルが（正しく）リンクされると、サンプルや特徴の部分集合が`assay`、`colData`、`rowRanges`の正しい部分集合に変わります。 さらに、ほとんどのBioconductorパッケージは同じコアデータインフラストラクチャに基づいて構築されているため、`SummarizedExperiment`オブジェクトを認識し、操作することができます。 さらに、ほとんどのBioconductorパッケージは同じコアデータインフラストラクチャの周りに構築されているため、`SummarizedExperiment`オブジェクトを認識し、操作できるようになります。 最も人気のある2つのRNA-seq統計分析パッケージは、統計結果用に追加のスロットがある`SummarizedExperiment`に類似した独自の拡張S4クラスを持っています：[DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#the-deseqdataset)の`DESeqDataSet`および[edgeR](https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/DGEList-class)の`DGEList`です。 統計分析に使用するものが何であれ、データを`SummarizedExperiment`に入れることから始めることができます。

## SummarizedExperimentを組み立てる

これらのオブジェクトから`SummarizedExperiment`を作成します。

- `count`オブジェクトは`assays`スロットに保存されます。
- サンプル情報を持つ`coldata`オブジェクトは、`colData`スロットに保存されます（_**サンプルメタデータ**_）
- 遺伝子を記述する`rowranges`オブジェクトは、`rowRanges`スロットに保存されます（_**特徴メタデータ**_）

それらを組み合わせる前に、サンプルと遺伝子が同じ順序であることを絶対に確認する必要があります！ `count`と`coldata`が同じ数のサンプルを持っていること、また`count`と`rowranges`が同じ数の遺伝子を持っていることはわかりましたが、同じ順序になっているかどうかを明示的に確認することはしていませんでした。 確認する簡単な方法：


``` r
all.equal(colnames(counts), rownames(coldata)) # samples
```

``` output
[1] TRUE
```

``` r
all.equal(rownames(counts), rownames(rowranges)) # genes
```

``` output
[1] TRUE
```

``` r
# 最初がTRUEでない場合は、このようにしてカウントのサンプル/列をコレクトします（これは最初がTRUEでも実行しても構いません）：

tempindex <- match(colnames(counts), rownames(coldata))
coldata <- coldata[tempindex, ]

# 再確認します:
all.equal(colnames(counts), rownames(coldata)) 
```

``` output
[1] TRUE
```

:::::::::::::::::::::::::::::::::::::::  challenge

アッセイ（例：`counts`）および遺伝子注釈テーブル（例：`rowranges`）内の特徴（すなわち遺伝子）が異なる場合、これらをどのように修正できますか？
コードを記述してください。

:::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::: solution


``` r
tempindex <- match(rownames(counts), rownames(rowranges))
rowranges <- rowranges[tempindex, ]

all.equal(rownames(counts), rownames(rowranges)) 
```

:::::::::::::::::::::::::::::::::::

サンプルと遺伝子が同じ順序になっていることを確認したら、`SummarizedExperiment`オブジェクトを作成します。


``` r
# 最後の確認:
stopifnot(rownames(rowranges) == rownames(counts), # features
          rownames(coldata) == colnames(counts)) # samples

se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    rowRanges = as(rowranges, "GRanges"),
    colData = coldata
)
```

遺伝子とサンプルが一致していることが非常に重要であるため、`SummarizedExperiment()`コンストラクタは内部で一致する遺伝子/サンプル数とサンプル/行名が一致することをチェックします。 そうでない場合、いくつかのエラーメッセージが表示されます：


``` r
# サンプル数の誤り:

bad1 <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    rowRanges = as(rowranges, "GRanges"),
    colData = coldata[1:3,]
)
```

``` error
Error in validObject(.Object): invalid class "SummarizedExperiment" object: 
    nb of cols in 'assay' (22) must equal nb of rows in 'colData' (3)
```


``` r
# 同じ数の遺伝子ですが異なる順序:

bad2 <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts)),
  rowRanges = as(rowranges[c(2:nrow(rowranges), 1),], "GRanges"),
  colData = coldata
)
```

``` error
Error in SummarizedExperiment(assays = list(counts = as.matrix(counts)), : the rownames and colnames of the supplied assay(s) must be NULL or identical
  to those of the RangedSummarizedExperiment object (or derivative) to
  construct
```

`SummarizedExperiment`のさまざまなデータスロットにアクセスする方法と、いくつかの操作を行う方法の簡単な概要：


``` r
# カウントにアクセス
head(assay(se))
```

``` output
             GSM2545336 GSM2545337 GSM2545338 GSM2545339 GSM2545340 GSM2545341
Xkr4               1891       2410       2159       1980       1977       1945
LOC105243853          0          0          1          4          0          0
LOC105242387        204        121        110        120        172        173
LOC105242467         12          5          5          5          2          6
Rp1                   2          2          0          3          2          1
Sox17               251        239        218        220        261        232
             GSM2545342 GSM2545343 GSM2545344 GSM2545345 GSM2545346 GSM2545347
Xkr4               1757       2235       1779       1528       1644       1585
LOC105243853          1          3          3          0          1          3
LOC105242387        177        130        131        160        180        176
LOC105242467          3          2          2          2          1          2
Rp1                   3          1          1          2          2          2
Sox17               179        296        233        271        205        230
             GSM2545348 GSM2545349 GSM2545350 GSM2545351 GSM2545352 GSM2545353
Xkr4               2275       1881       2584       1837       1890       1910
LOC105243853          1          0          0          1          1          0
LOC105242387        161        154        124        221        272        214
LOC105242467          2          4          7          1          3          1
Rp1                   3          6          5          3          5          1
Sox17               302        286        325        201        267        322
             GSM2545354 GSM2545362 GSM2545363 GSM2545380
Xkr4               1771       2315       1645       1723
LOC105243853          0          1          0          1
LOC105242387        124        189        223        251
LOC105242467          4          2          1          4
Rp1                   3          3          1          0
Sox17               273        197        310        246
```

``` r
dim(assay(se))
```

``` output
[1] 41786    22
```

``` r
# 上記は、私たちが今持っているのが1つのアッセイ、"counts"のために機能しています。
# しかし、アッセイが複数ある場合は、指定する必要があります。
# 例えば、

head(assay(se, "counts"))
```

``` output
             GSM2545336 GSM2545337 GSM2545338 GSM2545339 GSM2545340 GSM2545341
Xkr4               1891       2410       2159       1980       1977       1945
LOC105243853          0          0          1          4          0          0
LOC105242387        204        121        110        120        172        173
LOC105242467         12          5          5          5          2          6
Rp1                   2          2          0          3          2          1
Sox17               251        239        218        220        261        232
             GSM2545342 GSM2545343 GSM2545344 GSM2545345 GSM2545346 GSM2545347
Xkr4               1757       2235       1779       1528       1644       1585
LOC105243853          1          3          3          0          1          3
LOC105242387        177        130        131        160        180        176
LOC105242467          3          2          2          2          1          2
Rp1                   3          1          1          2          2          2
Sox17               179        296        233        271        205        230
             GSM2545348 GSM2545349 GSM2545350 GSM2545351 GSM2545352 GSM2545353
Xkr4               2275       1881       2584       1837       1890       1910
LOC105243853          1          0          0          1          1          0
LOC105242387        161        154        124        221        272        214
LOC105242467          2          4          7          1          3          1
Rp1                   3          6          5          3          5          1
Sox17               302        286        325        201        267        322
             GSM2545354 GSM2545362 GSM2545363 GSM2545380
Xkr4               1771       2315       1645       1723
LOC105243853          0          1          0          1
LOC105242387        124        189        223        251
LOC105242467          4          2          1          4
Rp1                   3          3          1          0
Sox17               273        197        310        246
```

``` r
# サンプル注釈にアクセス
colData(se)
```

``` output
DataFrame with 22 rows and 10 columns
                     title geo_accession     organism         age         sex
               <character>   <character>  <character> <character> <character>
GSM2545336 CNS_RNA-seq_10C    GSM2545336 Mus musculus     8 weeks      Female
GSM2545337 CNS_RNA-seq_11C    GSM2545337 Mus musculus     8 weeks      Female
GSM2545338 CNS_RNA-seq_12C    GSM2545338 Mus musculus     8 weeks      Female
GSM2545339 CNS_RNA-seq_13C    GSM2545339 Mus musculus     8 weeks      Female
GSM2545340 CNS_RNA-seq_14C    GSM2545340 Mus musculus     8 weeks        Male
...                    ...           ...          ...         ...         ...
GSM2545353  CNS_RNA-seq_3C    GSM2545353 Mus musculus     8 weeks      Female
GSM2545354  CNS_RNA-seq_4C    GSM2545354 Mus musculus     8 weeks        Male
GSM2545362  CNS_RNA-seq_5C    GSM2545362 Mus musculus     8 weeks      Female
GSM2545363  CNS_RNA-seq_6C    GSM2545363 Mus musculus     8 weeks        Male
GSM2545380  CNS_RNA-seq_9C    GSM2545380 Mus musculus     8 weeks      Female
             infection      strain        time      tissue     mouse
           <character> <character> <character> <character> <integer>
GSM2545336  InfluenzaA     C57BL/6        Day8  Cerebellum        14
GSM2545337 NonInfected     C57BL/6        Day0  Cerebellum         9
GSM2545338 NonInfected     C57BL/6        Day0  Cerebellum        10
GSM2545339  InfluenzaA     C57BL/6        Day4  Cerebellum        15
GSM2545340  InfluenzaA     C57BL/6        Day4  Cerebellum        18
...                ...         ...         ...         ...       ...
GSM2545353 NonInfected     C57BL/6        Day0  Cerebellum         4
GSM2545354 NonInfected     C57BL/6        Day0  Cerebellum         2
GSM2545362  InfluenzaA     C57BL/6        Day4  Cerebellum        20
GSM2545363  InfluenzaA     C57BL/6        Day4  Cerebellum        12
GSM2545380  InfluenzaA     C57BL/6        Day8  Cerebellum        19
```

``` r
dim(colData(se))
```

``` output
[1] 22 10
```

``` r
# 遺伝子注釈にアクセス
head(rowData(se))
```

``` output
DataFrame with 6 rows and 3 columns
                ENTREZID                product       gbkey
             <character>            <character> <character>
Xkr4              497097 X Kell blood group p..        mRNA
LOC105243853   105243853 uncharacterized LOC1..       ncRNA
LOC105242387   105242387 uncharacterized LOC1..       ncRNA
LOC105242467   105242467 lipoxygenase homolog..        mRNA
Rp1                19888 retinitis pigmentosa..        mRNA
Sox17              20671 SRY (sex determining..        mRNA
```

``` r
dim(rowData(se))
```

``` output
[1] 41786     3
```

``` r
# 性別、時間、マウスIDを表示するためのより良いサンプルIDを作成します：

se$Label <- paste(se$sex, se$time, se$mouse, sep = "_")
se$Label
```

``` output
 [1] "Female_Day8_14" "Female_Day0_9"  "Female_Day0_10" "Female_Day4_15"
 [5] "Male_Day4_18"   "Male_Day8_6"    "Female_Day8_5"  "Male_Day0_11"  
 [9] "Female_Day4_22" "Male_Day4_13"   "Male_Day8_23"   "Male_Day8_24"  
[13] "Female_Day0_8"  "Male_Day0_7"    "Male_Day4_1"    "Female_Day8_16"
[17] "Female_Day4_21" "Female_Day0_4"  "Male_Day0_2"    "Female_Day4_20"
[21] "Male_Day4_12"   "Female_Day8_19"
```

``` r
colnames(se) <- se$Label

# サンプルは性別と時間に基づいて並んでいません。
se$Group <- paste(se$sex, se$time, sep = "_")
se$Group
```

``` output
 [1] "Female_Day8" "Female_Day0" "Female_Day0" "Female_Day4" "Male_Day4"  
 [6] "Male_Day8"   "Female_Day8" "Male_Day0"   "Female_Day4" "Male_Day4"  
[11] "Male_Day8"   "Male_Day8"   "Female_Day0" "Male_Day0"   "Male_Day4"  
[16] "Female_Day8" "Female_Day4" "Female_Day0" "Male_Day0"   "Female_Day4"
[21] "Male_Day4"   "Female_Day8"
```

``` r
# これを順序を保持するファクターデータに変更し、seオブジェクトを再配置します：

se$Group <- factor(se$Group, levels = c("Female_Day0","Male_Day0", 
                                        "Female_Day4","Male_Day4",
                                        "Female_Day8","Male_Day8"))
se <- se[, order(se$Group)]
colData(se)
```

``` output
DataFrame with 22 rows and 12 columns
                         title geo_accession     organism         age
                   <character>   <character>  <character> <character>
Female_Day0_9  CNS_RNA-seq_11C    GSM2545337 Mus musculus     8 weeks
Female_Day0_10 CNS_RNA-seq_12C    GSM2545338 Mus musculus     8 weeks
Female_Day0_8  CNS_RNA-seq_27C    GSM2545348 Mus musculus     8 weeks
Female_Day0_4   CNS_RNA-seq_3C    GSM2545353 Mus musculus     8 weeks
Male_Day0_11   CNS_RNA-seq_20C    GSM2545343 Mus musculus     8 weeks
...                        ...           ...          ...         ...
Female_Day8_16  CNS_RNA-seq_2C    GSM2545351 Mus musculus     8 weeks
Female_Day8_19  CNS_RNA-seq_9C    GSM2545380 Mus musculus     8 weeks
Male_Day8_6    CNS_RNA-seq_17C    GSM2545341 Mus musculus     8 weeks
Male_Day8_23   CNS_RNA-seq_25C    GSM2545346 Mus musculus     8 weeks
Male_Day8_24   CNS_RNA-seq_26C    GSM2545347 Mus musculus     8 weeks
                       sex   infection      strain        time      tissue
               <character> <character> <character> <character> <character>
Female_Day0_9       Female NonInfected     C57BL/6        Day0  Cerebellum
Female_Day0_10      Female NonInfected     C57BL/6        Day0  Cerebellum
Female_Day0_8       Female NonInfected     C57BL/6        Day0  Cerebellum
Female_Day0_4       Female NonInfected     C57BL/6        Day0  Cerebellum
Male_Day0_11          Male NonInfected     C57BL/6        Day0  Cerebellum
...                    ...         ...         ...         ...         ...
Female_Day8_16      Female  InfluenzaA     C57BL/6        Day8  Cerebellum
Female_Day8_19      Female  InfluenzaA     C57BL/6        Day8  Cerebellum
Male_Day8_6           Male  InfluenzaA     C57BL/6        Day8  Cerebellum
Male_Day8_23          Male  InfluenzaA     C57BL/6        Day8  Cerebellum
Male_Day8_24          Male  InfluenzaA     C57BL/6        Day8  Cerebellum
                   mouse          Label       Group
               <integer>    <character>    <factor>
Female_Day0_9          9  Female_Day0_9 Female_Day0
Female_Day0_10        10 Female_Day0_10 Female_Day0
Female_Day0_8          8  Female_Day0_8 Female_Day0
Female_Day0_4          4  Female_Day0_4 Female_Day0
Male_Day0_11          11   Male_Day0_11 Male_Day0  
...                  ...            ...         ...
Female_Day8_16        16 Female_Day8_16 Female_Day8
Female_Day8_19        19 Female_Day8_19 Female_Day8
Male_Day8_6            6    Male_Day8_6 Male_Day8  
Male_Day8_23          23   Male_Day8_23 Male_Day8  
Male_Day8_24          24   Male_Day8_24 Male_Day8  
```

``` r
# 最後に、プロット内での順序を維持するためにLabel列もファクタにします:

se$Label <- factor(se$Label, levels = se$Label)
```

:::::::::::::::::::::::::::::::::::::::  challenge

1. `Infection`変数の各レベルに対して、サンプルは何個ですか？
2. `se_infected`と`se_noninfected`という名前の2つのオブジェクトを作成し、それぞれに感染サンプルと非感染サンプルのみを含む`se`のサブセットを含めます。
  その後、最初の500遺伝子の各オブジェクトの平均発現レベルを計算し、`summary()`関数を使用してこれらの遺伝子に基づく感染と非感染サンプルの発現レベルの分布を調べます。
3. インフルエンザAに感染した雌のマウスのサンプルは何個ありますか？

::::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::: solution


``` r
# 1
table(se$infection)
```

``` output

 InfluenzaA NonInfected 
         15           7 
```

``` r
# 2
se_infected <- se[, se$infection == "InfluenzaA"]
se_noninfected <- se[, se$infection == "NonInfected"]

means_infected <- rowMeans(assay(se_infected)[1:500, ])
means_noninfected <- rowMeans(assay(se_noninfected)[1:500, ])

summary(means_infected)
```

``` output
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.000e+00 1.333e-01 2.867e+00 7.641e+02 3.374e+02 1.890e+04 
```

``` r
summary(means_noninfected)
```

``` output
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.000e+00 1.429e-01 3.143e+00 7.710e+02 3.666e+02 2.001e+04 
```

``` r
# 3
ncol(se[, se$sex == "Female" & se$infection == "InfluenzaA" & se$time == "Day8"])
```

``` output
[1] 4
```

:::::::::::::::::::::::::::::::::::

## SummarizedExperimentを保存する

これが、私たちの`SummarizedExperiment`オブジェクトを作成するための少しのコードと時間でした。 ワークショップ全体で使用し続ける必要があるため、Rのメモリに戻すためにコンピュータ上の実際の単一ファイルとして保存することが有用です。 Rに特有のファイルを保存するには、`saveRDS()`関数を使用し、後で`readRDS()`関数を使用して再び読み込むことができます。


``` r
saveRDS(se, "data/GSE96870_se.rds")
rm(se) # オブジェクトを削除！
se <- readRDS("data/GSE96870_se.rds")
```

## データの由来と再現性

これで、RNA-SeqデータをRにインポートしてさまざまなパッケージによる分析で使用可能な形式の外部.rdsファイルを作成しました。 ただし、インターネットからダウンロードした3つのファイルから.rdsファイルを作成するために使用したコードの記録を保持する必要があります。 ファイルの由来はどうなっていますか？ つまり、それらはどこから来ており、どのように作成されたのですか？ 元々のカウントおよび遺伝子情報は、GEO公開データベースに預けられました。アクセッション番号は[GSE96870](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96870)です。 ただし、これらのカウントは、配列ベースコールや品質スコアを保持するファストqファイルで配列合わせ/定量化プログラムを実行することによって生成されたものであり、これらは特定のライブラリ調製法を使用して抽出されたRNAから収集したサンプルで生成されたものです。 ふぅ！

元の実験を実施した場合、理想的にはデータが生成された場所と方法の完全な記録を持つべきです。 しかし、公開データセットを利用している場合、最善の方法は、どの元のファイルがどこから来たか、そしてそれに対して行ってきた操作の記録を保持することです。 Rコードを使用してすべてを追跡することは、元の入力ファイルから全体の分析を再現可能にする素晴らしい方法です。 得られる正確な結果は、Rのバージョン、アドオンパッケージのバージョン、さらには使用しているオペレーティングシステムによって異なる可能性があるため、`sessionInfo()`を使用してすべての情報を追跡し、出力を記録するようにしてください（レッスンの最後に例を参照）。

:::::::::::::::::::::::::::::::::::::::  challenge

## チャレンジ: mRNA遺伝子をサブセットにする方法

以前は、mRNA遺伝子に対するサブセットを理論的に論じました。 現在、`SummarizedExperiment`オブジェクトを持っているため、`se`を新しいオブジェクト`se_mRNA`にサブセットするためのコードを書くことがはるかに簡単になります。このオブジェクトには、`rowData(se)$gbkey`がmRNAである遺伝子/行のみを含むものです。 コードを書くと、21,198のmRNA遺伝子を正しく取得したかを確認してください。

::::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::: solution


``` r
se_mRNA <- se[rowData(se)$gbkey == "mRNA" , ]
dim(se_mRNA)
```

``` output
[1] 21198    22
```

:::::::::::::::::::::::::::::::::::::::

## 遺伝子注釈

カウントデータを生成する人によっては、追加の遺伝子注釈の適切なファイルがないかもしれません。 遺伝子シンボルやENTREZID、あるいは他のデータベースのIDのみが存在するかもしれません。 遺伝子注釈の特性は、その注釈戦略と情報源によって異なります。 たとえば、RefSeqヒト遺伝子モデル（つまり、NCBIのEntrez）は、さまざまな研究でよくサポートされ、広く使用されています。 UCSC Known Genes データセットは、Swiss-Prot/TrEMBL (UniProt) のタンパク質データと、GenBankからの関連するmRNAデータに基づいており、UCSC Genome Browserの基盤として機能します。 Ensemblの遺伝子は、自動生成されたゲノムアノテーションと手動キュレーションの両方を含んでいます。

Bioconductorでの詳細情報は、[Annotation Workshop](https://jmacdon.github.io/Bioc2022Anno/articles/AnnotationWorkshop.html)の資料で見つけることができます。

Bioconductorには、遺伝子の追加注釈情報を取得するための多くのパッケージや関数があります。 利用可能なリソースについては、[エピソード7 遺伝子セット濃縮解析](https://carpentries-incubator.github.io/bioc-rnaseq/07-gene-set-analysis.html#gene-set-resources)で詳しく説明されています。

ここでは、遺伝子IDマッピング関数の1つである`mapIds`を紹介します：

```
mapIds(annopkg, keys, column, keytype, ..., multiVals)
```

どこで

- _annopkg_は、注釈パッケージです
- _keys_は、私たちが**知っている**IDです
- _column_は、私たちが**望む**値です
- _keytype_は、使用するキーのタイプです


``` r
mapIds(org.Mm.eg.db, keys = "497097", column = "SYMBOL", keytype = "ENTREZID")
```

``` output
'select()' returned 1:1 mapping between keys and columns
```

``` output
497097 
"Xkr4" 
```

`select()`関数とは異なり、`mapIds()`関数は、追加の引数`multiVals`を通じてキーと列の間の1:多のマッピングを処理します。
以下の例では、`hgu95av2.db`パッケージを使用してこの機能を示します。AffymetrixヒトゲノムU95セット注釈データ。


``` r
keys <- head(keys(hgu95av2.db, "ENTREZID"))
last <- function(x){x[[length(x)]]}

mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID")
```

``` output
'select()' returned 1:many mapping between keys and columns
```

``` output
       10       100      1000     10000 100008586     10001 
   "AAC2"    "ADA1"   "ACOGS"    "MPPH"     "AL4"   "ARC33" 
```

``` r
# 1:多のマッピングがある場合、デフォルトの動作は最初の一致を出力することでした。これは、上で定義した関数を使用して最後の一致を取得するように変更できます：

mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID", multiVals = last)
```

``` output
'select()' returned 1:many mapping between keys and columns
```

``` output
       10       100      1000     10000 100008586     10001 
   "NAT2"     "ADA"    "CDH2"    "AKT3" "GAGE12F"    "MED6" 
```

``` r
# または、すべての多くのマッピングを取得することができます：

mapIds(hgu95av2.db, keys = keys, column = "ALIAS", keytype = "ENTREZID", multiVals = "list")
```

``` output
'select()' returned 1:many mapping between keys and columns
```

``` output
$`10`
[1] "AAC2"  "NAT-2" "PNAT"  "NAT2" 

$`100`
[1] "ADA1" "ADA" 

$`1000`
[1] "ACOGS"  "ADHD8"  "ARVD14" "CD325"  "CDHN"   "CDw325" "NCAD"   "CDH2"  

$`10000`
[1] "MPPH"         "MPPH2"        "PKB-GAMMA"    "PKBG"         "PRKBG"       
[6] "RAC-PK-gamma" "RAC-gamma"    "STK-2"        "AKT3"        

$`100008586`
[1] "AL4"     "CT4.7"   "GAGE-7"  "GAGE-7B" "GAGE-8"  "GAGE7"   "GAGE7B" 
[8] "GAGE12F"

$`10001`
[1] "ARC33"     "NY-REN-28" "MED6"     
```

## セッション情報


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
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] hgu95av2.db_3.13.0          org.Hs.eg.db_3.21.0        
 [3] org.Mm.eg.db_3.21.0         AnnotationDbi_1.70.0       
 [5] SummarizedExperiment_1.38.1 Biobase_2.68.0             
 [7] MatrixGenerics_1.20.0       matrixStats_1.5.0          
 [9] GenomicRanges_1.60.0        GenomeInfoDb_1.44.0        
[11] IRanges_2.42.0              S4Vectors_0.46.0           
[13] BiocGenerics_0.54.0         generics_0.1.4             
[15] knitr_1.50                 

loaded via a namespace (and not attached):
 [1] Matrix_1.7-3            bit_4.6.0               jsonlite_2.0.0         
 [4] compiler_4.5.1          BiocManager_1.30.26     renv_1.1.5             
 [7] crayon_1.5.3            blob_1.2.4              Biostrings_2.76.0      
[10] png_0.1-8               fastmap_1.2.0           yaml_2.3.10            
[13] lattice_0.22-7          R6_2.6.1                XVector_0.48.0         
[16] S4Arrays_1.8.1          DelayedArray_0.34.1     GenomeInfoDbData_1.2.14
[19] DBI_1.2.3               rlang_1.1.6             KEGGREST_1.48.1        
[22] cachem_1.1.0            xfun_0.52               bit64_4.6.0-1          
[25] memoise_2.0.1           SparseArray_1.8.0       RSQLite_2.4.1          
[28] cli_3.6.5               grid_4.5.1              vctrs_0.6.5            
[31] evaluate_1.0.4          abind_1.4-8             httr_1.4.7             
[34] pkgconfig_2.0.3         tools_4.5.1             UCSC.utils_1.4.0       
```

::: keypoints

- 使用される遺伝子発現定量ツールによって、出力を`SummarizedExperiment`または`DGEList`オブジェクトに読み込む方法が異なります（多くはBioconductorパッケージで配布されています）。
- EnsemblやEntrez IDなどの安定した遺伝子識別子は、RNA-seq分析全体で主要な識別子として使用されるべきで、解釈を容易にするために遺伝子シンボルを追加する必要があります。
  :::
