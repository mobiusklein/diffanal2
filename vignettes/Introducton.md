---
title: "Introduction and Overview of Usage"
author: "Joshua Klein"
date: "2014-05-19"
output: html_document
---
<!--
%!\SweaveUTF8
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{Introduction and Overview of Usage}
-->


# Introduction to diffanal2

## What is diffanal2?

Based on the original `diffanal.R` script, `diffanal2` supports performing differential gene expression analysis of microarray data using asymptotic and permutation t-tests, as well as a new feature to support using `limma` - *Linear Models for Microarray Data*. It also supports accepting an arbitrary scoring function, so you can extend it without needing to rewrite the guts of the program.

## A Brief Example

Start by loading a microarray dataset. We'll use a lymphoma sample study with 116 participants.

```r
data(lymphoma.2010)
lymphoma
```

```
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 19081 features, 116 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: MS_D_1001 MS_D_1002 ... old_c514 (116 total)
##   varLabels: SampleID Class ... COO (7 total)
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation: HG_U133plus2
```
Each sample in the study was clinically characterized several times.

```r
head(pData(lymphoma))
```

```
##            SampleID Class meta.1 best.10_13     cc gold.standard  COO
## MS_D_1001 MS_D_1001    HR     HR         HR     HR            HR <NA>
## MS_D_1002 MS_D_1002   BCR    BCR        BCR OxPhos          <NA>  GCB
## MS_D_1003 MS_D_1003   BCR    BCR        BCR    BCR           BCR  GCB
## MS_D_1004 MS_D_1004    HR OxPhos       <NA>     HR          <NA> <NA>
## MS_D_1005 MS_D_1005    HR     HR         HR     HR            HR  ABC
## MS_D_1006 MS_D_1006    HR     HR         HR    BCR          <NA> <NA>
```
For now, let's assume the column **`meta.1`** as being the phenotype of interest for this dataset. While `diffanal2` can operate on a matrix of expression data and a phenotype vector or data.frame, it knows how to unpack an `ExpressionSet` object and will assume me mean for the phenotype to be taken from the `pData` member of the object unless told otherwise. We can specify the predictor phenotype either as the first parameter of a `formula` object passed to the `model.formula` parameter, or as a single string character vector to the `predictor` parameter. 

To start small, we can use a regular t-test for our scoring function. We'll also turn on the verbose option to see each step in the pipeline. This command will verify that the data is clean, without any NA values in the phenotype, removing any samples that do, perform variance filtering using the `mad` statistic to discard genes with extreme variability, keeping the top 5000 (see the `ngenes.filter` parameter to change this), and then partition and test the data.

```r
library('diffanal2')
options(verbose=T)
t.test_meta.1 <- df2.t.test(lymphoma, model.formula = ~ meta.1)
```

```
## Cleaning data, removing samples with NA phenotype data
## 0  samples dropped
## Variation filtering based on mad .. 
## done.
## 
## Selecting top 5000 by mad .. 
## done, 5000 genes selected.
## 
## Calculating summary statistics
```

For a multinomial phenotype like this one, the test will run in one-vs-all mode by default, comparing all of the memebers of one class to all members of all other classes, for each class present in the phenotype selected. The results are returned as an instance of `diffanal.results data.frame`, an S3 class built on top of a data.frame with additional methods.


```r
head(t.test_meta.1)
```

```
##                            gene BCR..Not.BCR..p.value HR..Not.HR..p.value
## ENSG00000121774 ENSG00000121774             4.814e-18           0.1020927
## ENSG00000099783 ENSG00000099783             3.452e-17           0.1338365
## ENSG00000162298 ENSG00000162298             3.411e-17           0.1193124
## ENSG00000129351 ENSG00000129351             1.140e-16           0.0027089
## ENSG00000099246 ENSG00000099246             2.190e-16           0.0007024
## ENSG00000132182 ENSG00000132182             2.783e-16           0.2912822
##                 OxPhos..Not.OxPhos..p.value BCR..Not.BCR..adj.p.value
## ENSG00000121774                   4.897e-07                 2.407e-14
## ENSG00000099783                   4.176e-10                 5.754e-14
## ENSG00000162298                   7.395e-09                 5.754e-14
## ENSG00000129351                   1.380e-07                 1.424e-13
## ENSG00000099246                   3.201e-07                 2.190e-13
## ENSG00000132182                   8.146e-14                 2.319e-13
##                 HR..Not.HR..adj.p.value OxPhos..Not.OxPhos..adj.p.value
## ENSG00000121774                0.143298                       2.898e-06
## ENSG00000099783                0.181843                       6.779e-09
## ENSG00000162298                0.164433                       7.665e-08
## ENSG00000129351                0.005749                       9.480e-07
## ENSG00000099246                0.001751                       1.991e-06
## ENSG00000132182                0.358722                       5.737e-12
##                 BCR..Not.BCR..t.score HR..Not.HR..t.score
## ENSG00000121774                -10.40               1.655
## ENSG00000099783                -10.22               1.513
## ENSG00000162298                -10.12               1.572
## ENSG00000129351                -10.44               3.067
## ENSG00000099246                 10.41              -3.503
## ENSG00000132182                -10.21               1.061
##                 OxPhos..Not.OxPhos..t.score BCR..mean HR..mean
## ENSG00000121774                       5.477    1163.8    836.1
## ENSG00000099783                       6.984    1985.4   1478.2
## ENSG00000162298                       6.433     678.8    377.5
## ENSG00000129351                       5.743     451.0    294.4
## ENSG00000099246                      -5.461     292.3    381.5
## ENSG00000132182                       8.647     250.9    150.5
##                 OxPhos..mean Not.BCR..mean Not.HR..mean Not.OxPhos..mean
## ENSG00000121774        729.8         780.9        944.1           1004.3
## ENSG00000099783       1274.7        1372.5       1625.6           1738.5
## ENSG00000162298        281.7         327.7        477.7            532.1
## ENSG00000129351        271.1         282.3        359.9            374.8
## ENSG00000099246        393.5         387.7        343.5            335.7
## ENSG00000132182        107.4         128.1        178.2            202.1
##                 BCR..std.dev HR..std.dev OxPhos..std.dev Not.BCR..std.dev
## ENSG00000121774       151.07      186.31          187.75           193.35
## ENSG00000099783       336.49      270.04          246.88           276.20
## ENSG00000162298       211.46      187.37          136.48           168.86
## ENSG00000129351        98.66       55.15           67.94            62.82
## ENSG00000099246        39.83       55.70           55.05            55.33
## ENSG00000132182        83.05       51.75           33.63            48.17
##                 Not.HR..std.dev Not.OxPhos..std.dev BCR..median HR..median
## ENSG00000121774          276.45              235.42      1149.6      850.4
## ENSG00000099783          462.08              396.88      1890.6     1530.3
## ENSG00000162298          266.49              249.96       664.8      334.8
## ENSG00000129351          123.47              112.25       452.6      290.0
## ENSG00000099246           69.86               65.63       288.4      374.2
## ENSG00000132182           95.62               85.62       227.5      144.1
##                 OxPhos..median Not..BCR..median Not..HR..median
## ENSG00000121774          740.4            782.4           972.5
## ENSG00000099783         1257.9           1384.9          1598.2
## ENSG00000162298          250.8            291.3           445.1
## ENSG00000129351          274.7            287.3           346.1
## ENSG00000099246          381.5            380.8           338.3
## ENSG00000132182          101.3            119.5           147.4
##                 Not..OxPhos..median BCR..Not.BCR..fold.change
## ENSG00000121774              1005.8                    1.4904
## ENSG00000099783              1703.2                    1.4466
## ENSG00000162298               502.2                    2.0711
## ENSG00000129351               352.9                    1.5976
## ENSG00000099246               328.4                    0.7539
## ENSG00000132182               186.3                    1.9585
##                 HR..Not.HR..fold.change OxPhos..Not.OxPhos..fold.change
## ENSG00000121774                  0.8856                          0.7267
## ENSG00000099783                  0.9093                          0.7332
## ENSG00000162298                  0.7903                          0.5293
## ENSG00000129351                  0.8180                          0.7234
## ENSG00000099246                  1.1105                          1.1721
## ENSG00000132182                  0.8447                          0.5314
##                 BCR..mad HR..mad OxPhos..mad Not..BCR..mad Not..HR..mad
## ENSG00000121774   180.71  188.69      166.54        172.06       310.18
## ENSG00000099783   332.17  215.52      263.93        301.93       499.18
## ENSG00000162298   207.51  134.11      120.37        134.53       299.86
## ENSG00000129351   113.00   47.29       61.52         51.02       130.92
## ENSG00000099246    43.45   56.76       55.52         56.98        64.08
## ENSG00000132182    69.72   42.40       32.74         41.38        83.98
##                 Not..OxPhos..mad
## ENSG00000121774           230.72
## ENSG00000099783           313.44
## ENSG00000162298           276.36
## ENSG00000129351           112.60
## ENSG00000099246            68.06
## ENSG00000132182            76.17
```

The data.frame shown has a column for each class compared for each summary statistic (`mean`,`median`, `std.dev`, `mad`), and a mean `fold.change`, `t-score`, `p-value` and FDR-adjusted `p-value` for each class-class comparison done. The table is also sorted by FDR-adjusted `p-value`, smallest to largest. In this case, there were six classes (`BCR, HR, OxPhos, Not.BCR, Not.HR,` and `Not.OxPhos`) and three comparisons (`BCR-Not.BCR, HR-Not.HR` and `OxPhos-Not.OxPhos`). Each of these types of statistics is considered a "column.group", and can be accessed using the `column.groups` function to get the list of column groups in the table and the corresponding columns, or given the name of a group, say, "means", all columns from that group using the `group` function. 


```r
column.groups(t.test_meta.1)
```

```
## $p.values
## [1] "BCR..Not.BCR..p.value"       "HR..Not.HR..p.value"        
## [3] "OxPhos..Not.OxPhos..p.value"
## 
## $adj.p.values
## [1] "BCR..Not.BCR..adj.p.value"       "HR..Not.HR..adj.p.value"        
## [3] "OxPhos..Not.OxPhos..adj.p.value"
## 
## $t.scores
## [1] "BCR..Not.BCR..t.score"       "HR..Not.HR..t.score"        
## [3] "OxPhos..Not.OxPhos..t.score"
## 
## $means
## [1] "BCR..mean"        "HR..mean"         "OxPhos..mean"    
## [4] "Not.BCR..mean"    "Not.HR..mean"     "Not.OxPhos..mean"
## 
## $std.devs
## [1] "BCR..std.dev"        "HR..std.dev"         "OxPhos..std.dev"    
## [4] "Not.BCR..std.dev"    "Not.HR..std.dev"     "Not.OxPhos..std.dev"
## 
## $medians
## [1] "BCR..median"         "HR..median"          "OxPhos..median"     
## [4] "Not..BCR..median"    "Not..HR..median"     "Not..OxPhos..median"
## 
## $fold.changes
## [1] "BCR..Not.BCR..fold.change"       "HR..Not.HR..fold.change"        
## [3] "OxPhos..Not.OxPhos..fold.change"
## 
## $mads
## [1] "BCR..mad"         "HR..mad"          "OxPhos..mad"     
## [4] "Not..BCR..mad"    "Not..HR..mad"     "Not..OxPhos..mad"
## 
## $descriptives
## [1] "gene"
```

```r
head(group(t.test_meta.1, "means"))
```

```
##                 BCR..mean HR..mean OxPhos..mean Not.BCR..mean Not.HR..mean
## ENSG00000121774    1163.8    836.1        729.8         780.9        944.1
## ENSG00000099783    1985.4   1478.2       1274.7        1372.5       1625.6
## ENSG00000162298     678.8    377.5        281.7         327.7        477.7
## ENSG00000129351     451.0    294.4        271.1         282.3        359.9
## ENSG00000099246     292.3    381.5        393.5         387.7        343.5
## ENSG00000132182     250.9    150.5        107.4         128.1        178.2
##                 Not.OxPhos..mean
## ENSG00000121774           1004.3
## ENSG00000099783           1738.5
## ENSG00000162298            532.1
## ENSG00000129351            374.8
## ENSG00000099246            335.7
## ENSG00000132182            202.1
```

There are a lot of columns in this table, and usually, only one of the comparisons is interesting, the one that minimizes or maximizes the statistic. To assign each gene to a class, we can use the `unique.diffanal.results` (as an S3 method, we only need to write `unique`), specifying a column group and a scoring function to choose the best class for each gene from the set of scores in each group. The default column group is `adj.p.values` and the default scoring function is `base::which.min`. 


```r
ut.test_meta.1 <- unique(t.test_meta.1, criteria.column = "adj.p.values", score.fn = which.min)
head(ut.test_meta.1)
```

```
##                            gene        class class..p.value
## ENSG00000121774 ENSG00000121774 BCR..Not.BCR      4.814e-18
## ENSG00000099783 ENSG00000099783 BCR..Not.BCR      3.452e-17
## ENSG00000162298 ENSG00000162298 BCR..Not.BCR      3.411e-17
## ENSG00000129351 ENSG00000129351 BCR..Not.BCR      1.140e-16
## ENSG00000099246 ENSG00000099246 BCR..Not.BCR      2.190e-16
## ENSG00000132182 ENSG00000132182 BCR..Not.BCR      2.783e-16
##                 class..adj.p.value class..t.score class..mean
## ENSG00000121774          2.407e-14         -10.40      1163.8
## ENSG00000099783          5.754e-14         -10.22      1985.4
## ENSG00000162298          5.754e-14         -10.12       678.8
## ENSG00000129351          1.424e-13         -10.44       451.0
## ENSG00000099246          2.190e-13          10.41       292.3
## ENSG00000132182          2.319e-13         -10.21       250.9
##                 class..std.dev class..median class..fold.change class..mad
## ENSG00000121774         151.07        1149.6             1.4904     180.71
## ENSG00000099783         336.49        1890.6             1.4466     332.17
## ENSG00000162298         211.46         664.8             2.0711     207.51
## ENSG00000129351          98.66         452.6             1.5976     113.00
## ENSG00000099246          39.83         288.4             0.7539      43.45
## ENSG00000132182          83.05         227.5             1.9585      69.72
```

All of the statistic columns have been reduced to just `class..<statistic>`, and a new `class` column was added, showing which class that gene was assigned to. Note that since the gene order hasn't changed, many of the genes will appear to be sorted by class, but this is a side effect of the initial ordering, not of the `unique` function.

## Using `limma` with diffanal2

The t-test strategy is simple, and that has advantages and disadvantages. We can also fit a gene-wise linear model using the `limma` package within `diffanal2`. The call is very similar, but instead of using the t-test wrapper, `df2.t.test`, we use the limma wrapper, `df2.limma`. This function takes most of the same arguments that the t-test version does, but it also handles confounders (more on that later.). Consult ?


```r
limma_meta.1 <- df2.limma(lymphoma, model.formula= ~meta.1)
```

```
## Cleaning data, removing samples with NA phenotype data
## 0  samples dropped
## Variation filtering based on mad .. 
## done.
## 
## Selecting top 5000 by mad .. 
## done, 5000 genes selected.
## 
## Calculating summary statistics
```

```r
head(limma_meta.1)
```

```
##                            gene BCR..Not.BCR..p.value HR..Not.HR..p.value
## ENSG00000054118 ENSG00000054118             8.617e-09           0.0733121
## ENSG00000213020 ENSG00000213020             3.674e-08           0.0002891
## ENSG00000196584 ENSG00000196584             6.939e-07           0.0005985
## ENSG00000153310 ENSG00000153310             2.706e-06           0.0029840
## ENSG00000186517 ENSG00000186517             1.333e-05           0.2915297
## ENSG00000168802 ENSG00000168802             2.573e-05           0.1312186
##                 OxPhos..Not.OxPhos..p.value BCR..Not.BCR..adj.p.value
## ENSG00000054118                   0.0000241                 4.308e-05
## ENSG00000213020                   0.0359422                 9.184e-05
## ENSG00000196584                   0.0950319                 1.157e-03
## ENSG00000153310                   0.0638451                 3.383e-03
## ENSG00000186517                   0.0006691                 1.333e-02
## ENSG00000168802                   0.0050177                 1.689e-02
##                 HR..Not.HR..adj.p.value OxPhos..Not.OxPhos..adj.p.value
## ENSG00000054118                  0.6192                         0.02629
## ENSG00000213020                  0.1205                         0.33466
## ENSG00000196584                  0.1468                         0.44869
## ENSG00000153310                  0.2720                         0.39362
## ENSG00000186517                  0.7951                         0.09637
## ENSG00000168802                  0.6768                         0.19133
##                 BCR..Not.BCR..lods HR..Not.HR..lods
## ENSG00000054118              9.462          -4.2630
## ENSG00000213020              8.171           0.1916
## ENSG00000196584              5.558          -0.4153
## ENSG00000153310              4.350          -1.7416
## ENSG00000186517              2.939          -5.2079
## ENSG00000168802              2.359          -4.6819
##                 OxPhos..Not.OxPhos..lods BCR..Not.BCR..t.score
## ENSG00000054118                   2.3814                 6.203
## ENSG00000213020                  -3.7938                -5.896
## ENSG00000196584                  -4.5395                -5.248
## ENSG00000153310                  -4.2399                -4.933
## ENSG00000186517                  -0.4883                 4.547
## ENSG00000168802                  -2.1910                 4.382
##                 HR..Not.HR..t.score OxPhos..Not.OxPhos..t.score BCR..mean
## ENSG00000054118              -1.807                      -4.399     372.4
## ENSG00000213020               3.737                       2.122     114.9
## ENSG00000196584               3.528                       1.683     106.8
## ENSG00000153310               3.033                       1.871     101.8
## ENSG00000186517              -1.060                      -3.495     355.6
## ENSG00000168802              -1.520                      -2.860     483.5
##                 HR..mean OxPhos..mean Not.BCR..mean Not.HR..mean
## ENSG00000054118    294.1        272.9         283.3        318.8
## ENSG00000213020    167.9        157.2         162.5        134.4
## ENSG00000196584    152.1        140.9         146.4        122.7
## ENSG00000153310    165.0        153.3         159.0        125.0
## ENSG00000186517    284.9        259.1         271.7        303.6
## ENSG00000168802    420.6        407.9         414.2        444.1
##                 Not.OxPhos..mean BCR..std.dev HR..std.dev OxPhos..std.dev
## ENSG00000054118            331.0        1.117       1.121           1.116
## ENSG00000213020            138.9        1.117       1.121           1.116
## ENSG00000196584            127.4        1.117       1.121           1.116
## ENSG00000153310            129.6        1.117       1.121           1.116
## ENSG00000186517            318.3        1.117       1.121           1.116
## ENSG00000168802            451.0        1.117       1.121           1.116
##                 Not.BCR..std.dev Not.HR..std.dev Not.OxPhos..std.dev
## ENSG00000054118            1.118           1.117               1.119
## ENSG00000213020            1.118           1.117               1.119
## ENSG00000196584            1.118           1.117               1.119
## ENSG00000153310            1.118           1.117               1.119
## ENSG00000186517            1.118           1.117               1.119
## ENSG00000168802            1.118           1.117               1.119
##                 BCR..median HR..median OxPhos..median Not..BCR..median
## ENSG00000054118       342.2      276.0          307.8            291.4
## ENSG00000213020       151.0      160.6          141.0            145.7
## ENSG00000196584       136.0      143.6          120.9            132.4
## ENSG00000153310       166.3      149.0          117.4            133.9
## ENSG00000186517       318.3      320.3          249.1            279.6
## ENSG00000168802       440.8      409.4          457.3            434.9
##                 Not..HR..median Not..OxPhos..median
## ENSG00000054118           325.2               315.5
## ENSG00000213020           143.9               153.9
## ENSG00000196584           123.5               140.1
## ENSG00000153310           143.9               160.2
## ENSG00000186517           270.7               320.1
## ENSG00000168802           449.8               422.5
##                 BCR..Not.BCR..fold.change HR..Not.HR..fold.change
## ENSG00000054118                    1.3145                  0.9224
## ENSG00000213020                    0.7074                  1.2490
## ENSG00000196584                    0.7294                  1.2398
## ENSG00000153310                    0.6404                  1.3201
## ENSG00000186517                    1.3089                  0.9384
## ENSG00000168802                    1.1672                  0.9471
##                 OxPhos..Not.OxPhos..fold.change BCR..mad HR..mad
## ENSG00000054118                          0.8247    47.17   51.33
## ENSG00000213020                          1.1318    54.78   50.98
## ENSG00000196584                          1.1058    51.96   39.84
## ENSG00000153310                          1.1829    72.88   66.40
## ENSG00000186517                          0.8141    72.90   85.33
## ENSG00000168802                          0.9046    67.09   59.64
##                 OxPhos..mad Not..BCR..mad Not..HR..mad Not..OxPhos..mad
## ENSG00000054118       69.44         62.93        67.06            70.27
## ENSG00000213020       49.09         52.28        52.30            54.86
## ENSG00000196584       27.83         41.50        41.14            47.42
## ENSG00000153310       72.86         84.16        93.66            78.99
## ENSG00000186517       34.93         75.41        70.55            77.32
## ENSG00000168802       72.89         72.84        77.65            65.99
##                 F.score F.p.value
## ENSG00000054118  20.518 2.294e-08
## ENSG00000213020  17.727 1.881e-07
## ENSG00000196584  14.246 2.899e-06
## ENSG00000153310  12.337 1.373e-05
## ENSG00000186517  11.438 2.901e-05
## ENSG00000168802   9.962 1.011e-04
```

The result table for a limma strategy is similar to the t-test strategy, with a few key differences. There is an additional comparison statistic, `lods`, standing for log-odds of differential expression, which `limma::eBayes` computes. The class means and standard errors are also extracted from `limma::lmFit`, and the complementary class means are backsolved from `limma::contrast.fit`'s fold change coefficient. Additionally, the table has two additional un-classed statistics, `F.score` and `F.p.value`, which are not part of the other strategies by default (see `do.F.stat` parameter in `do.diffanal2` for how to add it).

The limma strategy can also calculate the influence of a confounding variable has on the mean expression of a gene, and use that information to calculate the "true" mean gene expression as a result of the predictor variable. A confounder can be specified by adding additional terms to the `model.formula` parameter, or by specifying a list of character strings to the `confounder` parameter. While there is no reason to believe so, what if we believed that the Class label in the phenotype data was a confounder? To test this, we would run:


```r
limma_meta.1xClass <- df2.limma(lymphoma, model.formula= ~ meta.1 + Class)
```

```
## Cleaning data, removing samples with NA phenotype data
## 0  samples dropped
## Variation filtering based on mad .. 
## done.
## 
## Selecting top 5000 by mad .. 
## done, 5000 genes selected.
## 
## Calculating summary statistics
```

```r
head(limma_meta.1xClass)
```

```
##                            gene BCR..Not.BCR..p.value HR..Not.HR..p.value
## ENSG00000180957 ENSG00000180957             3.353e-07           0.0061995
## ENSG00000184787 ENSG00000184787             1.953e-06           0.0008755
## ENSG00000213020 ENSG00000213020             2.215e-06           0.0019137
## ENSG00000173875 ENSG00000173875             3.936e-06           0.0002097
## ENSG00000120690 ENSG00000120690             1.142e-05           0.0770862
## ENSG00000088448 ENSG00000088448             1.782e-05           0.3224916
##                 OxPhos..Not.OxPhos..p.value BCR..Not.BCR..adj.p.value
## ENSG00000180957                   0.0044158                  0.001676
## ENSG00000184787                   0.0564685                  0.003692
## ENSG00000213020                   0.0365055                  0.003692
## ENSG00000173875                   0.1685358                  0.004919
## ENSG00000120690                   0.0034584                  0.008911
## ENSG00000088448                   0.0004827                  0.008911
##                 HR..Not.HR..adj.p.value OxPhos..Not.OxPhos..adj.p.value
## ENSG00000180957                 0.28330                          0.2066
## ENSG00000184787                 0.15614                          0.3567
## ENSG00000213020                 0.19528                          0.3131
## ENSG00000173875                 0.09367                          0.5046
## ENSG00000120690                 0.63705                          0.1973
## ENSG00000088448                 0.85360                          0.1508
##                 BCR..Not.BCR..lods HR..Not.HR..lods
## ENSG00000180957              6.279          -2.3407
## ENSG00000184787              4.693          -0.7247
## ENSG00000213020              4.580          -1.3752
## ENSG00000173875              4.064           0.4751
## ENSG00000120690              3.110          -4.3175
## ENSG00000088448              2.712          -5.2922
##                 OxPhos..Not.OxPhos..lods BCR..Not.BCR..t.score
## ENSG00000180957                  -2.0596                -5.418
## ENSG00000184787                  -4.0603                -5.014
## ENSG00000213020                  -3.7303                -4.984
## ENSG00000173875                  -4.8432                -4.848
## ENSG00000120690                  -1.8612                -4.589
## ENSG00000088448                  -0.2412                -4.478
##                 HR..Not.HR..t.score OxPhos..Not.OxPhos..t.score BCR..mean
## ENSG00000180957              2.7884                       2.904     462.8
## ENSG00000184787              3.4171                       1.927     275.9
## ENSG00000213020              3.1765                       2.116     118.1
## ENSG00000173875              3.8292                       1.386     217.6
## ENSG00000120690              1.7838                       2.986     684.8
## ENSG00000088448              0.9936                       3.593     387.6
##                 HR..mean OxPhos..mean Not.BCR..mean Not.HR..mean
## ENSG00000180957    780.4        799.5         789.9        608.3
## ENSG00000184787    452.7        421.4         436.8        341.0
## ENSG00000213020    187.7        179.4         183.5        145.5
## ENSG00000173875    470.7        387.8         427.2        290.5
## ENSG00000120690   1006.6       1095.2        1049.9        866.0
## ENSG00000088448    552.6        659.4         603.6        505.5
##                 Not.OxPhos..mean BCR..std.dev HR..std.dev OxPhos..std.dev
## ENSG00000180957            601.0        1.127       1.254           1.287
## ENSG00000184787            353.5        1.127       1.254           1.287
## ENSG00000213020            148.9        1.127       1.254           1.287
## ENSG00000173875            320.0        1.127       1.254           1.287
## ENSG00000120690            830.3        1.127       1.254           1.287
## ENSG00000088448            462.8        1.127       1.254           1.287
##                 Not.BCR..std.dev Not.HR..std.dev Not.OxPhos..std.dev
## ENSG00000180957             1.27           1.207                1.19
## ENSG00000184787             1.27           1.207                1.19
## ENSG00000213020             1.27           1.207                1.19
## ENSG00000173875             1.27           1.207                1.19
## ENSG00000120690             1.27           1.207                1.19
## ENSG00000088448             1.27           1.207                1.19
##                 BCR..median HR..median OxPhos..median Not..BCR..median
## ENSG00000180957       592.9      548.0          500.7            535.1
## ENSG00000184787       345.5      353.4          285.0            322.0
## ENSG00000213020       151.0      160.6          141.0            145.7
## ENSG00000173875       297.5      305.4          237.4            271.6
## ENSG00000120690       816.2      779.4          730.1            747.0
## ENSG00000088448       443.5      476.1          416.5            429.6
##                 Not..HR..median Not..OxPhos..median
## ENSG00000180957           551.7               584.3
## ENSG00000184787           317.6               349.5
## ENSG00000213020           143.9               153.9
## ENSG00000173875           264.7               301.2
## ENSG00000120690           747.0               789.3
## ENSG00000088448           426.4               459.8
##                 BCR..Not.BCR..fold.change HR..Not.HR..fold.change
## ENSG00000180957                    0.5859                   1.283
## ENSG00000184787                    0.6318                   1.328
## ENSG00000213020                    0.6433                   1.290
## ENSG00000173875                    0.5093                   1.620
## ENSG00000120690                    0.6522                   1.162
## ENSG00000088448                    0.6421                   1.093
##                 OxPhos..Not.OxPhos..fold.change BCR..mad HR..mad
## ENSG00000180957                           1.330    75.31   75.92
## ENSG00000184787                           1.192    59.51   44.07
## ENSG00000213020                           1.205    54.78   50.98
## ENSG00000173875                           1.212    89.41   57.81
## ENSG00000120690                           1.319   202.34  160.83
## ENSG00000088448                           1.425   141.81  174.30
##                 OxPhos..mad Not..BCR..mad Not..HR..mad Not..OxPhos..mad
## ENSG00000180957      118.82         96.16        92.94            77.57
## ENSG00000184787       83.76         82.82        80.21            56.01
## ENSG00000213020       49.09         52.28        52.30            54.86
## ENSG00000173875       94.66         82.27        90.77            78.14
## ENSG00000120690      183.71        159.45       181.75           163.11
## ENSG00000088448      156.16        177.42       155.91           153.97
##                 F.score F.p.value
## ENSG00000180957   14.74 2.001e-06
## ENSG00000184787   13.37 6.013e-06
## ENSG00000213020   12.93 8.591e-06
## ENSG00000173875   13.40 5.886e-06
## ENSG00000120690   10.59 5.987e-05
## ENSG00000088448   10.73 5.334e-05
```

To get an idea of how the addition of a confounder changes the mean calculation, try this plot:

```r
gs <- intersect(row.names(limma_meta.1), row.names(t.test_meta.1))
means.meta.1 <- group(limma_meta.1[gs,], "means")
means.meta.1xClass <- group(limma_meta.1xClass[gs,], "means")
means.frame <- cbind(m1 = means.meta.1, m2 = means.meta.1xClass)
for(i in seq(1:3)){
  plot(means.meta.1[,i] ~ means.meta.1xClass[,i], 
       ylab=paste("meta.1", colnames(means.meta.1)[i]), 
       xlab=paste("meta.1xClass", colnames(means.meta.1xClass)[i],
                  main = "Comparing difference in Means"))
  abline(0,1)
}
```

![plot of chunk change-in-mean](figure/change-in-mean1.png) ![plot of chunk change-in-mean](figure/change-in-mean2.png) ![plot of chunk change-in-mean](figure/change-in-mean3.png) 
The line through the origin denotes 

## Running permutation t-tests

A a core part of `diffanal` was the permutation test, which let it you permute the class/phenotype labels to test the robustness of the class results by using the `perm.t.test` strategy. This test also allows you to specify confounders, and will only generate label permutations that swap samples within the same confounder class. At this time, only binary confounders are supported. The number of permutations ran is controlled by the `nperm` parameter, and augmented by the `exhaustive` parameter, which limits the number of permutations to generate to the maximum number of unique permutations available. 

Because each class comparison is reperformed many, many times, this strategy can run each class comparison in parallel, controlling the number of parallel processors used by changing the `ncores` parameter. If `ncores` is 1, the default, the comparisons run sequentially. This multiprocessing is done using `plyr` and `doMC`. `doMC` makes initializing the environment of the child processes trivial because they are `fork`s of the parent process. 


```r
perm.t.meta.1 <- df2.perm.t.test(lymphoma, model.formula= ~ meta.1 + Class, 
                                 nperm = 100, ncores=3)
```

```
## Cleaning data, removing samples with NA phenotype data
## 0  samples dropped
## Variation filtering based on mad .. 
## done.
## 
## Selecting top 5000 by mad .. 
## done, 5000 genes selected.
## 
## Running tests in parallel
## Calculating summary statistics
## done
```

```r
head(perm.t.meta.1)
```

```
##                                      gene BCR..Not.BCR..asymp.p.value
## ENSG00000138755           ENSG00000138755                   1.738e-03
## ENSG00000204287           ENSG00000204287                   1.634e-03
## AFFX-HSAC07/X00351_3 AFFX-HSAC07/X00351_3                   7.258e-08
## ENSG00000034510           ENSG00000034510                   2.238e-03
## AFFX-HSAC07/X00351_M AFFX-HSAC07/X00351_M                   2.739e-09
## ENSG00000123416           ENSG00000123416                   8.677e-05
##                      HR..Not.HR..asymp.p.value
## ENSG00000138755                      3.176e-09
## ENSG00000204287                      9.230e-01
## AFFX-HSAC07/X00351_3                 1.900e-02
## ENSG00000034510                      1.688e-12
## AFFX-HSAC07/X00351_M                 4.551e-01
## ENSG00000123416                      1.242e-02
##                      OxPhos..Not.OxPhos..asymp.p.value
## ENSG00000138755                              6.569e-02
## ENSG00000204287                              1.497e-04
## AFFX-HSAC07/X00351_3                         2.962e-09
## ENSG00000034510                              1.216e-02
## AFFX-HSAC07/X00351_M                         6.008e-07
## ENSG00000123416                              1.612e-01
##                      BCR..Not.BCR..adj.asymp.p.value
## ENSG00000138755                            5.972e-03
## ENSG00000204287                            5.674e-03
## AFFX-HSAC07/X00351_3                       1.226e-06
## ENSG00000034510                            7.397e-03
## AFFX-HSAC07/X00351_M                       7.737e-08
## ENSG00000123416                            4.842e-04
##                      HR..Not.HR..adj.asymp.p.value
## ENSG00000138755                          4.609e-08
## ENSG00000204287                          9.378e-01
## AFFX-HSAC07/X00351_3                     3.246e-02
## ENSG00000034510                          1.592e-10
## AFFX-HSAC07/X00351_M                     5.237e-01
## ENSG00000123416                          2.234e-02
##                      OxPhos..Not.OxPhos..adj.asymp.p.value
## ENSG00000138755                                  1.038e-01
## ENSG00000204287                                  4.546e-04
## AFFX-HSAC07/X00351_3                             3.517e-08
## ENSG00000034510                                  2.337e-02
## AFFX-HSAC07/X00351_M                             3.469e-06
## ENSG00000123416                                  2.234e-01
##                      BCR..Not.BCR..asymp.t.score HR..Not.HR..asymp.t.score
## ENSG00000138755                            3.258                  -6.51963
## ENSG00000204287                            3.256                   0.09719
## AFFX-HSAC07/X00351_3                      -5.788                  -2.38033
## ENSG00000034510                            3.132                  -7.94107
## AFFX-HSAC07/X00351_M                      -6.472                  -0.75087
## ENSG00000123416                           -4.145                   2.57252
##                      OxPhos..Not.OxPhos..asymp.t.score
## ENSG00000138755                                  1.864
## ENSG00000204287                                 -3.923
## AFFX-HSAC07/X00351_3                             7.040
## ENSG00000034510                                  2.587
## AFFX-HSAC07/X00351_M                             5.560
## ENSG00000123416                                  1.412
##                      BCR..Not.BCR..perm.p.value HR..Not.HR..perm.p.value
## ENSG00000138755                        0.009901                 0.059406
## ENSG00000204287                        0.009901                 0.970297
## AFFX-HSAC07/X00351_3                   0.009901                 0.039604
## ENSG00000034510                        0.009901                 0.009901
## AFFX-HSAC07/X00351_M                   0.009901                 0.366337
## ENSG00000123416                        0.009901                 0.079208
##                      OxPhos..Not.OxPhos..perm.p.value
## ENSG00000138755                              0.950495
## ENSG00000204287                              0.009901
## AFFX-HSAC07/X00351_3                         0.009901
## ENSG00000034510                              0.504950
## AFFX-HSAC07/X00351_M                         0.009901
## ENSG00000123416                              0.049505
##                      BCR..Not.BCR..adj.perm.p.value
## ENSG00000138755                             0.06165
## ENSG00000204287                             0.06165
## AFFX-HSAC07/X00351_3                        0.06165
## ENSG00000034510                             0.06165
## AFFX-HSAC07/X00351_M                        0.06165
## ENSG00000123416                             0.06165
##                      HR..Not.HR..adj.perm.p.value
## ENSG00000138755                           0.12853
## ENSG00000204287                           0.97891
## AFFX-HSAC07/X00351_3                      0.09575
## ENSG00000034510                           0.03791
## AFFX-HSAC07/X00351_M                      0.48470
## ENSG00000123416                           0.15976
##                      OxPhos..Not.OxPhos..adj.perm.p.value
## ENSG00000138755                                   0.97487
## ENSG00000204287                                   0.03261
## AFFX-HSAC07/X00351_3                              0.03261
## ENSG00000034510                                   0.62587
## AFFX-HSAC07/X00351_M                              0.03261
## ENSG00000123416                                   0.10511
##                      BCR..Not.BCR..max.t.score HR..Not.HR..max.t.score
## ENSG00000138755                        3.09544                 -3.0772
## ENSG00000204287                        3.00984                  1.4849
## AFFX-HSAC07/X00351_3                  -1.36882                  0.3212
## ENSG00000034510                        2.92104                 -3.1923
## AFFX-HSAC07/X00351_M                  -1.53082                  1.6948
## ENSG00000123416                       -0.05415                  3.3254
##                      OxPhos..Not.OxPhos..max.t.score BCR..mean HR..mean
## ENSG00000138755                               4.4145      2684     6928
## ENSG00000204287                               0.6707      8986    10909
## AFFX-HSAC07/X00351_3                          5.1771     10850     9693
## ENSG00000034510                               4.4517      6923     9708
## AFFX-HSAC07/X00351_M                          4.4769      6477     4886
## ENSG00000123416                               2.2805     10233     8245
##                      OxPhos..mean Not.BCR..mean Not.HR..mean
## ENSG00000138755              3372          5081         3032
## ENSG00000204287             12702         11840        10867
## AFFX-HSAC07/X00351_3         7034          8312         8918
## ENSG00000034510              7188          8399         7057
## AFFX-HSAC07/X00351_M         3071          3943         4752
## ENSG00000123416              8530          8393         9371
##                      Not.OxPhos..mean BCR..std.dev HR..std.dev
## ENSG00000138755                  4750         2168        3205
## ENSG00000204287                  9922         3166        2657
## AFFX-HSAC07/X00351_3            10287         2240        1779
## ENSG00000034510                  8279         1403        1644
## AFFX-HSAC07/X00351_M             5702         1885        1654
## ENSG00000123416                  9265         2445        2400
##                      OxPhos..std.dev Not.BCR..std.dev Not.HR..std.dev
## ENSG00000138755                 3317             3702            2813
## ENSG00000204287                 3475             3218            3797
## AFFX-HSAC07/X00351_3            1944             2286            2832
## ENSG00000034510                 2439             2437            1988
## AFFX-HSAC07/X00351_M            1590             1851            2436
## ENSG00000123416                 1803             2101            2296
##                      Not.OxPhos..std.dev BCR..median HR..median
## ENSG00000138755                     3446        1806       6627
## ENSG00000204287                     3066        9005      11201
## AFFX-HSAC07/X00351_3                2098       10793       9338
## ENSG00000034510                     2064        7357       9443
## AFFX-HSAC07/X00351_M                1938        6557       4830
## ENSG00000123416                     2606       10002       7742
##                      OxPhos..median Not..BCR..median Not..HR..median
## ENSG00000138755                2080             4790            1904
## ENSG00000204287               11790            11589           11282
## AFFX-HSAC07/X00351_3           7128             8408            8856
## ENSG00000034510                6776             8802            7177
## AFFX-HSAC07/X00351_M           3178             4041            4391
## ENSG00000123416                8763             8159            9346
##                      Not..OxPhos..median BCR..Not.BCR..fold.change
## ENSG00000138755                     4399                    0.5281
## ENSG00000204287                    10706                    0.7589
## AFFX-HSAC07/X00351_3                9844                    1.3054
## ENSG00000034510                     8012                    0.8243
## AFFX-HSAC07/X00351_M                5561                    1.6427
## ENSG00000123416                     9343                    1.2192
##                      HR..Not.HR..fold.change
## ENSG00000138755                       2.2848
## ENSG00000204287                       1.0038
## AFFX-HSAC07/X00351_3                  1.0869
## ENSG00000034510                       1.3757
## AFFX-HSAC07/X00351_M                  1.0281
## ENSG00000123416                       0.8799
##                      OxPhos..Not.OxPhos..fold.change BCR..mad HR..mad
## ENSG00000138755                               0.7100     1588    2844
## ENSG00000204287                               1.2802     3527    2356
## AFFX-HSAC07/X00351_3                          0.6838     2199    1498
## ENSG00000034510                               0.8682     1390    1557
## AFFX-HSAC07/X00351_M                          0.5385     1861    1253
## ENSG00000123416                               0.9207     1504    2263
##                      OxPhos..mad Not..BCR..mad Not..HR..mad
## ENSG00000138755             1939          4293         1719
## ENSG00000204287             2257          2195         3214
## AFFX-HSAC07/X00351_3        1621          2076         2813
## ENSG00000034510             2448          2276         1828
## AFFX-HSAC07/X00351_M        1370          1764         2569
## ENSG00000123416             1566          1986         2003
##                      Not..OxPhos..mad
## ENSG00000138755                  4079
## ENSG00000204287                  2684
## AFFX-HSAC07/X00351_3             2128
## ENSG00000034510                  2031
## AFFX-HSAC07/X00351_M             1973
## ENSG00000123416                  2566
```

The permutation test strategy outputs the asymptotic t statistic, resulting p-value, and FDR p-value, as well as the permutation test's `max-t-score`, and the `permutation p-value` and the associated FDR. The results should be otherwise identical to the normal t-test strategy.


```r
gs2 <- intersect(row.names(perm.t.meta.1), row.names(t.test_meta.1))
cols <- colnames(group(perm.t.meta.1[gs2,], "means"))
for(i in seq(1:3)){
  plot(group(perm.t.meta.1[gs2,], "means")[,i], 
       group(t.test_meta.1[gs2,], "means")[,i], ylab= "t-test", 
       xlab= "permutation t-test",
       main = cols[i])
  abline(0,1)
}
```

![plot of chunk perm.t-vs-t](figure/perm.t-vs-t1.png) ![plot of chunk perm.t-vs-t](figure/perm.t-vs-t2.png) ![plot of chunk perm.t-vs-t](figure/perm.t-vs-t3.png) 

## Adding annotations to `diffanal.results`

To add annotations from biomaRt to a table, call the following function, (Note that it requires an internet connection)

```r
annotate.diffanal.results(table = limma_meta.1xClass,
                          dataset ="hsapiens_gene_ensembl",
                          mapping ="ensg",
                          symbol.idx ="hgnc_symbol",
                          na.rm = F)
```

```
## Request to BioMart web service failed. Verify if you are still connected to the internet.  Alternatively the BioMart web service is temporarily down.  Check http://www.biomart.org and verify if this website is available.
```

```
## Error: XML content does not seem to be XML:
```

```r
head(group(limma_meta.1xClass, "descriptives"))
```

```
##                            gene
## ENSG00000180957 ENSG00000180957
## ENSG00000184787 ENSG00000184787
## ENSG00000213020 ENSG00000213020
## ENSG00000173875 ENSG00000173875
## ENSG00000120690 ENSG00000120690
## ENSG00000088448 ENSG00000088448
```
The gene symbol and a short description are added to each row as part of the descriptives group. 
