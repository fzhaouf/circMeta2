
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- README.md is generated from README.Rmd. Please edit that file -->

# circMeta2

Circular RNA (circRNA) is one kind of single-stranded, covalently closed
noncoding RNA, which exerts important biological functions by acting as
transcriptional regulators, microRNA sponges and protein templates. One
important analysis of circRNAs is performing differential expression
analysis of circRNAs across different biological conditions. Considering
the data characteristics of circRNAs such as low read counts, the DE
analysis of circRNAs imposes a unique challenge. To address the
challenge, we develop a novel computational pipeline named circMeta2 to
perform differential expression analysis of alternative back-splicing
circRNAs in clusters to improve the statistical power by leveraging the
often-overlooked additive effects of individual circRNAs in the clusters
identified by the alternative back-splicing events, which may help
improve the statistical power and downstream biological findings.

# Installation

You can install circMeta2 from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fzhaouf/circMeta2")
```

# samll sample size example

two samples from human cerebellum and two samples from human frontal
cortex are used for demontrating small sample size usage. circRNAs are
called using ’CIRCexplorer2 and output files are stored in data folder.

``` r

library(circMeta2)
# read circexplorer2 or ciri2 output files using file path, specify number of sample for each 
# condition using conditions parameter and the default is 2 samples for each condition.
circexplorers <- c(
  'path/to/cerebellum/SRR3192427.txt',
  'path/to/cerebellum/SRR3192428.txt',
  'path/to/frontalcortex/SRR3192424.txt',
  'path/to/frontalcortex/SRR3192425.txt'
)

circ.obj = makecircObj(samplefiles = circexplorers, conditions = c(2,2), circ.method=c('CIRCexplorer2'))
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
#>     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
#>     table, tapply, union, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: GenomeInfoDb
#> Process ~/Library/CloudStorage/OneDrive-UniversityofFlorida/circRNA_project/Rpackage/circRNA/cerebellum/CIRCexplorer/SRR3192427.txt
#> Process ~/Library/CloudStorage/OneDrive-UniversityofFlorida/circRNA_project/Rpackage/circRNA/cerebellum/CIRCexplorer/SRR3192428.txt
#> Process ~/Library/CloudStorage/OneDrive-UniversityofFlorida/circRNA_project/Rpackage/circRNA/frontalcortex/CIRCexplorer/SRR3192424.txt
#> Process ~/Library/CloudStorage/OneDrive-UniversityofFlorida/circRNA_project/Rpackage/circRNA/frontalcortex/CIRCexplorer/SRR3192425.txt

#  clustering A5BS and A3BS circ-clusters
circ.obj = getCircCluster(circObj = circ.obj)
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:GenomicRanges':
#> 
#>     intersect, setdiff, union
#> The following object is masked from 'package:GenomeInfoDb':
#> 
#>     intersect
#> The following objects are masked from 'package:IRanges':
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following objects are masked from 'package:BiocGenerics':
#> 
#>     combine, intersect, setdiff, union
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

# individual circRNA DE using either pois Z for small sample or GLM for large sample with covariates
circ.obj = circRNADE(circObj = circ.obj, DEmethod = 'Pois')

# circ-cluster DE
results = circClusterDE(circObj=circ.obj, circ.cutoff=2, DEmethod='Meta')
head(results$A5BS.cluster)
#>    juncid numcircs        m0        m1      pvalue         fc       fdr
#> 6       6        3 2.4083128 5.1912220 0.030929968  2.1555431 0.3482280
#> 7       7        2 0.0200000 0.7863781 0.224627447 39.3189025 0.5778582
#> 13     13        2 0.0200000 0.9471272 0.173450086 47.3563596 0.5673436
#> 14     14        3 0.7191564 4.2020724 0.021024814  5.8430576 0.2909171
#> 15     15        3 0.0300000 1.8649071 0.055562415 62.1635705 0.4454473
#> 17     17        2 8.6886731 2.2984600 0.005971206  0.2645352 0.1388121
head(results$A3BS.cluster)
#>    juncid numcircs       m0        m1     pvalue         fc       fdr
#> 5       5        2 0.020000 0.7056718 0.25145167 35.2835879 0.5816162
#> 6       6        2 1.520000 0.8063888 0.30567594  0.5305189 0.6110769
#> 7       7        3 1.752891 1.2199202 0.01921387  0.6959476 0.2817073
#> 9       9        2 0.020000 0.4228680 0.36938390 21.1434008 0.6126097
#> 10     10        2 0.020000 1.3306479 0.10552125 66.5323969 0.5322011
#> 11     11        2 1.020000 0.9071058 0.53652615  0.8893194 0.6126097
```

# large sample size example

AD circRNA data from GRanges.

``` r
data("BM10.circs", package = "circMeta2")
data("metainfo", package = "circMeta2")

# read circexplorer2 or ciri2 output files using file path, specify number of sample for each 
# condition using conditions parameter and the default is 2 samples for each condtion.
circ.obj = makecircObjfromGRanges(GRanges = BM10.circs,metadata=metainfo)

#  clustering A5BS and A3BS circ-clusters
circ.obj = getCircCluster(circObj = circ.obj)

# call inidivual circRNA DE 
circ.obj = circRNADE(circObj = circ.obj, DEmethod = 'GLM', formula_str = "readNumber ~ condid + age + sex")

# circ-cluster DE
results = circClusterDE(circObj=circ.obj, circ.cutoff=2, DEmethod='Meta')
head(results$A5BS.cluster)
#>    juncid numcircs        m0       m1     pvalue        fc       fdr
#> 2       2        2  5.676516 5.680127 0.30571477 1.0006362 0.9839426
#> 4       4        2  5.152040 5.518598 0.32189768 1.0711481 0.9839426
#> 10     10        2  6.382077 6.338035 0.51963998 0.9930990 0.9839426
#> 12     12        2  6.192094 6.216405 0.41310838 1.0039261 0.9839426
#> 15     15        2 10.057018 9.542208 0.20357676 0.9488108 0.9839426
#> 16     16        2  6.550191 6.162054 0.03144616 0.9407442 0.9839426
head(results$A3BS.cluster)
#>    juncid numcircs       m0       m1    pvalue        fc       fdr
#> 4       4        2 7.584117 7.588470 0.5585926 1.0005740 0.9801751
#> 6       6        2 4.375104 4.594721 0.5294981 1.0501968 0.9801751
#> 10     10        2 5.281847 5.700707 0.2588774 1.0793019 0.9801751
#> 14     14        2 6.191477 6.643974 0.3300785 1.0730839 0.9801751
#> 18     18        2 5.803474 5.885271 0.5774544 1.0140943 0.9801751
#> 19     19        2 6.052055 5.638473 0.1428980 0.9316626 0.9801751
```
