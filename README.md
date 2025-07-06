# circMeta2

Circular RNA (circRNA) is one kind of single-stranded, covalently closed noncoding RNA, which exerts important biological functions by acting as transcriptional regulators, microRNA sponges and protein templates. One important analysis of circRNAs is performing differential expression analysis of circRNAs across different biological conditions. Considering the data characteristics of circRNAs such as low read counts, the DE analysis of circRNAs imposes a unique challenge. To address the challenge, we develop a novel computational pipeline named circMeta2 to perform differential expression analysis of alternative back-splicing circRNAs in clusters to improve the statistical power by leveraging the often-overlooked additive effects of individual circRNAs in the clusters identified by the alternative back-splicing events, which may help improve the statistical power and downstream biological findings.

# Installation

You can install circMeta2 from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fzhaouf/circMeta2")
```

# Reference and Case Group Definition

In circMeta2, experimental groups are defined by the order of samples in the `conditions` parameter:

- **Reference Group (condid=1)**: The first group of samples (control/baseline condition)
- **Case Group (condid=2)**: The second group of samples (treatment/disease condition)

**Example**: `conditions = c(2,2)` means:
- First 2 samples → Reference Group (condid=1) 
- Next 2 samples → Case Group (condid=2)

**Results Interpretation**:
- `m0`: Expected expression in reference group
- `m1`: Expected expression in case group  
- `log2fc`: Log2 fold change - positive values indicate up-regulation in case group

# Small sample size example

Two samples from human cerebellum (reference group) and two samples from human frontal cortex (case group) are used for demonstrating small sample size usage. circRNAs are called using CIRCexplorer2 and output files are stored in data folder.

``` r
library(circMeta2)

# read circexplorer2 or ciri2 output files using file path, specify number of sample for each 
# condition using conditions parameter and 2 samples for each condition.
circexplorers <- c(
  'path/to/cerebellum/SRR3192427.txt',
  'path/to/cerebellum/SRR3192428.txt',
  'path/to/frontalcortex/SRR3192424.txt',
  'path/to/frontalcortex/SRR3192425.txt'
)

circ.obj <- makecircObj(
  samplefiles = circexplorers,
  conditions = c(2, 2),
  circ.method = 'CIRCexplorer2',
  cutoff = 2
)
#> Process /path/to/cerebellum/SRR3192427.txt
#> Process /path/to/cerebellum/SRR3192428.txt
#> Process /path/to/frontalcortex/SRR3192424.txt
#> Process /path/to/frontalcortex/SRR3192425.txt

#  clustering A5BS and A3BS circ-clusters
circ.obj <- getCircCluster(circObj = circ.obj)

# individual circRNA DE using either pois Z for small sample or GLM for large sample with covariates
circ.obj <- circRNADE(circObj = circ.obj, DEmethod = 'Pois')

# circ-cluster DE
circ.obj <- circClusterDE(circObj = circ.obj, circ.cutoff = 2, DEmethod = 'Meta')
#> [==========================================] 100% 00:01:13
#> [==========================================] 100% 00:01:15

# Get results with simple accessor functions
A5BS_DE <- getClusterDE(circ.obj, cluster_type = "A5BS")
head(A5BS_DE)
#>    juncid numcircs        m0        m1         fc      pvalue       fdr cluster_type
#>        6        3 2.4083128 5.1912220  2.1555431 0.030929968 0.3482280         A5BS
#>        7        2 0.0200000 0.7863781 39.3189025 0.224627447 0.5778582         A5BS
#>       13        2 0.0200000 0.9471272 47.3563596 0.173450086 0.5673436         A5BS
#>       14        3 0.7191564 4.2020724  5.8430576 0.021024814 0.2909171         A5BS
#>       15        3 0.0300000 1.8649071 62.1635705 0.055562415 0.4454473         A5BS
#>       17        2 8.6886731 2.2984600  0.2645352 0.005971206 0.1388121         A5BS

A3BS_DE <- getClusterDE(circ.obj, cluster_type = "A3BS")
head(A3BS_DE)
#>    juncid numcircs       m0        m1         fc     pvalue       fdr cluster_type
#>        5        2 0.020000 0.7056718 35.2835879 0.25145167 0.5816162         A3BS
#>        6        2 1.520000 0.8063888  0.5305189 0.30567594 0.6110769         A3BS
#>        7        3 1.752891 1.2199202  0.6959476 0.01921387 0.2817073         A3BS
#>        9        2 0.020000 0.4228680 21.1434008 0.36938390 0.6126097         A3BS
#>       10        2 0.020000 1.3306479 66.5323969 0.10552125 0.5322011         A3BS
#>       11        2 1.020000 0.9071058  0.8893194 0.53652615 0.6126097         A3BS

# Get detailed information about specific clusters
getClusterAnnotation(circ.obj, cluster_type = "A5BS", juncid = c(6, 14, 17))
#>   juncid circid  chr     start       end width host_gene cluster_type
#>       6  11852 chr1 100335955 100379294 43340       AGL         A5BS
#>       6  20710 chr1 100335955 100343384  7430       AGL         A5BS
#>       6  37799 chr1 100335955 100347247 11293       AGL         A5BS
#>      14  20712 chr1 100459092 100483371 24280   SLC35A3         A5BS
#>      14  37803 chr1 100459092 100464971  5880   SLC35A3         A5BS
#>      14  37804 chr1 100459092 100480976 21885   SLC35A3         A5BS
#>      17  11853 chr1 100515464 100535241 19778     HIAT1         A5BS
#>      17  20713 chr1 100515464 100527505 12042     HIAT1         A5BS
```

# Large sample size example

AD circRNA data from GRanges with patient metadata including age and sex as covariates.

**Data Structure**: 
- **Reference group (condid=1)**: Control samples  
- **Case group (condid=2)**: AD patients
- **Covariates**: Age and sex to control for confounding factors

``` r
library(circMeta2)

data("BM10.circs", package = "circMeta2")
data("metainfo", package = "circMeta2")

# read circexplorer2 or ciri2 output files using file path, specify number of sample for each 
# condition using conditions parameter and the default is 2 samples for each condition.
circ.obj = makecircObjfromGRanges(GRanges = BM10.circs, metadata=metainfo)

#  clustering A5BS and A3BS circ-clusters
circ.obj = getCircCluster(circObj = circ.obj)

# call individual circRNA DE 
circ.obj = circRNADE(circObj = circ.obj, DEmethod = 'GLM', formula_str = "readNumber ~ condid + age + sex")
#> [==========================================] 100% 00:01:51

# circ-cluster DE
circ.obj = circClusterDE(circObj=circ.obj, circ.cutoff=2, DEmethod='Meta')
#> [==========================================] 100% 00:00:51
#> [==========================================] 100% 00:00:48

# Get results with simple accessor functions
A5BS_DE =  getClusterDE(circ.obj, cluster_type = "A5BS")
head(A5BS_DE)
#>    juncid numcircs         m0        m1        fc      pvalue       fdr cluster_type
#>      2        2 0.16117133 0.1321821 0.8201338 0.429049514 0.6532540         A5BS
#>      4        2 0.07259335 0.0887072 1.2219743 0.357621474 0.6018229         A5BS
#>     12        2 2.17063192 1.6688186 0.7688169 0.007636536 0.1063564         A5BS
A3BS_DE = getClusterDE(circ.obj, cluster_type = "A3BS")
head(A3BS_DE)
#>    juncid numcircs         m0         m1        fc      pvalue       fdr cluster_type
#>      10        2 14.5608643 29.5406839 2.0287727 0.008826275 0.1152943         A3BS
#>      14        2  0.5767828  1.0399432 1.8030068 0.087975056 0.3138943         A3BS
#>      18        2  0.3669290  0.8523791 2.3230081 0.106600377 0.3371806         A3BS

# Get detailed information about specific clusters
getClusterAnnotation(circ.obj, cluster_type = "A5BS", juncid = c(17))
  juncid circid  chr     start       end width host_gene cluster_type
#>     17  15157 chr1  10177516  10195244 17729     UBE4B         A5BS
#>     17  29014 chr1  10177516  10205097 27582     UBE4B         A5BS
#>     17  30156 chr1  10177516  10190871 13356     UBE4B         A5BS
#>     17  34214 chr1  10177516  10197263 19748     UBE4B         A5BS
```

# Understanding the Results

**Reference vs Case Groups:**
- `m0`: Expected expression in reference group (e.g., control samples)
- `m1`: Expected expression in case group (e.g., disease samples)  
- `log2fc`: Log2 fold change - positive values indicate upregulation in case group

**Individual vs Cluster Analysis:**
- Individual circRNA analysis tests each circRNA separately
- Cluster analysis combines related circRNAs sharing alternative splicing sites

# Main Functions

- `makecircObj()` / `makecircObjfromGRanges()`: Create analysis object
- `getCircCluster()`: Identify alternative splicing clusters  
- `circRNADE()`: Individual circRNA differential expression
- `circClusterDE()`: Cluster-level differential expression
- `getClusterAnnotation()`: detailed information about clusters
- `getCircRNADE()`, `getClusterDE()`: Access results

