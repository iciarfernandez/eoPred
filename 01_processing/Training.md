# Introduction

In this Rmd, I process the training group by following a pipeline of
sample and probe filtering. Data and metadata parsing code are not shown
in this script. I also BMIQ normalize and batch correct this group for
dataset effects. At each data transformation step, I run PCA to
visualize the data in a way that allows me to understand the effect of
the normalization and batch correction.

# Download data from GEO

For this project, I downloaded all of the datasets that I was
considering for use and assembled them into a
[MethylSet](https://rdrr.io/bioc/minfi/man/MethylSet-class.html), which
is comprised of *raw* data.

In the **Data** section, I subset to the four datasets that are included
in the training group, and apply the processing pipeline to those
samples in the subsequent sections. The four datasets are:

-   [GSE100197](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100197)

-   [GSE103253](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103253)

-   [GSE57767](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57767)

-   [GSE73375](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73375)

We thank all dataset authors for facilitating this work by making their
data public on GEO.

# Setup

In this section, I load the required packages and objects.

``` r
# PACKAGES

suppressPackageStartupMessages(library(tidyverse))
```

    ## Warning: package 'tidyverse' was built under R version 4.1.3

    ## Warning: package 'ggplot2' was built under R version 4.1.3

    ## Warning: package 'tibble' was built under R version 4.1.3

    ## Warning: package 'tidyr' was built under R version 4.1.3

    ## Warning: package 'readr' was built under R version 4.1.3

    ## Warning: package 'purrr' was built under R version 4.1.3

    ## Warning: package 'dplyr' was built under R version 4.1.3

    ## Warning: package 'stringr' was built under R version 4.1.3

    ## Warning: package 'forcats' was built under R version 4.1.3

``` r
suppressPackageStartupMessages(library(minfi))
```

    ## Warning: package 'S4Vectors' was built under R version 4.1.3

    ## Warning: package 'matrixStats' was built under R version 4.1.3

    ## Warning: package 'locfit' was built under R version 4.1.3

``` r
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(tictoc))
```

    ## Warning: package 'tictoc' was built under R version 4.1.3

``` r
suppressPackageStartupMessages(library(wateRmelon))
```

    ## Warning: package 'limma' was built under R version 4.1.3

    ## Warning: package 'scales' was built under R version 4.1.3

    ## No methods found in package 'RSQLite' for request: 'dbListFields' when loading 'lumi'

``` r
suppressPackageStartupMessages(library(ENmix))
```

    ## Warning: package 'ENmix' was built under R version 4.1.3

``` r
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(ewastools))
suppressPackageStartupMessages(library(plomics))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(impute))
suppressPackageStartupMessages(library(sva))
```

    ## Warning: package 'mgcv' was built under R version 4.1.3

    ## Warning: package 'nlme' was built under R version 4.1.3

``` r
suppressPackageStartupMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressPackageStartupMessages(library(IlluminaHumanMethylation450kmanifest))

# OBJECTS

mset <- readRDS(here("Data", "01_Raw", "complete_mset.rds")) # contains all eight datasets included in the study, directly downloaded and assembled as an mset from GEO by taking the methylated and unmethylated intensities

detp <- readRDS(here("Data", "01_Raw", "detp.rds")) # detection p-value (*not available for all datasets)

snps <- readRDS(here("Data", "01_Raw", "snps.rds")) # single nucleotide polymorphisms (*not available for all datasets)

bc <- readRDS(here("Data", "01_Raw", "beadcount.rds")) # beadcounts (*not available for all datasets)

# PROBE INFO

probeInfo <- cbind(IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations,
                   IlluminaHumanMethylation450kanno.ilmn12.hg19::Manifest,
                   IlluminaHumanMethylation450kanno.ilmn12.hg19::Other) %>% as.data.frame()

probeInfo <- rownames_to_column(probeInfo, var="probeID")

chrXY <- probeInfo %>% filter(chr %in% c("chrX", "chrY"))
chrX <- probeInfo %>% filter(chr %in% c("chrX"))
chrY <- probeInfo %>% filter(chr %in% c("chrY"))

zhouAnno <- readRDS("Z:/Amy/Data/Probe Info (Illumina Manifests)/Zhou/HM450.hg19.manifest.rds") %>% as.data.frame
```

# Data

I assemble the training group by selecting the datasets that will be
included in this group and subsetting the complete MethylSet, which
contains all the samples from the raw data of the eight datasets
included in this study.

``` r
pDat <- as.data.frame(pData(mset))

trainMeta <- pDat %>% 
  filter(GEO_Dataset %in% c("GSE100197", "GSE103253", "GSE57767", "GSE73375")) %>%
  filter(!Condition == "Replicate") %>% 
  filter(!Condition == "IUGR") %>%
  mutate(Class = case_when(Condition == "EOPE" ~ "EarlyPE",
                           Condition == "PE" & GA_RPC < 34 ~ "EarlyPE",
                           Condition == "Non-PE Preterm" ~ "Normotensive",
                           PE == "No" & GA_RPC < 37 ~ "Normotensive")) %>%
  mutate(Class = replace_na(Class, "Other")) %>%
  filter(!Class == "Other")

dim(trainMeta) # 89 by 23
```

    ## [1] 89 23

``` r
tabyl(trainMeta, Class) # 47 Early PE and 52 Normotensive
```

    ##         Class  n   percent
    ##       EarlyPE 47 0.5280899
    ##  Normotensive 42 0.4719101

``` r
train <- mset[,colnames(mset) %in% trainMeta$Sample_ID]
dim(train) # 485512 by 89
```

    ## [1] 485512     89

``` r
all(colnames(train) == trainMeta$Sample_ID) # TRUE
```

    ## [1] TRUE

# Sample QC

## Intersample correlation

``` r
train_beta <- getBeta(train)

all(colnames(train_beta) == trainMeta$Sample_ID) # TRUE
```

    ## [1] TRUE

``` r
cor <- cor(na.omit(train_beta))
sample_cor <- apply(cor, 1, mean)
sample_cor <- data.frame(mean_sample_cor = sample_cor, 
                         Sample_ID = as.factor(names(sample_cor)))
head(sample_cor)
```

    ##            mean_sample_cor  Sample_ID
    ## GSM1388241       0.9692632 GSM1388241
    ## GSM1388242       0.9714751 GSM1388242
    ## GSM1388243       0.9726865 GSM1388243
    ## GSM1388259       0.9723081 GSM1388259
    ## GSM1388262       0.9675114 GSM1388262
    ## GSM1388266       0.9634210 GSM1388266

``` r
all(sample_cor$Sample_ID == trainMeta$Sample_ID) # TRUE
```

    ## [1] TRUE

``` r
trainMeta$meanSScor <- sample_cor$mean_sample_cor

ggplot(trainMeta, aes(x = Sample_ID, y = meanSScor)) +
  geom_point(alpha=0.8) +
  geom_label(data = trainMeta %>% filter(meanSScor < 0.95), 
             aes(label = Sample_ID), vjust = 1.5, col = 'grey') +
  
  scale_color_manual(values=c("#117733", "#CC6677")) +
  scale_y_continuous(limits = c(0.9, 1), breaks = seq(0.9, 1, 0.01)) +
  
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust=0.5)) +
  labs(x = 'Sample', y = 'Mean correlation', col = 'Platform') +
  geom_hline(yintercept = 0.95, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = 0.96, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = 0.97, linetype = 'dashed', color = 'grey') +
  ggtitle("Mean Intersample Correlation - Normalized only")
```

![](00_Training_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
trainMeta %>%
  filter(meanSScor < 0.95)
```

    ##             Sample_ID Subject_ID GEO_Dataset    Sex Condition    IUGR Preterm
    ## GSM1892032 GSM1892032         Q4    GSE73375 Female        PE Unknown     Yes
    ## GSM1892059 GSM1892059        Q17    GSE73375 Female        PE Unknown     Yes
    ##            Term  PE GA   GA_RPC        Ethnicity Maternal_Age Chip  Row Column
    ## GSM1892032   No Yes 35 33.52367 African American           24 <NA> <NA>   <NA>
    ## GSM1892059   No Yes 30 29.63957 African American           32 <NA> <NA>   <NA>
    ##            Predicted_ethnicity_nothresh Predicted_ethnicity Prob_African
    ## GSM1892032                      African             African    0.9968682
    ## GSM1892059                      African             African    0.9761206
    ##             Prob_Asian Prob_Caucasian Highest_Prob   Class meanSScor
    ## GSM1892032 0.001783388    0.001348389    0.9968682 EarlyPE 0.9437336
    ## GSM1892059 0.014626659    0.009252736    0.9761206 EarlyPE 0.9411874

``` r
# flag samples
trainMeta <- 
  trainMeta %>%
  mutate(FAIL_Intersample_Cor = ifelse(meanSScor < 0.95, TRUE, FALSE))
```

## Sex

In this section, I check the match between genetic (XX/XY) sex and
reported sex. Amy Inkster wrote the following code to adapt the
`ewastools` check sex functions so that they can be used with an mset
rather than using the `read_idats` from `ewastools`.

``` r
# ALL CREDIT FOR CODE BELOW IS ATTRIBUTED TO EWASTOOLS PACKAGE AUTHORS (Heiss & Just)
# For ewastools package info please visit https://github.com/hhhh5/ewastools

# Original ewastools paper:
# Heiss, J., Just, A. Identifying mislabeled and contaminated DNA methylation
# microarray data: an extended quality control toolset with examples from GEO. 
# Clin Epigenet 10, 73 (2018). https://doi.org/10.1186/s13148-018-0504-1

# function equivalent to ewastools::check_sex(), pull sample XY intensity norm to auto

mset_check_sex = function(mset){
  
  if(!"tidyverse" %in% installed.packages()) {
    install.packages("tidyverse")
  }
  
  
  if(!"BiocManager" %in% installed.packages()) {
    install.packages("BiocManager")
  }
  
  if(!"minfi" %in% installed.packages()) {
    BiocManager::install("minfi")
  }
  
  require(tidyverse)
  require(minfi)
  
  # define M, U, and probe anno objects
  M = minfi::getMeth(mset) %>% as.data.frame()
  U = minfi::getUnmeth(mset) %>% as.data.frame()
  anno = minfi::getAnnotation(mset) %>% as.data.frame() # annotation with "chr" column, entries chr1 - chrY, probeIDs as rownames
    
  
  # select allosomal probes
  chrX = dimnames(anno[anno$chr=='chrX',])[[1]]
  chrY = dimnames(anno[anno$chr=='chrY',])[[1]]
  
  # compute the total intensities
  chrX = colMeans(M[chrX,,drop=FALSE]+U[chrX,,drop=FALSE],na.rm=TRUE)
  chrY = colMeans(M[chrY,,drop=FALSE]+U[chrY,,drop=FALSE],na.rm=TRUE)
  
  # compute the average total intensity across all autosomal probes
  autosomes = dimnames(anno[!anno$chr %in% c("chrX","chrY"),])[[1]]
  autosomes = colMeans(M[autosomes,,drop=FALSE]+U[autosomes,,drop=FALSE],na.rm=TRUE)
  
  # normalize total intensities
  chrX_nonorm = chrX
  chrY_nonorm = chrY
  
  chrX = chrX/autosomes
  chrY = chrY/autosomes
  
  # define a df with the output
  xy_int = data.frame(Sample_ID = dimnames(mset)[[2]],
                      Array = annotation(mset)[[1]],
                      X = chrX,
                      Y = chrY,
                      X_nonorm = chrX_nonorm,
                      Y_nonorm = chrY_nonorm)
  
    return((xy_int))
}

mset_predict_sex = function(sex_int,male,female){

    # compute the robust Hodges-Lehmann estimator for the total intensity for X chr probes
    cutX = outer(sex_int$X[male],sex_int$X[female],"+")
    cutX = median(cutX)/2

    # ... likewise for Y chr probes
    cutY = outer(sex_int$Y[male],sex_int$Y[female],"+")
    cutY = median(cutY)/2

    # Prediction based on in which quadrant (cutX/cutY) samples fall
    DNAme_Sex = rep(NA,times=length(sex_int$X))
    DNAme_Sex[sex_int$X>=cutX & sex_int$Y<=cutY] =  "XX"
    DNAme_Sex[sex_int$X<=cutX & sex_int$Y>=cutY] =  "XY"
    factor(DNAme_Sex,levels=c("XY","XX"),labels=c("XY","XX"))
    
    sex_preds = cbind(sex_int, DNAme_Sex)
    
    results = list("Sex_Preds" = sex_preds, "cutX" = cutX, "cutY" = cutY) 
    
    return(results)
}
```

``` r
pred_sex <- mset_check_sex(train)

male <- trainMeta$Sex == "Male"
female <- trainMeta$Sex == "Female"

pred_sex <- mset_predict_sex(pred_sex, male, female)
pred_sex_df <- pred_sex$Sex_Preds

trainMeta <- trainMeta %>% left_join(pred_sex_df %>% dplyr::select(Sample_ID, X, Y, DNAme_Sex), by="Sample_ID")
table(trainMeta$Sex, trainMeta$DNAme_Sex) # two disagreements
```

    ##         
    ##          XX XY
    ##   Female 38  2
    ##   Male    0 48

``` r
trainMeta <- trainMeta %>% mutate(FAIL_DNAme_Sex =
                                    ifelse((Sex == "Male" & DNAme_Sex == "XY") |
                                   (Sex == "Female" & DNAme_Sex == "XX"), 
                                 FALSE, TRUE))
table(trainMeta$FAIL_DNAme_Sex) 
```

    ## 
    ## FALSE  TRUE 
    ##    86     2

``` r
trainMeta %>% filter(FAIL_DNAme_Sex == T) 
```

    ##    Sample_ID Subject_ID GEO_Dataset    Sex      Condition    IUGR Preterm Term
    ## 1 GSM1892055        Q13    GSE73375 Female             PE Unknown     Yes   No
    ## 2 GSM2759051       <NA>   GSE103253 Female Non-PE Preterm Unknown     Yes   No
    ##    PE GA   GA_RPC Ethnicity Maternal_Age       Chip  Row Column
    ## 1 Yes 22 27.39576     Other           24       <NA> <NA>   <NA>
    ## 2  No NA 32.31346      <NA>           NA 9344730069  R06    C02
    ##   Predicted_ethnicity_nothresh Predicted_ethnicity Prob_African  Prob_Asian
    ## 1                      African             African   0.99447380 0.001762925
    ## 2                    Caucasian           Caucasian   0.02442828 0.038016082
    ##   Prob_Caucasian Highest_Prob        Class meanSScor FAIL_Intersample_Cor
    ## 1    0.003763276    0.9944738      EarlyPE 0.9553686                FALSE
    ## 2    0.937555640    0.9375556 Normotensive 0.9656976                FALSE
    ##           X         Y DNAme_Sex FAIL_DNAme_Sex
    ## 1 0.8955116 0.7280090        XY           TRUE
    ## 2 0.8683297 0.7986473        XY           TRUE

``` r
pred_sex_df %>%
  
  ggplot(aes(x=X, y=Y, color=DNAme_Sex)) +
  
  geom_point() +
  scale_color_manual(values=c("#997700", "#6699CC") ) +
  geom_vline(xintercept = pred_sex$cutX) +
  geom_hline(yintercept = pred_sex$cutY) +
  
  labs(title="DNAme-Based Sex",
       color="Sex",
       x= "chrX Normalized Total Fluo Intensity",
       y="chrY Normalized Total Fluo Intensity") +
  
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5))
```

![](00_Training_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
# remove sex-mismatched samples
dim(trainMeta) # 89
```

    ## [1] 89 29

``` r
trainMeta <- trainMeta %>% filter(!Sample_ID %in% c("GSM1892055", "GSM2759051"))
dim(trainMeta) # 87
```

    ## [1] 87 29

``` r
train <- train[,colnames(train) %in% trainMeta$Sample_ID]

dim(train) # 87
```

    ## [1] 485512     87

## SNP identity with `ewastools`

``` r
# ALL CREDIT FOR CODE BELOW IS ATTRIBUTED TO EWASTOOLS PACKAGE AUTHORS (Heiss & Just)
# For ewastools package info please visit https://github.com/hhhh5/ewastools

# Original ewastools paper:
# Heiss, J., Just, A. Identifying mislabeled and contaminated DNA methylation
# microarray data: an extended quality control toolset with examples from GEO. 
# Clin Epigenet 10, 73 (2018). https://doi.org/10.1186/s13148-018-0504-1

# function equivalent to ewastools::check_sex(), pull sample XY intensity norm to auto

# note that not all samples had SNP information

dim(snps) # 297
```

    ## [1]  65 419

``` r
trainSNP <- snps[,colnames(snps) %in% trainMeta$Sample_ID]
dim(trainSNP) # 65 by 68
```

    ## [1] 65 87

``` r
trainGenotypes <- call_genotypes(trainSNP)

# label samples with donor ID and number of times present in data genetically

donors <- colnames(trainSNP)
samples <- colnames(trainSNP)

check_snp_agreement(trainGenotypes, donor_ids = donors, sample_ids = samples) # 0
```

    ## NULL

## SNP contamination with `ewastools`

``` r
dim(trainMeta) # 87
```

    ## [1] 87 29

``` r
trainMeta_snps <- trainMeta %>% filter(Sample_ID %in% colnames(trainSNP))
dim(trainMeta_snps) # 68
```

    ## [1] 87 29

``` r
trainMeta_snps <- trainMeta_snps %>%
  dplyr::mutate(SNP_Outlier = snp_outliers(trainGenotypes))
 
ggplot(trainMeta_snps, aes(x = Sample_ID, y = SNP_Outlier, col = GEO_Dataset)) +
  geom_point(alpha=0.8) +
  geom_label(data = trainMeta_snps %>% filter(SNP_Outlier > -2.5), 
             aes(label = Sample_ID), vjust = 1.5, col = 'grey') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust=0.5)) +
  labs(x = 'Sample', y = 'Mean correlation', col = 'Dataset') +
  geom_hline(yintercept = -4, linetype = 'dashed', color = 'grey') +
  ggtitle("Mean SNP Outlier Posterior Probability")
```

    ## Warning: Removed 19 rows containing missing values (`geom_point()`).

![](00_Training_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
train_SNPoutliers <- trainMeta_snps %>% filter(SNP_Outlier > -4)

# i will wait to see if they also fail for failed probes and remove if 2x fail
trainMeta <- 
  trainMeta %>%
  mutate(FAIL_Contamination_Check = ifelse(Sample_ID %in% train_SNPoutliers$Sample_ID, TRUE, FALSE))
```

## Sample quality with `minfi`

``` r
qc <- getQC(train)
plotQC(qc) 
```

![](00_Training_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
train <- addQC(train, qc) # adds uMed and mMed rows to the phenodata of the mSet
as.data.frame(pData(train)) %>%
  filter(mMed < 10.5) 
```

    ##             Sample_ID Subject_ID GEO_Dataset  Sex Condition    IUGR Preterm
    ## GSM1892056 GSM1892056        Q14    GSE73375 Male        PE Unknown     Yes
    ##            Term  PE GA   GA_RPC        Ethnicity Maternal_Age Chip  Row Column
    ## GSM1892056   No Yes 27 30.09195 African American           34 <NA> <NA>   <NA>
    ##            Predicted_ethnicity_nothresh Predicted_ethnicity Prob_African
    ## GSM1892056                      African             African    0.9824547
    ##             Prob_Asian Prob_Caucasian Highest_Prob     mMed     uMed
    ## GSM1892056 0.001918254       0.015627    0.9824547 9.834471 10.91439

# Probe QC

``` r
## POOR QUALITY ##

# detection p value
train_detp <- detp[,colnames(detp) %in% trainMeta$Sample_ID]
dim(train_detp) 
```

    ## [1] 485512     87

``` r
dim(trainMeta) 
```

    ## [1] 87 30

``` r
# beadcount
train_bc <- bc[,colnames(bc) %in% trainMeta$Sample_ID]
dim(train_bc) # 50 - not all have beadcounts
```

    ## [1] 485512     50

``` r
head(train_bc)
```

    ##            GSM2674413 GSM2674414 GSM2674415 GSM2674416 GSM2674417 GSM2674418
    ## cg00050873         13         17          6         15         17         13
    ## cg00212031          9         16         15         10         21         18
    ## cg00213748         14          8          8         13         15         14
    ## cg00214611         18          9         18         19         14         12
    ## cg00455876         21         12         11         17         14         11
    ## cg01707559          9          9         10          7         14         11
    ##            GSM2674419 GSM2674420 GSM2674421 GSM2674422 GSM2674423 GSM2674424
    ## cg00050873         18         21         16         13          9         15
    ## cg00212031         13         11         15         13         14         19
    ## cg00213748         13          7         10         13          7         17
    ## cg00214611         16         22         14          5          9         15
    ## cg00455876         11         20         14         16         13          8
    ## cg01707559         10          5          7         11          5         16
    ##            GSM2674425 GSM2674426 GSM2674427 GSM2674428 GSM2674429 GSM2674430
    ## cg00050873         22          6         11          8          3         20
    ## cg00212031         13         10         12         17          5         14
    ## cg00213748         11          9          9         12          9          9
    ## cg00214611         10         10          9         12          5         22
    ## cg00455876         13         NA         10         15          7         11
    ## cg01707559         13         14         10          7          7         11
    ##            GSM2674431 GSM2674432 GSM2674433 GSM2674434 GSM2674464 GSM2674465
    ## cg00050873         17         21         16         17         20         18
    ## cg00212031         11         11         13         14         14          8
    ## cg00213748         21         15          9          8         13         10
    ## cg00214611         18         10         18         12         12         12
    ## cg00455876         11          8         13         20         15         16
    ## cg01707559          5         21          7          6         16          9
    ##            GSM2674466 GSM2674467 GSM2674468 GSM2674469 GSM2674470 GSM2674471
    ## cg00050873         16         10         13         22         16         18
    ## cg00212031         14         16         10         25         16          7
    ## cg00213748         12         13         14         14         12          9
    ## cg00214611         12         18         18         14         10         10
    ## cg00455876         16         12         12         16         13          6
    ## cg01707559         11         10         15          8         10          8
    ##            GSM2674472 GSM2674473 GSM2674474 GSM2674475 GSM2674476 GSM2674477
    ## cg00050873          9         16         20         11         14          5
    ## cg00212031         11         14         14         15         10          8
    ## cg00213748          5         12         13         11          6          6
    ## cg00214611          4         14         16         10         12          7
    ## cg00455876          8         20          9         11         16         11
    ## cg01707559          9         12          8         11         16         12
    ##            GSM2674478 GSM2674479 GSM2674480 GSM2674481 GSM2674482 GSM2674483
    ## cg00050873         16         21         19         22         18         13
    ## cg00212031         10         15         18         15          9         12
    ## cg00213748          6         13         11         15         11          9
    ## cg00214611          9         12         20         15         11         14
    ## cg00455876         15         11         13         10         20         13
    ## cg01707559         10          9         14         10         14         19
    ##            GSM2674484 GSM2674485 GSM2674486 GSM2674487 GSM2674499 GSM2674503
    ## cg00050873         11         20         14         16         16         18
    ## cg00212031          8         20          8         10         12         14
    ## cg00213748         10         15         14         13         11         13
    ## cg00214611          8         15         11         18         13          9
    ## cg00455876         12          9         15         14         10         19
    ## cg01707559         12         12         10         19         11         11
    ##            GSM2674511 GSM2674514
    ## cg00050873         22          5
    ## cg00212031         19          4
    ## cg00213748         14         11
    ## cg00214611         12          5
    ## cg00455876         11         13
    ## cg01707559         13         13

``` r
empty <- matrix(,nrow = 485512, ncol = 37)
head(empty)
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
    ## [1,]   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA
    ## [2,]   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA
    ## [3,]   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA
    ## [4,]   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA
    ## [5,]   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA
    ## [6,]   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA
    ##      [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
    ## [1,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## [2,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## [3,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## [4,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## [5,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## [6,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ##      [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37]
    ## [1,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## [2,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## [3,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## [4,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## [5,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## [6,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA

``` r
rownames(empty) <- rownames(train_bc)
head(empty)
```

    ##            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
    ## cg00050873   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA
    ## cg00212031   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA
    ## cg00213748   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA
    ## cg00214611   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA
    ## cg00455876   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA
    ## cg01707559   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA
    ##            [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
    ## cg00050873    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## cg00212031    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## cg00213748    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## cg00214611    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## cg00455876    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## cg01707559    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ##            [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35]
    ## cg00050873    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## cg00212031    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## cg00213748    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## cg00214611    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## cg00455876    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ## cg01707559    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
    ##            [,36] [,37]
    ## cg00050873    NA    NA
    ## cg00212031    NA    NA
    ## cg00213748    NA    NA
    ## cg00214611    NA    NA
    ## cg00455876    NA    NA
    ## cg01707559    NA    NA

``` r
samples_with_detp <- colnames(train_detp) # all samples
samples_with_bc <- colnames(train_bc) # only 50 out of the 87

length(samples_with_detp_without_bc <- setdiff(samples_with_detp, samples_with_bc)) # 37, makes sense
```

    ## [1] 37

``` r
colnames(empty) <- samples_with_detp_without_bc
head(empty)
```

    ##            GSM1388241 GSM1388242 GSM1388243 GSM1388259 GSM1388262 GSM1388266
    ## cg00050873         NA         NA         NA         NA         NA         NA
    ## cg00212031         NA         NA         NA         NA         NA         NA
    ## cg00213748         NA         NA         NA         NA         NA         NA
    ## cg00214611         NA         NA         NA         NA         NA         NA
    ## cg00455876         NA         NA         NA         NA         NA         NA
    ## cg01707559         NA         NA         NA         NA         NA         NA
    ##            GSM1388269 GSM1388271 GSM1892029 GSM1892032 GSM1892034 GSM1892044
    ## cg00050873         NA         NA         NA         NA         NA         NA
    ## cg00212031         NA         NA         NA         NA         NA         NA
    ## cg00213748         NA         NA         NA         NA         NA         NA
    ## cg00214611         NA         NA         NA         NA         NA         NA
    ## cg00455876         NA         NA         NA         NA         NA         NA
    ## cg01707559         NA         NA         NA         NA         NA         NA
    ##            GSM1892046 GSM1892047 GSM1892048 GSM1892056 GSM1892057 GSM1892059
    ## cg00050873         NA         NA         NA         NA         NA         NA
    ## cg00212031         NA         NA         NA         NA         NA         NA
    ## cg00213748         NA         NA         NA         NA         NA         NA
    ## cg00214611         NA         NA         NA         NA         NA         NA
    ## cg00455876         NA         NA         NA         NA         NA         NA
    ## cg01707559         NA         NA         NA         NA         NA         NA
    ##            GSM2759050 GSM2759052 GSM2759053 GSM2759054 GSM2759056 GSM2759057
    ## cg00050873         NA         NA         NA         NA         NA         NA
    ## cg00212031         NA         NA         NA         NA         NA         NA
    ## cg00213748         NA         NA         NA         NA         NA         NA
    ## cg00214611         NA         NA         NA         NA         NA         NA
    ## cg00455876         NA         NA         NA         NA         NA         NA
    ## cg01707559         NA         NA         NA         NA         NA         NA
    ##            GSM2759061 GSM2759062 GSM2759086 GSM2759089 GSM2759090 GSM2759091
    ## cg00050873         NA         NA         NA         NA         NA         NA
    ## cg00212031         NA         NA         NA         NA         NA         NA
    ## cg00213748         NA         NA         NA         NA         NA         NA
    ## cg00214611         NA         NA         NA         NA         NA         NA
    ## cg00455876         NA         NA         NA         NA         NA         NA
    ## cg01707559         NA         NA         NA         NA         NA         NA
    ##            GSM2759092 GSM2759093 GSM2759094 GSM2759095 GSM2759096 GSM2759098
    ## cg00050873         NA         NA         NA         NA         NA         NA
    ## cg00212031         NA         NA         NA         NA         NA         NA
    ## cg00213748         NA         NA         NA         NA         NA         NA
    ## cg00214611         NA         NA         NA         NA         NA         NA
    ## cg00455876         NA         NA         NA         NA         NA         NA
    ## cg01707559         NA         NA         NA         NA         NA         NA
    ##            GSM2759101
    ## cg00050873         NA
    ## cg00212031         NA
    ## cg00213748         NA
    ## cg00214611         NA
    ## cg00455876         NA
    ## cg01707559         NA

``` r
# substitute NAs of the samples i do not have beadcounts for for any number >3, or else they will be dropped for no reason
empty[is.na(empty)] <- 10
head(empty)
```

    ##            GSM1388241 GSM1388242 GSM1388243 GSM1388259 GSM1388262 GSM1388266
    ## cg00050873         10         10         10         10         10         10
    ## cg00212031         10         10         10         10         10         10
    ## cg00213748         10         10         10         10         10         10
    ## cg00214611         10         10         10         10         10         10
    ## cg00455876         10         10         10         10         10         10
    ## cg01707559         10         10         10         10         10         10
    ##            GSM1388269 GSM1388271 GSM1892029 GSM1892032 GSM1892034 GSM1892044
    ## cg00050873         10         10         10         10         10         10
    ## cg00212031         10         10         10         10         10         10
    ## cg00213748         10         10         10         10         10         10
    ## cg00214611         10         10         10         10         10         10
    ## cg00455876         10         10         10         10         10         10
    ## cg01707559         10         10         10         10         10         10
    ##            GSM1892046 GSM1892047 GSM1892048 GSM1892056 GSM1892057 GSM1892059
    ## cg00050873         10         10         10         10         10         10
    ## cg00212031         10         10         10         10         10         10
    ## cg00213748         10         10         10         10         10         10
    ## cg00214611         10         10         10         10         10         10
    ## cg00455876         10         10         10         10         10         10
    ## cg01707559         10         10         10         10         10         10
    ##            GSM2759050 GSM2759052 GSM2759053 GSM2759054 GSM2759056 GSM2759057
    ## cg00050873         10         10         10         10         10         10
    ## cg00212031         10         10         10         10         10         10
    ## cg00213748         10         10         10         10         10         10
    ## cg00214611         10         10         10         10         10         10
    ## cg00455876         10         10         10         10         10         10
    ## cg01707559         10         10         10         10         10         10
    ##            GSM2759061 GSM2759062 GSM2759086 GSM2759089 GSM2759090 GSM2759091
    ## cg00050873         10         10         10         10         10         10
    ## cg00212031         10         10         10         10         10         10
    ## cg00213748         10         10         10         10         10         10
    ## cg00214611         10         10         10         10         10         10
    ## cg00455876         10         10         10         10         10         10
    ## cg01707559         10         10         10         10         10         10
    ##            GSM2759092 GSM2759093 GSM2759094 GSM2759095 GSM2759096 GSM2759098
    ## cg00050873         10         10         10         10         10         10
    ## cg00212031         10         10         10         10         10         10
    ## cg00213748         10         10         10         10         10         10
    ## cg00214611         10         10         10         10         10         10
    ## cg00455876         10         10         10         10         10         10
    ## cg01707559         10         10         10         10         10         10
    ##            GSM2759101
    ## cg00050873         10
    ## cg00212031         10
    ## cg00213748         10
    ## cg00214611         10
    ## cg00455876         10
    ## cg01707559         10

``` r
train_bc <- cbind(train_bc, empty)
head(train_bc)
```

    ##            GSM2674413 GSM2674414 GSM2674415 GSM2674416 GSM2674417 GSM2674418
    ## cg00050873         13         17          6         15         17         13
    ## cg00212031          9         16         15         10         21         18
    ## cg00213748         14          8          8         13         15         14
    ## cg00214611         18          9         18         19         14         12
    ## cg00455876         21         12         11         17         14         11
    ## cg01707559          9          9         10          7         14         11
    ##            GSM2674419 GSM2674420 GSM2674421 GSM2674422 GSM2674423 GSM2674424
    ## cg00050873         18         21         16         13          9         15
    ## cg00212031         13         11         15         13         14         19
    ## cg00213748         13          7         10         13          7         17
    ## cg00214611         16         22         14          5          9         15
    ## cg00455876         11         20         14         16         13          8
    ## cg01707559         10          5          7         11          5         16
    ##            GSM2674425 GSM2674426 GSM2674427 GSM2674428 GSM2674429 GSM2674430
    ## cg00050873         22          6         11          8          3         20
    ## cg00212031         13         10         12         17          5         14
    ## cg00213748         11          9          9         12          9          9
    ## cg00214611         10         10          9         12          5         22
    ## cg00455876         13         NA         10         15          7         11
    ## cg01707559         13         14         10          7          7         11
    ##            GSM2674431 GSM2674432 GSM2674433 GSM2674434 GSM2674464 GSM2674465
    ## cg00050873         17         21         16         17         20         18
    ## cg00212031         11         11         13         14         14          8
    ## cg00213748         21         15          9          8         13         10
    ## cg00214611         18         10         18         12         12         12
    ## cg00455876         11          8         13         20         15         16
    ## cg01707559          5         21          7          6         16          9
    ##            GSM2674466 GSM2674467 GSM2674468 GSM2674469 GSM2674470 GSM2674471
    ## cg00050873         16         10         13         22         16         18
    ## cg00212031         14         16         10         25         16          7
    ## cg00213748         12         13         14         14         12          9
    ## cg00214611         12         18         18         14         10         10
    ## cg00455876         16         12         12         16         13          6
    ## cg01707559         11         10         15          8         10          8
    ##            GSM2674472 GSM2674473 GSM2674474 GSM2674475 GSM2674476 GSM2674477
    ## cg00050873          9         16         20         11         14          5
    ## cg00212031         11         14         14         15         10          8
    ## cg00213748          5         12         13         11          6          6
    ## cg00214611          4         14         16         10         12          7
    ## cg00455876          8         20          9         11         16         11
    ## cg01707559          9         12          8         11         16         12
    ##            GSM2674478 GSM2674479 GSM2674480 GSM2674481 GSM2674482 GSM2674483
    ## cg00050873         16         21         19         22         18         13
    ## cg00212031         10         15         18         15          9         12
    ## cg00213748          6         13         11         15         11          9
    ## cg00214611          9         12         20         15         11         14
    ## cg00455876         15         11         13         10         20         13
    ## cg01707559         10          9         14         10         14         19
    ##            GSM2674484 GSM2674485 GSM2674486 GSM2674487 GSM2674499 GSM2674503
    ## cg00050873         11         20         14         16         16         18
    ## cg00212031          8         20          8         10         12         14
    ## cg00213748         10         15         14         13         11         13
    ## cg00214611          8         15         11         18         13          9
    ## cg00455876         12          9         15         14         10         19
    ## cg01707559         12         12         10         19         11         11
    ##            GSM2674511 GSM2674514 GSM1388241 GSM1388242 GSM1388243 GSM1388259
    ## cg00050873         22          5         10         10         10         10
    ## cg00212031         19          4         10         10         10         10
    ## cg00213748         14         11         10         10         10         10
    ## cg00214611         12          5         10         10         10         10
    ## cg00455876         11         13         10         10         10         10
    ## cg01707559         13         13         10         10         10         10
    ##            GSM1388262 GSM1388266 GSM1388269 GSM1388271 GSM1892029 GSM1892032
    ## cg00050873         10         10         10         10         10         10
    ## cg00212031         10         10         10         10         10         10
    ## cg00213748         10         10         10         10         10         10
    ## cg00214611         10         10         10         10         10         10
    ## cg00455876         10         10         10         10         10         10
    ## cg01707559         10         10         10         10         10         10
    ##            GSM1892034 GSM1892044 GSM1892046 GSM1892047 GSM1892048 GSM1892056
    ## cg00050873         10         10         10         10         10         10
    ## cg00212031         10         10         10         10         10         10
    ## cg00213748         10         10         10         10         10         10
    ## cg00214611         10         10         10         10         10         10
    ## cg00455876         10         10         10         10         10         10
    ## cg01707559         10         10         10         10         10         10
    ##            GSM1892057 GSM1892059 GSM2759050 GSM2759052 GSM2759053 GSM2759054
    ## cg00050873         10         10         10         10         10         10
    ## cg00212031         10         10         10         10         10         10
    ## cg00213748         10         10         10         10         10         10
    ## cg00214611         10         10         10         10         10         10
    ## cg00455876         10         10         10         10         10         10
    ## cg01707559         10         10         10         10         10         10
    ##            GSM2759056 GSM2759057 GSM2759061 GSM2759062 GSM2759086 GSM2759089
    ## cg00050873         10         10         10         10         10         10
    ## cg00212031         10         10         10         10         10         10
    ## cg00213748         10         10         10         10         10         10
    ## cg00214611         10         10         10         10         10         10
    ## cg00455876         10         10         10         10         10         10
    ## cg01707559         10         10         10         10         10         10
    ##            GSM2759090 GSM2759091 GSM2759092 GSM2759093 GSM2759094 GSM2759095
    ## cg00050873         10         10         10         10         10         10
    ## cg00212031         10         10         10         10         10         10
    ## cg00213748         10         10         10         10         10         10
    ## cg00214611         10         10         10         10         10         10
    ## cg00455876         10         10         10         10         10         10
    ## cg01707559         10         10         10         10         10         10
    ##            GSM2759096 GSM2759098 GSM2759101
    ## cg00050873         10         10         10
    ## cg00212031         10         10         10
    ## cg00213748         10         10         10
    ## cg00214611         10         10         10
    ## cg00455876         10         10         10
    ## cg01707559         10         10         10

``` r
all(colnames(train_detp) == trainMeta$Sample_ID) 
```

    ## [1] TRUE

``` r
all(colnames(train_bc) == trainMeta$Sample_ID) 
```

    ## [1] FALSE

``` r
# reorder
train_bc <- train_bc[,colnames(train_detp)]
all(colnames(train_bc) == trainMeta$Sample_ID) 
```

    ## [1] TRUE

``` r
# after dropping samples, update sex vectors
female <- trainMeta$Sex == "Female"
male <- trainMeta$Sex == "Male"

sum(female) 
```

    ## [1] 39

``` r
sum(male) 
```

    ## [1] 48

``` r
# for female Y chromosome, set detp to 0 (they are not failing so not > 0.05)
train_detp[rownames(train_detp) %in% chrY$probeID, female] <- 0

# create a failed probes matrix
# if detP > 0.01 OR bc<3 (NA) OR missing value, sum per cell will be > 0
fp <- (train_detp>0.01) + is.na(train_bc) + is.na(train)
fp <- fp > 0

table(rowSums(fp) > ncol(fp)*0.05) # 1151 failed probes in > 1% samples
```

    ## 
    ##  FALSE   TRUE 
    ## 484361   1151

``` r
fail_probes <- fp[ rowSums(fp) > ncol(fp)*0.05 ,]
dim(fail_probes) 
```

    ## [1] 1151   87

``` r
# remove failed probes in > 5% samples
dim(train) # 485512  
```

    ## [1] 485512     87

``` r
train <- train[!(rownames(train) %in% rownames(fail_probes)),]
dim(train) # 484361 
```

    ## [1] 484361     87

``` r
# identify (& remove) samples with > 5% failed probes
table(colSums(fp) > nrow(fp)*0.05) # zero failing samples
```

    ## 
    ## FALSE 
    ##    87

``` r
## CROSS-HYBRIDIZING ##

# remove CH and polymorphic with zhou
zhou_mask <- zhouAnno %>% filter(MASK_general == TRUE)
dim(zhou_mask)    
```

    ## [1] 60466    57

``` r
dim(train)      
```

    ## [1] 484361     87

``` r
table(rownames(train) %in% rownames(zhou_mask)) 
```

    ## 
    ##  FALSE   TRUE 
    ## 424409  59952

``` r
train <- train[!(rownames(train) %in% rownames(zhou_mask)),]
dim(train) 
```

    ## [1] 424409     87

``` r
# price probes
## read in files
Price450K <- read.csv(here::here("Documents", "excel-files", "annotations", "450KMagdaAnnotation.csv"))
dim(Price450K) 
```

    ## [1] 485512     61

``` r
# pull ids of cross hyb probes
Price_CH_XY <- Price450K[grep("YES", Price450K$XY_Hits), ] 
Price_CH_Auto <- Price450K[grep("YES", Price450K$Autosomal_Hits), ] 
dim(train[which((rownames(train) %in% Price_CH_XY$IlmnID) | rownames(train) %in% Price_CH_Auto$IlmnID), ])
```

    ## [1] 9053   87

``` r
# remove XY cross hybridizing 
dim(train) 
```

    ## [1] 424409     87

``` r
train <- train[!(rownames(train) %in% Price_CH_XY$IlmnID), ] 
dim(train) 
```

    ## [1] 422846     87

``` r
# remove autosomal cross hybridizing
train <- train[!(rownames(train) %in% Price_CH_Auto$IlmnID), ]
dim(train) 
```

    ## [1] 415356     87

``` r
## NON-VARIABLE

# need beta values for this step
betas_train <- getBeta(train)

# function by Rachel Edgar to use on beta values
Variation <- function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}

# call variability in our data
ref_range_train <- lapply(1:nrow(betas_train), function(x) Variation(betas_train[x,]))

# how many probes have a variability range less than 0.05 from the 10th to 90th centile in our dataset
length(which(ref_range_train < 0.05))  
```

    ## [1] 76742

``` r
# create new dataset containing these to overlap with Rachel's list
dim(beta_invariable_train <- betas_train[which(ref_range_train < 0.05),]) 
```

    ## [1] 76742    87

``` r
# now read in database of probes we know to be invariable in many studies, if also invariable in our dataset, remove
Edgar_Invariable <- read.table(here::here("Documents", "EdgarInvariable.txt"), header = TRUE, sep = ",", dec = ".")
dim(Edgar_Invariable) 
```

    ## [1] 101367      3

``` r
head(Edgar_Invariable)
```

    ##   X        CpG   RefRange
    ## 1 1 cg00000108 0.03311130
    ## 2 2 cg00000622 0.03898778
    ## 3 3 cg00000658 0.04814935
    ## 4 4 cg00000721 0.04994386
    ## 5 5 cg00000734 0.02556291
    ## 6 6 cg00000957 0.04776670

``` r
# which of these invariable probes that are in our dataset
dim(beta_invariable_train <- beta_invariable_train[which(rownames(beta_invariable_train) %in% Edgar_Invariable$CpG), ]) 
```

    ## [1] 39600    87

``` r
# filter to match
dim(train) 
```

    ## [1] 415356     87

``` r
train <- train[which(!(rownames(train) %in% rownames(beta_invariable_train))), ]
dim(train) 
```

    ## [1] 375756     87

# Finish QC

Based on our previous observation that GSE73375 has a large
cohort-specific effect that may be due to technical variation (see [Yuan
et al.,
2019](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0296-3#Sec1),
Figure S1), we modified our criteria to remove samples that failed one
or more checks in this dataset such that any samples that failed any
single check would be removed. As a result, we removed GSM1892032,
GSM1892056, and GSM1892059, which failed two checks, but also
GSE1892047, even though it failed a single check.

``` r
# remove XY probes and probes that are not in the EPIC array

dim(train <- train[!rownames(train) %in% chrXY$probeID,]) 
```

    ## [1] 366085     87

``` r
dim(train <- train[rownames(train) %in% rownames(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Manifest),]) 
```

    ## [1] 341281     87

``` r
trainBetas <- getBeta(train)
trainMVals <- lumi::beta2m(trainBetas)

dim(trainBetas) 
```

    ## [1] 341281     87

``` r
range(trainBetas, na.rm=TRUE) 
```

    ## [1] 0 1

``` r
dim(trainMVals) 
```

    ## [1] 341281     87

``` r
range(trainMVals, na.rm=TRUE, finite=TRUE) 
```

    ## [1] -14.97791  14.07848

``` r
sum(is.na(trainBetas)) 
```

    ## [1] 621

``` r
# remove samples that failed qc

dim(train) 
```

    ## [1] 341281     87

``` r
train <- train[,!(colnames(train) %in% c("GSM1892032", "GSM1892047", "GSM1892056", "GSM1892059"))]
dim(train) 
```

    ## [1] 341281     83

# PCA (filtered+raw data)

``` r
pca <- prcomp(t(na.omit(trainBetas)), center=TRUE)
pca_summary <- summary(pca)

loadings <- pca$x
dim(loadings) 
```

    ## [1] 87 87

``` r
rownames(loadings) <- rownames(trainMeta)
head(loadings) 
```

    ##          PC1        PC2        PC3       PC4         PC5       PC6        PC7
    ## 1 -11.862128  17.231064 -10.444621 -3.404631  3.13151863 -1.128859   2.260235
    ## 2  -8.303677  -1.165100  -1.865031  5.029454  4.91235902 -8.750651 -12.922790
    ## 3  -1.437409  -2.554667   1.441985 -9.605772 -0.07209896 -6.410210   4.089717
    ## 4 -10.924293  -3.255333 -10.810382  5.736919  4.54571787 -1.590479  -2.097833
    ## 5 -20.548282   5.317323 -10.792096  6.176827  0.72226342  7.627564   3.169436
    ## 6 -21.367698 -15.851891   2.197787  1.809435 -0.23405543  5.952576 -10.685800
    ##          PC8       PC9       PC10      PC11       PC12       PC13       PC14
    ## 1 -3.9738294 6.4746224 -4.8319048  2.093341 -1.2238647  0.5717859 -0.4474445
    ## 2 -3.0866476 3.6258505 -1.4818535  4.563624  5.1520445  1.7285921 -0.8825030
    ## 3 -4.8586754 6.5187024 -9.2134839  6.633660 -1.8728845  3.1792594  2.6249386
    ## 4 -1.4769340 4.8922445 -9.8666083  3.723540  1.8228194  1.9258482  1.0589702
    ## 5 -0.2648404 5.7634241 -0.8966766  8.986639 -1.2297205  4.6298611  4.5150219
    ## 6  2.3918289 0.7955746 -3.6547647 17.402641 -0.8870158 -3.6940668 -5.3817734
    ##         PC15       PC16       PC17       PC18       PC19       PC20       PC21
    ## 1  0.6493504 -2.6904411  0.4457176 -3.1973885  0.5020024 -1.3477481  2.4664353
    ## 2  1.6789090 -1.9489005  0.3683349 -0.4525661 -2.4944546  2.0488806 -2.8660131
    ## 3  1.4699661 -1.0537008 -0.6309062  3.8593057 -1.9824655 -1.3768145 -0.5839516
    ## 4 -0.5493597 -0.1194807 -3.0523804  0.8825151 -0.7273029  0.6740684  2.6061938
    ## 5 -2.3984402 -6.3946155 -6.0771836 -9.6566675 -0.1033290 -1.8864243  4.8365188
    ## 6 -1.4014888  0.3250363  9.1261381  5.6921644 10.2283733 -6.0810489  2.9084659
    ##         PC22        PC23       PC24        PC25       PC26      PC27       PC28
    ## 1 -1.9906854 -2.12624841  1.8500195 -0.03586114 -1.3325786 -4.681696  0.8962674
    ## 2  0.6523628  2.02879714 -2.4909581  0.61073664 -2.5173823 -3.758324 -1.1607698
    ## 3 -0.2086995 -0.08796067 -1.1000446  1.89573054 -0.1839938 -2.368896  4.7665122
    ## 4  1.0142125 -1.52485162 -1.1695369  2.51780568 -1.9646148 -1.421372  0.9605085
    ## 5  8.7569650 -0.80809616 -1.5351412 10.91910190  6.4443354  2.641001 -0.4497065
    ## 6 -9.9993572  1.37927304 -0.3921779  2.72188187  4.6550639  2.595546 -1.7925805
    ##         PC29       PC30       PC31       PC32       PC33        PC34       PC35
    ## 1  1.1312530  1.4435091  1.4449645 -0.2238490 -1.9969742 -2.02499293  4.6975423
    ## 2  6.0525229  0.4293851  5.5738011  3.0683798 -0.2821479  0.18316433  1.8554732
    ## 3 -0.1065443  7.6080214  0.8426656  1.5604340 -3.3707093 -3.54911742 -1.3583982
    ## 4  0.9894405  3.7183559  2.3241108  0.3070837 -2.2806746 -3.03678285  0.6076417
    ## 5  2.8313755 -5.4442742 -3.2482712 -1.3595064  4.4861177  1.61983196 -1.4895537
    ## 6 -4.3556735  1.8118263 -2.8685214  0.3885169  1.5746626  0.05526557  0.5763697
    ##         PC36       PC37       PC38       PC39        PC40       PC41      PC42
    ## 1 -1.4778517 -0.7151208 -0.2954355 -1.0633314 -2.22354556  6.2074390  1.070118
    ## 2 -0.9665176  4.0547941 -2.2869082 -1.5841091  1.40724503 -0.3948309 -2.172605
    ## 3  4.4694043  0.8704302  2.3827778  3.4100542 -1.85567333 -0.9953614 -6.814053
    ## 4  1.8422000 -1.2060680  0.9038264  2.7618349  0.81442115  2.9882296 -2.552524
    ## 5  1.0016599  5.0988413 -2.7849449 -0.8594524 -0.09931148 -4.7503534  1.839019
    ## 6 -1.9201883  1.3974375 -0.8472293 -3.0174997  2.01640315 -0.1893929  1.335534
    ##          PC43       PC44       PC45       PC46         PC47      PC48
    ## 1  1.77180002  0.1235929 -1.9752839 -1.0890222 -2.143301594 -1.399941
    ## 2 -2.12337983  6.8248934 -2.2773634  0.2166827  4.968281941  1.822640
    ## 3 -3.66721520  1.4671714  1.5042392 -0.3346529 -0.004790732  1.377123
    ## 4 -0.01131173 -2.5120620  2.8456280  2.7770156  0.041772555 -3.111091
    ## 5  0.63801311  1.2264507  0.1655110 -0.1171945 -0.080845823  1.077356
    ## 6  1.52242820 -0.9048020 -0.4885068  0.6536793 -0.286972746 -1.858874
    ##          PC49       PC50       PC51       PC52       PC53       PC54       PC55
    ## 1 -0.57653650  0.4666478  2.6018915  2.7676876  2.2162061 -0.6885552 -0.4647893
    ## 2  4.75735713  3.6620278 -2.8877333 -6.5286126 -0.1186810  6.0635785 -5.6461149
    ## 3 -0.05791213 -0.8257687 -4.6867810  5.5353115  1.8070917 -0.7836241  4.3873070
    ## 4 -4.60447185  1.6671510  0.7050339 -5.3528383 -3.0218530 -3.7509671  0.6074960
    ## 5 -0.42702863  0.1744645 -0.3295374  1.1800795 -0.8201858 -0.5313423  1.0004628
    ## 6  0.33208667 -0.4192195 -0.6888417 -0.6046416  0.2774852  0.2501045  0.2140663
    ##          PC56       PC57       PC58       PC59        PC60       PC61
    ## 1 -1.92637619 -6.2370335  4.1266721  3.5687037  2.75629154  3.1450912
    ## 2 -0.07944708  0.3716956 -0.7732280  2.1211646  0.08575383  3.6741368
    ## 3 -6.33479313  2.9887290 -3.1129544  3.4238458 -1.69708063 -2.1531610
    ## 4  7.96935638 -5.7778285  0.4072716 -5.6240745  2.29813229  1.2926851
    ## 5 -0.40081194  0.9561776 -1.0688204  0.6251475 -0.12269896 -1.5214949
    ## 6 -0.69093418  0.5368902 -0.2494238 -0.7145838 -0.06111362  0.6762006
    ##         PC62       PC63       PC64       PC65        PC66        PC67
    ## 1  2.4289444 -0.4830342 -6.5135697 -3.7957674 -3.72797574  1.67252307
    ## 2  1.6308464 -0.2599255  1.5934253 -0.5677761  4.00001943  0.86665085
    ## 3 -0.6484177 -2.2735460  1.0440746  0.5064854  0.74588599 -4.69906449
    ## 4 -3.2862100  4.7671568 -0.6432479  2.4498535 -0.02661778  0.04091004
    ## 5 -0.3622046 -0.5443908 -1.0008203  0.2330153 -0.22502761  0.63171139
    ## 6  0.5109043 -0.4650242  0.1282726 -0.1976373  0.02150952 -0.23179671
    ##         PC68        PC69       PC70       PC71        PC72       PC73
    ## 1  2.4464517  6.15575712 -1.2820355 -5.2983412 -0.38546723 -1.4426516
    ## 2  2.5214511 -0.27361910 -0.1437868  0.6213541 -0.64435250 -0.5123122
    ## 3 -0.6537093 -0.81708409 -0.9425759  0.5703949 -0.06777019  2.2161245
    ## 4 -5.1044870 -3.71581060 -0.2238732  1.3155964 -1.69910406  1.9773606
    ## 5  0.1621482 -0.01672417 -0.4313230  0.4861293  0.27514866 -1.0803352
    ## 6 -0.3820613 -0.28526265  0.2325251  0.4621548 -0.37333237 -0.2106025
    ##          PC74       PC75        PC76         PC77        PC78        PC79
    ## 1 -2.21971697  2.4148760  1.25526920 -1.169957597 -0.20234239 -0.33033974
    ## 2  0.94639356  0.4388258  1.56720962 -0.437814791  0.09346123  0.22839956
    ## 3 -1.12549671  0.5112281  0.13893496 -1.085326408 -0.33843668 -0.23946379
    ## 4  1.05271973 -0.2231722  0.31049549 -1.609035609  0.42843593  2.00718626
    ## 5 -0.02818242 -0.1559652  0.01549557 -0.210730873 -0.06107676 -0.09971771
    ## 6 -0.44533395  0.4005139 -0.43732436 -0.007000956  0.37565096  0.03109227
    ##         PC80       PC81       PC82        PC83        PC84       PC85
    ## 1 -1.3434257  0.9584146 -0.5043083  1.21476374 -0.63991043  1.8924520
    ## 2 -1.0731194  0.7797537  0.1218304  0.06356568  0.39646466  0.1193992
    ## 3  1.6958523 -0.2998315  0.8939818  0.21596367  0.42271088  0.9467526
    ## 4  0.1454098 -1.2125188  0.3487004  0.92514164  0.16214901 -0.6735271
    ## 5 -0.2105762  0.7458720 -0.3541831 -0.25783787  0.01025978  0.1819603
    ## 6 -0.4347175  0.5148263  0.1817878 -0.26929732 -0.28033540  0.1330625
    ##            PC86          PC87
    ## 1 -0.6560511462  2.280084e-13
    ## 2  0.0005625914 -6.302608e-14
    ## 3  0.1479122184  8.037242e-14
    ## 4  0.1447169735 -1.988161e-14
    ## 5  0.0852696683  2.099958e-13
    ## 6 -0.0605243805 -6.414666e-14

``` r
fviz_eig(pca, geom = "bar", bar_width = 0.5, 
              addlabels = TRUE, hjust = -0.3,
              barfill = "#BF98A0", barcolor = "#BF98A0") + 
  ggtitle("PCA Variances of Training Data (Filtered, Pre-Normalization)") +
  theme_minimal() +
  labs(x = "Principal Components", y = "% of explained variances")
```

![](00_Training_files/figure-markdown_github/pca%20filt-1.png)

``` r
# rename columns
pca_scores <- pca$x %>% as_tibble() %>% mutate(Sample_Name = trainMeta$Sample_ID) 

# r squared
pca_cor_rs <- lmmatrix(dep = pca_scores[,1:10],
                        ind = as.data.frame(trainMeta) %>% 
                          dplyr::select(GEO_Dataset,
                                        Sex, 
                                        PE,
                                        Class,
                                        GA_RPC,
                                        Predicted_ethnicity_nothresh),
                       metric = 'Rsquared')

pca_plot_rs <- pca_cor_rs %>% as.data.frame() %>% 
  mutate(indep = rownames(pca_cor_rs)) %>%
 pivot_longer(c(-indep), names_to = "PC", values_to = "rsqr" ) %>%
  mutate(
  PC = factor(PC, levels = colnames(pca_cor_rs))) %>% as_tibble()

# Plot PCA 
p1_r <- ggplot(pca_plot_rs, aes(x = PC, y = indep, fill = rsqr)) +
  geom_tile(col = 'lightgrey') + theme_classic() +
  scale_x_discrete(expand = c(0, 0), labels = 1:20) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = '', fill = 'R Squared') +
  scale_fill_gradientn(colors = c("#ff4d00", "#f98010", "#f09623", "#e9b448")) +
  ggtitle("PCA - Training Data (Filtered, Pre-Normalization)") 

p1_r
```

![](00_Training_files/figure-markdown_github/pca%20filt-2.png)

``` r
as.data.frame(pca_cor_rs) %>%
  filter(PC1 > 0.5) # Dataset has an R-squared of > 0.5 with PC1, no other variables do
```

    ##  [1] PC1  PC2  PC3  PC4  PC5  PC6  PC7  PC8  PC9  PC10
    ## <0 rows> (or 0-length row.names)

``` r
# p value
pca_cor_pval <- lmmatrix(dep = pca_scores[,1:10],
                        ind = as.data.frame(trainMeta) %>% 
                          dplyr::select(GEO_Dataset,
                                        Sex, 
                                        PE,
                                        Class,
                                        GA_RPC,
                                        Predicted_ethnicity_nothresh),
                       metric = 'Pvalue')

pca_plot_pval <- pca_cor_pval %>% as.data.frame() %>% 
  mutate(indep = rownames(pca_cor_pval)) %>%
 pivot_longer(c(-indep), names_to = "PC", values_to = "pval" ) %>%
  mutate(
  PC = factor(PC, levels = colnames(pca_cor_pval))) %>% as_tibble()

# Plot PCA 
p1_r <- ggplot(pca_plot_pval, aes(x = PC, y = indep, fill = pval)) +
  geom_tile(col = 'lightgrey') + theme_classic() +
  scale_x_discrete(expand = c(0, 0), labels = 1:20) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = '', fill = 'P Value') +
  scale_fill_gradientn(colors = c("#0047AB", "#6495ED", "#ADD8E6")) +
  ggtitle("PCA - Training Data (Filtered, Pre-Normalization)") 

p1_r # also strongly associated to Dataset, Sex, Ethnicity, but not so much class
```

![](00_Training_files/figure-markdown_github/pca%20filt-3.png)

``` r
as.data.frame(pca_cor_pval) %>%
  filter(PC1 < 0.01) # Dataset is the only variable has a p-val < 0.01 to PC1 *and* PC2
```

    ##                      PC1       PC2          PC3        PC4        PC5       PC6
    ## GEO_Dataset 5.787261e-12 0.1201757 8.193549e-09 2.5728e-06 0.03340006 0.5966683
    ##                   PC7        PC8          PC9       PC10
    ## GEO_Dataset 0.2965272 0.05821421 2.592911e-05 0.00734379

``` r
as.data.frame(pca_cor_pval) %>%
  filter(PC1 < 0.05) # No other variables are actually significantly associated to PC1 with a p-val of less than 0.05
```

    ##                      PC1       PC2          PC3        PC4        PC5       PC6
    ## GEO_Dataset 5.787261e-12 0.1201757 8.193549e-09 2.5728e-06 0.03340006 0.5966683
    ##                   PC7        PC8          PC9       PC10
    ## GEO_Dataset 0.2965272 0.05821421 2.592911e-05 0.00734379

``` r
# plot betas by type
probeDesign <- data.frame(Type = getProbeType(train, withColor = FALSE))
probeID <- rownames(train)
probeDesign$Name <- probeID
any(is.na(probeDesign)) 
```

    ## [1] FALSE

``` r
plotBetasByType(train[,1], main="Training Data - Raw (Filtered)")
```

![](00_Training_files/figure-markdown_github/pca%20filt-4.png)

``` r
# clear environment
rm(pca)
rm(pca_cor_pval)
rm(pca_cor_rs)
rm(pca_plot_pval)
rm(pca_plot_rs)
rm(pca_scores)
rm(pca_summary)
gc()
```

    ##              used   (Mb) gc trigger    (Mb)   max used    (Mb)
    ## Ncells   17031694  909.6   26591544  1420.2   26591544  1420.2
    ## Vcells 1103360201 8418.0 1860859972 14197.3 1781159761 13589.2

# Normalization

``` r
dim(train) # 341,281 by 83
```

    ## [1] 341281     83

``` r
## normalize
set.seed(2022)
tictoc::tic()
trainBMIQ <- wateRmelon::BMIQ(train)
```

    ## Warning: package 'RPMM' was built under R version 4.1.3

    ## Warning: package 'cluster' was built under R version 4.1.3

``` r
tictoc::toc() 
```

    ## 414.01 sec elapsed

``` r
# 398.54 sec
```

# PCA (filtered+normalized data)

``` r
plotBetasByType(as.matrix(trainBMIQ)[,1], probeTypes = as.data.frame(probeDesign), main="Training Data - BMIQ Normalized + Filtered")
```

![](00_Training_files/figure-markdown_github/pca%20filt_norm-1.png)

``` r
dim(trainBMIQ)
```

    ## [1] 341281     83

``` r
dim(trainMeta)
```

    ## [1] 87 30

``` r
trainMeta <- trainMeta %>% filter(Sample_ID %in% colnames(trainBMIQ))

all(colnames(trainBMIQ) == trainMeta$Sample_ID)
```

    ## [1] TRUE

``` r
## run PCA again
pca <- prcomp(t(na.omit(trainBMIQ)), center = TRUE)
pca_summary <- summary(pca)

loadings <- pca$x
dim(loadings) 
```

    ## [1] 83 83

``` r
rownames(loadings) <- trainMeta$Sample_ID
head(loadings) 
```

    ##                   PC1       PC2        PC3        PC4        PC5        PC6
    ## GSM1388241 -0.5966416 24.797918  -2.136828  3.8703274 -0.2048727  0.8914226
    ## GSM1388242  1.9926001  0.862125  -3.777975 -4.5674421 10.6297191 15.2988610
    ## GSM1388243 -4.4085977 -2.988381   3.819879  8.5066095 -5.6347196  4.0353714
    ## GSM1388259  4.7539747  5.484147 -11.824309 -3.8843330  0.7548755  3.5017748
    ## GSM1388262  5.8632650 14.039268  -8.204104 -3.0424204  1.6997968 -9.6426769
    ## GSM1388266 27.5333200 -9.503669  -3.331792 -0.9459334  2.8106612  3.2115292
    ##                   PC7        PC8        PC9      PC10      PC11       PC12
    ## GSM1388241 -5.8106910  0.3124159  -8.438158 -3.698408  1.005354 -1.4689790
    ## GSM1388242  0.1309417  2.3420614  -7.228406  1.349075 -4.177903 -0.4893841
    ## GSM1388243 -6.7923494 -2.0690571 -13.305194 -3.373267  1.186385 -3.3690004
    ## GSM1388259  0.2012108 -4.7634912 -11.486700 -2.814492 -2.049200 -1.1818918
    ## GSM1388262  2.5950177 -0.6017148 -12.190026 -1.939474  1.639502  0.4391159
    ## GSM1388266 14.4510683  2.2817950 -13.853694 11.224227 10.013035  3.5481949
    ##                  PC13       PC14      PC15       PC16       PC17        PC18
    ## GSM1388241  0.9913692 -1.0732261  1.075526 -4.5388658  0.4814486  -1.2901250
    ## GSM1388242  0.5687350 -0.2500278  4.092640  0.2597929  1.5724564   0.8872229
    ## GSM1388243 -3.5604005 -2.1140688  1.543135  2.9643777  3.7672127   2.7207674
    ## GSM1388259  0.1193609  1.0797710 -1.353116  2.3495402  2.9693758  -1.7414593
    ## GSM1388262 -5.7260019  2.1578866  2.098083 -5.9945832  5.3614076 -17.4690716
    ## GSM1388266  0.0237771 -2.3882877  3.921308 -2.2747025 -7.7569626   8.1655147
    ##                   PC19       PC20       PC21       PC22       PC23      PC24
    ## GSM1388241  0.80557231 -2.9119355  2.4280320  2.1774171 -0.5632014 -1.638398
    ## GSM1388242 -3.12002530  3.0733395  0.5854590 -0.2714133 -0.2887433  3.006159
    ## GSM1388243  0.90651893  0.8076546  0.9267566  2.1887599 -1.0178861  1.671895
    ## GSM1388259  0.06782724 -1.2848133  0.6003416  3.1889027  1.6472046  1.880886
    ## GSM1388262  0.79204895  0.8346637 -7.5480667 -5.5188462  6.2972577  8.094009
    ## GSM1388266 11.00869190 -7.1954511  7.9644562 -3.2134858 -0.6633405 -0.242198
    ##                 PC25      PC26       PC27       PC28       PC29       PC30
    ## GSM1388241 -5.456323  1.107155  2.4708446  0.8291942 -2.8297314 -0.7787061
    ## GSM1388242 -3.983834  3.977930  3.1514607  2.9189660 -1.2984581  4.9446107
    ## GSM1388243 -4.038549 -5.320636 -2.8352510  4.4857335 -5.3125615 -2.2340419
    ## GSM1388259 -1.810942 -1.103026 -0.1404099  3.9231736 -3.2049747  0.9617049
    ## GSM1388262  6.386270 -2.139798  1.4787695 -4.7972929  0.9534349 -3.1063815
    ## GSM1388266  6.280030 -6.611774  0.5828661  0.7272589  1.1756647 -2.0329737
    ##                  PC31        PC32       PC33        PC34       PC35       PC36
    ## GSM1388241  2.5585946  1.03977764  2.1284568  3.86747357 -3.8565616 -1.8948518
    ## GSM1388242  4.5666910  2.31730767  1.7406810  0.08316909  4.5975507 -3.2922725
    ## GSM1388243  4.7973369 -0.04434816  6.1164499 -2.78474784  0.7829240  3.6661130
    ## GSM1388259  1.2244910 -0.39752384  3.0474495  2.42686768  0.7273056  1.7688000
    ## GSM1388262  0.3751222  5.27999028 -0.8150882 -2.90831451  4.4118400 -0.2271115
    ## GSM1388266 -0.5575302  0.94826413  0.5271929  0.41684650 -1.0241820 -1.7731255
    ##                  PC37       PC38       PC39       PC40        PC41       PC42
    ## GSM1388241  0.1102995  0.9457548  4.0174271 -5.0425330 -0.08738941 -3.3510346
    ## GSM1388242 -2.8245309  0.4971691  2.4791012  8.1607617 -2.91242759 -1.7519642
    ## GSM1388243  6.5528622  3.2467659  0.5299488  2.3735092 -1.75912801  4.4156461
    ## GSM1388259  4.4718401 -1.1081165  3.3252732 -1.9327297  0.51163890  2.6100618
    ## GSM1388262 -4.0710595  1.3947425 -3.4825478  1.9134704  0.13894104 -2.0476527
    ## GSM1388266 -3.7907283 -1.5036599  0.7582995  0.5369075  1.44050507 -0.5735387
    ##                   PC43        PC44        PC45       PC46       PC47       PC48
    ## GSM1388241 -1.13707732 -1.84537950  2.26780271  2.4542381  2.7720213  4.1303146
    ## GSM1388242 -2.75855082 -0.08706722 -3.65407035 -5.4976623 -7.6172177 -0.4239709
    ## GSM1388243  0.06911233 -0.54293456 -0.71260441 -1.2916116  0.7772960 -0.5738529
    ## GSM1388259  3.40690978  1.12004823 -1.19427265  2.6031914  3.6822786 -0.5258218
    ## GSM1388262  0.70155906  0.17139212  0.32995359 -0.7748844  0.2154959 -0.5940791
    ## GSM1388266  0.38653158  1.45135383 -0.05635601  2.0187102 -0.7003429 -0.6373096
    ##                  PC49       PC50        PC51       PC52       PC53        PC54
    ## GSM1388241 -0.1014985 -1.5830575 -1.67081923 -0.2462929 -1.8475744 -6.98739303
    ## GSM1388242 -1.2946619  7.0524973  3.97744742 -6.3041787 -2.2614410 -0.57711873
    ## GSM1388243 -3.7941263 -6.9582886  0.54859445  0.2622000  5.9332630 -4.12706571
    ## GSM1388259 -1.5225935  4.9398349 -0.47144868  1.9234475 -5.6501384  5.41546986
    ## GSM1388262  0.1849848 -1.5969880  0.35145342  0.8263612  0.9668774 -0.08349462
    ## GSM1388266 -0.4757602  0.2618105 -0.04518976  0.3286175  0.3886136 -0.17811159
    ##                  PC55       PC56       PC57       PC58       PC59       PC60
    ## GSM1388241  0.3817185  3.8957806  6.8063288 -0.7212147 -8.7740942 -0.1675551
    ## GSM1388242 -1.8409473 -0.5577046  3.9991277 -3.5098152  0.7090711  1.1005578
    ## GSM1388243 -1.5830436 -6.2595266 -1.6668565  1.8571332  4.5159921 -1.2542235
    ## GSM1388259  8.5597867  5.2642780 -1.8792294 -0.9734167  0.1675779  7.4088652
    ## GSM1388262  0.2667805 -1.1172251 -0.5322325  1.5080783 -0.4995387 -0.9167520
    ## GSM1388266 -0.6928714 -0.2253038 -0.7891829 -1.1232297 -0.1234114  0.3924916
    ##                  PC61       PC62        PC63       PC64        PC65       PC66
    ## GSM1388241 -0.4961262  2.6749930  3.06808924  2.6479690  0.06436105  5.4436758
    ## GSM1388242 -1.6160642 -2.0763437 -3.66534799  1.0602725 -0.02317186  2.4225997
    ## GSM1388243 -2.7954670 -2.7440179  1.48851003 -4.9805631  0.33754105  0.6130625
    ## GSM1388259  5.9348712  0.5015109 -0.00942793 -1.6465817 -2.36083833 -7.2410999
    ## GSM1388262 -0.4369721  0.8344967 -0.26090836  0.2410517 -0.48808657  0.1188297
    ## GSM1388266 -0.8196476 -0.2432387  0.13029984  0.2206395 -0.09751085 -0.1184362
    ##                    PC67       PC68       PC69        PC70       PC71
    ## GSM1388241 -6.632369741 -1.6508587 -1.8292437  2.23969662  1.4316465
    ## GSM1388242  0.992023995 -0.6861065 -0.1545418 -1.84432113  1.0197031
    ## GSM1388243 -0.002644713  1.1795888  1.2897219  1.19421518  0.4208884
    ## GSM1388259  0.827564902  2.1848175  1.2928050 -0.69277959  0.8792031
    ## GSM1388262 -0.185897400  0.3031881 -0.8844711  0.07249975 -0.4190479
    ## GSM1388266  0.621285306 -0.1787277 -0.7963839  0.45747414  0.1204492
    ##                   PC72        PC73       PC74       PC75       PC76        PC77
    ## GSM1388241 -1.06017656 -0.34189091 -0.2981583  0.2918537  0.4541353  0.93191961
    ## GSM1388242 -1.40990819 -0.40375406  0.4442235  0.1553045  0.7565776  0.97594653
    ## GSM1388243 -0.84286684  0.09249829 -1.2831208  0.9904707 -0.6462082 -2.00379918
    ## GSM1388259 -1.20744157  1.69746442  0.1191747 -1.8766140 -1.2535686  0.20201168
    ## GSM1388262  0.04897722  0.06637090  0.0741864  0.0542601  0.7001201  0.07350267
    ## GSM1388266  0.53080149  0.20546244  0.2259146 -0.3275368  0.6933152  0.33172843
    ##                  PC78       PC79       PC80        PC81        PC82
    ## GSM1388241 -0.4113412  0.9926320 -1.7456765  2.20413476 -0.50194600
    ## GSM1388242  0.3855084  0.3104372  0.6383076  0.45172039  0.01597822
    ## GSM1388243  0.8807227  0.3901961 -0.3569140  0.89649782  0.19930496
    ## GSM1388259 -0.2783207  0.8508070  1.0766373 -0.81876984 -0.14461389
    ## GSM1388262 -0.2829751 -0.3327596 -0.0379121  0.17775972  0.13635983
    ## GSM1388266  0.1375513 -0.3007597 -0.4068040 -0.04019409  0.06968866
    ##                     PC83
    ## GSM1388241 -1.172687e-13
    ## GSM1388242  6.935690e-14
    ## GSM1388243 -7.372654e-14
    ## GSM1388259 -3.229434e-14
    ## GSM1388262  9.276598e-15
    ## GSM1388266 -6.542189e-14

``` r
fviz_eig(pca, geom = "bar", bar_width = 0.5, 
              addlabels = TRUE, hjust = -0.3,
              barfill = "#BF98A0", barcolor = "#BF98A0") + 
  ggtitle("PCA Variances of Training Data (Filtered, BMIQ Normalized)") +
  theme_minimal() +
  labs(x = "Principal Components", y = "% of explained variances")
```

![](00_Training_files/figure-markdown_github/pca%20filt_norm-2.png)

``` r
# rename columns
pca_scores <- pca$x %>% as_tibble() %>% mutate(Sample_Name = trainMeta$Sample_ID) 

# r squared
pca_cor_rs <- lmmatrix(dep = pca_scores[,1:10],
                        ind = as.data.frame(trainMeta) %>% 
                          dplyr::select(GEO_Dataset,
                                        Sex, 
                                        PE,
                                        Class,
                                        GA_RPC,
                                        Predicted_ethnicity_nothresh),
                       metric = 'Rsquared')

pca_plot_rs <- pca_cor_rs %>% as.data.frame() %>% 
  mutate(indep = rownames(pca_cor_rs)) %>%
 pivot_longer(c(-indep), names_to = "PC", values_to = "rsqr" ) %>%
  mutate(
  PC = factor(PC, levels = colnames(pca_cor_rs))) %>% as_tibble()

# Plot PCA 
p1_r <- ggplot(pca_plot_rs, aes(x = PC, y = indep, fill = rsqr)) +
  geom_tile(col = 'lightgrey') + theme_classic() +
  scale_x_discrete(expand = c(0, 0), labels = 1:20) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = '', fill = 'R Squared') +
  scale_fill_gradientn(colors = c("#ff4d00", "#f98010", "#f09623", "#e9b448")) +
  ggtitle("PCA - Training Data (Filtered, BMIQ Normalized)") 

p1_r
```

![](00_Training_files/figure-markdown_github/pca%20filt_norm-3.png)

``` r
as.data.frame(pca_cor_rs) %>%
  filter(PC1 > 0.5) 
```

    ##  [1] PC1  PC2  PC3  PC4  PC5  PC6  PC7  PC8  PC9  PC10
    ## <0 rows> (or 0-length row.names)

``` r
as.data.frame(pca_cor_rs) %>%
  filter(PC1 > 0.4) 
```

    ##                   PC1       PC2      PC3      PC4       PC5        PC6
    ## GEO_Dataset 0.4985784 0.3050988 0.017986 0.239499 0.1647143 0.01697783
    ##                   PC7        PC8       PC9       PC10
    ## GEO_Dataset 0.0962021 0.05085486 0.6058007 0.05907701

``` r
# p value
pca_cor_pval <- lmmatrix(dep = pca_scores[,1:10],
                        ind = as.data.frame(trainMeta) %>% 
                          dplyr::select(GEO_Dataset,
                                        Sex, 
                                        PE,
                                        Class,
                                        GA_RPC,
                                        Predicted_ethnicity_nothresh),
                       metric = 'Pvalue')

pca_plot_pval <- pca_cor_pval %>% as.data.frame() %>% 
  mutate(indep = rownames(pca_cor_pval)) %>%
 pivot_longer(c(-indep), names_to = "PC", values_to = "pval" ) %>%
  mutate(
  PC = factor(PC, levels = colnames(pca_cor_pval))) %>% as_tibble()

# Plot PCA 
p1_r <- ggplot(pca_plot_pval, aes(x = PC, y = indep, fill = pval)) +
  geom_tile(col = 'lightgrey') + theme_classic() +
  scale_x_discrete(expand = c(0, 0), labels = 1:20) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = '', fill = 'P Value') +
  scale_fill_gradientn(colors = c("#0047AB", "#6495ED", "#ADD8E6")) +
  ggtitle("PCA - Training Data (Filtered, BMIQ Normalized)") 

p1_r
```

![](00_Training_files/figure-markdown_github/pca%20filt_norm-4.png)

``` r
as.data.frame(pca_cor_pval) %>%
  filter(PC1 < 0.05) 
```

    ##                      PC1          PC2         PC3          PC4         PC5
    ## GEO_Dataset 7.362520e-12 2.316133e-06 0.695523927 7.312562e-05 0.002516051
    ## PE          1.451264e-02 5.253989e-04 0.006215008 1.368382e-01 0.448479563
    ## Class       1.451264e-02 5.253989e-04 0.006215008 1.368382e-01 0.448479563
    ##                      PC6        PC7       PC8          PC9      PC10
    ## GEO_Dataset 7.146380e-01 0.04517812 0.2458038 6.028803e-16 0.1838360
    ## PE          8.346316e-05 0.33965151 0.5184442 2.277890e-01 0.1382215
    ## Class       8.346316e-05 0.33965151 0.5184442 2.277890e-01 0.1382215

``` r
as.data.frame(pca_cor_pval) %>%
  filter(PC1 < 0.01) 
```

    ##                     PC1          PC2       PC3          PC4         PC5
    ## GEO_Dataset 7.36252e-12 2.316133e-06 0.6955239 7.312562e-05 0.002516051
    ##                  PC6        PC7       PC8          PC9     PC10
    ## GEO_Dataset 0.714638 0.04517812 0.2458038 6.028803e-16 0.183836

# Batch correction

``` r
## combat correct

# create model matrix to protect variables
combat_mod <- model.matrix(~as.factor(Class) + 
                             as.factor(Sex) + 
                             GA_RPC + 
                             Prob_Asian + 
                             Prob_African, 
                           data = trainMeta)

head(combat_mod) 
```

    ##   (Intercept) as.factor(Class)Normotensive as.factor(Sex)Male   GA_RPC
    ## 1           1                            1                  1 35.87492
    ## 2           1                            1                  0 36.28301
    ## 3           1                            1                  1 36.07620
    ## 4           1                            0                  1 33.62615
    ## 5           1                            0                  1 29.30807
    ## 6           1                            0                  0 30.92914
    ##    Prob_Asian Prob_African
    ## 1 0.012693618    0.9574839
    ## 2 0.423882904    0.5489429
    ## 3 0.016452763    0.9486234
    ## 4 0.002674077    0.9965530
    ## 5 0.005656520    0.9895070
    ## 6 0.006810293    0.9892902

``` r
sum(is.na(trainBMIQ)) # 611
```

    ## [1] 611

``` r
# impute NAs to be able to run combat
set.seed(2022)

cluster <- makeCluster(16)
registerDoParallel(cluster)

tic()
trainBMIQ_noNA <- impute.knn(as.matrix(trainBMIQ), maxp = 15000)$data
```

    ## Cluster size 341281 broken into 171729 169552 
    ## Cluster size 171729 broken into 68801 102928 
    ## Cluster size 68801 broken into 37052 31749 
    ## Cluster size 37052 broken into 17691 19361 
    ## Cluster size 17691 broken into 9962 7729 
    ## Done cluster 9962 
    ## Done cluster 7729 
    ## Done cluster 17691 
    ## Cluster size 19361 broken into 9663 9698 
    ## Done cluster 9663 
    ## Done cluster 9698 
    ## Done cluster 19361 
    ## Done cluster 37052 
    ## Cluster size 31749 broken into 15585 16164 
    ## Cluster size 15585 broken into 7248 8337 
    ## Done cluster 7248 
    ## Done cluster 8337 
    ## Done cluster 15585 
    ## Cluster size 16164 broken into 10240 5924 
    ## Done cluster 10240 
    ## Done cluster 5924 
    ## Done cluster 16164 
    ## Done cluster 31749 
    ## Done cluster 68801 
    ## Cluster size 102928 broken into 63383 39545 
    ## Cluster size 63383 broken into 33534 29849 
    ## Cluster size 33534 broken into 17241 16293 
    ## Cluster size 17241 broken into 5492 11749 
    ## Done cluster 5492 
    ## Done cluster 11749 
    ## Done cluster 17241 
    ## Cluster size 16293 broken into 14450 1843 
    ## Done cluster 14450 
    ## Done cluster 1843 
    ## Done cluster 16293 
    ## Done cluster 33534 
    ## Cluster size 29849 broken into 15151 14698 
    ## Cluster size 15151 broken into 3829 11322 
    ## Done cluster 3829 
    ## Done cluster 11322 
    ## Done cluster 15151 
    ## Done cluster 14698 
    ## Done cluster 29849 
    ## Done cluster 63383 
    ## Cluster size 39545 broken into 22224 17321 
    ## Cluster size 22224 broken into 12553 9671 
    ## Done cluster 12553 
    ## Done cluster 9671 
    ## Done cluster 22224 
    ## Cluster size 17321 broken into 9787 7534 
    ## Done cluster 9787 
    ## Done cluster 7534 
    ## Done cluster 17321 
    ## Done cluster 39545 
    ## Done cluster 102928 
    ## Done cluster 171729 
    ## Cluster size 169552 broken into 110814 58738 
    ## Cluster size 110814 broken into 27441 83373 
    ## Cluster size 27441 broken into 15278 12163 
    ## Cluster size 15278 broken into 12171 3107 
    ## Done cluster 12171 
    ## Done cluster 3107 
    ## Done cluster 15278 
    ## Done cluster 12163 
    ## Done cluster 27441 
    ## Cluster size 83373 broken into 53129 30244 
    ## Cluster size 53129 broken into 24512 28617 
    ## Cluster size 24512 broken into 13130 11382 
    ## Done cluster 13130 
    ## Done cluster 11382 
    ## Done cluster 24512 
    ## Cluster size 28617 broken into 15431 13186 
    ## Cluster size 15431 broken into 7613 7818 
    ## Done cluster 7613 
    ## Done cluster 7818 
    ## Done cluster 15431 
    ## Done cluster 13186 
    ## Done cluster 28617 
    ## Done cluster 53129 
    ## Cluster size 30244 broken into 22156 8088 
    ## Cluster size 22156 broken into 10676 11480 
    ## Done cluster 10676 
    ## Done cluster 11480 
    ## Done cluster 22156 
    ## Done cluster 8088 
    ## Done cluster 30244 
    ## Done cluster 83373 
    ## Done cluster 110814 
    ## Cluster size 58738 broken into 31030 27708 
    ## Cluster size 31030 broken into 15848 15182 
    ## Cluster size 15848 broken into 10081 5767 
    ## Done cluster 10081 
    ## Done cluster 5767 
    ## Done cluster 15848 
    ## Cluster size 15182 broken into 8837 6345 
    ## Done cluster 8837 
    ## Done cluster 6345 
    ## Done cluster 15182 
    ## Done cluster 31030 
    ## Cluster size 27708 broken into 14051 13657 
    ## Done cluster 14051 
    ## Done cluster 13657 
    ## Done cluster 27708 
    ## Done cluster 58738 
    ## Done cluster 169552

``` r
toc() 
```

    ## 16.45 sec elapsed

``` r
sum(is.na(trainBMIQ_noNA)) # 0
```

    ## [1] 0

``` r
stopCluster(cluster)
stopImplicitCluster()

# convert to mvals
mvals <- lumi::beta2m(trainBMIQ_noNA)
dim(mvals) 
```

    ## [1] 341281     83

``` r
range(mvals) 
```

    ## [1] -26.14073  24.77860

``` r
all(colnames(mvals) == trainMeta$Sample_ID) 
```

    ## [1] TRUE

``` r
tic()
combat_mvals <- sva::ComBat(dat=mvals,
                            batch=trainMeta$GEO_Dataset,
                            mod=combat_mod,
                            par.prior = T,
                            prior.plots = F)
toc() 
```

    ## 58.59 sec elapsed

``` r
all(trainMeta$Sample_ID == colnames(trainBMIQ))  
```

    ## [1] TRUE

``` r
all(trainMeta$Sample_ID == colnames(mvals)) 
```

    ## [1] TRUE

``` r
dim(combat_mvals) 
```

    ## [1] 341281     83

``` r
range(combat_mvals) 
```

    ## [1] -19.22162  22.22131

``` r
trainComBat <- lumi::m2beta(combat_mvals)

dim(trainComBat) 
```

    ## [1] 341281     83

``` r
range(trainComBat) 
```

    ## [1] 1.635743e-06 9.999998e-01

# PCA (filtered+normalized+batch corrected data)

``` r
plotBetasByType(as.matrix(trainComBat)[,1], probeTypes = as.data.frame(probeDesign), main="Training Data - Batch Corrected + BMIQ Normalized + Filtered") 
```

![](00_Training_files/figure-markdown_github/pca%20filt_norm_batchcorrected-1.png)

``` r
# now check pca
pca <- prcomp(t(na.omit(trainComBat)), n = 10)
```

    ## Warning: In prcomp.default(t(na.omit(trainComBat)), n = 10) :
    ##  extra argument 'n' will be disregarded

``` r
pca_summary <- summary(pca)

loadings <- pca$x
dim(loadings) 
```

    ## [1] 83 83

``` r
rownames(loadings) <- rownames(trainMeta)
head(loadings) 
```

    ##           PC1         PC2         PC3        PC4         PC5         PC6
    ## 1 -21.1386816   5.2054131  -0.5166723  0.9120434  -5.6014310   7.3446481
    ## 2  -1.8013097  -3.8940127  -1.9011885 -3.5723923  16.6085017  -3.5120567
    ## 3  -0.4670007 -13.6073757   7.7395143 -4.6430490  -4.9716163   7.6406378
    ## 4  -3.0340684  -0.2236467 -10.5659364 -5.2578817  -0.4621531  -0.2366149
    ## 5  -9.2911654   6.3584434  -7.6196127 -1.9115768 -11.4957820  -7.8818635
    ## 6  23.0680954  12.5449739  -2.9081244 -0.5267981   2.7292622 -10.8770107
    ##         PC7        PC8        PC9      PC10       PC11       PC12       PC13
    ## 1 -2.206391  2.0668549  -4.433549  7.282831 -3.4220633  3.0597038 -0.5116367
    ## 2  3.863458  5.9376477  -2.449366  3.138519 -1.9748568  1.3190813 -4.5953687
    ## 3 -3.040114  6.5690974  -7.787524  5.036557 -2.1108239 -3.1855452 -1.5841164
    ## 4 -3.103837 10.6191946  -4.439183  5.489028 -0.9375366 -2.1639377  1.6652839
    ## 5 -3.164166  2.2318920  -5.237035  5.602731 -1.2571816 -4.2949973 -3.9433207
    ## 6 -4.642267  0.3645654 -15.598791 -2.341621  0.5073873 -0.5270459  3.2280717
    ##          PC14      PC15       PC16      PC17        PC18      PC19      PC20
    ## 1   2.3008124 -4.408393  0.1388014 -1.987154   5.3570821 -1.168872 -3.302048
    ## 2   1.2110217  2.316774 -4.1315755  1.861779   0.9196818 -1.305010 -1.563090
    ## 3   0.2483641  2.691900 -3.2986713 -1.380452   0.9294441 -1.671877 -3.138875
    ## 4   1.5162442  2.048768 -1.5275178 -1.088731   0.7742403 -1.165909 -1.523211
    ## 5  10.3796784 -6.313460  4.6068477  7.631406 -11.6742066 -7.894054  8.533010
    ## 6 -10.1635818 -8.168188  5.2926374 -9.883659   5.7147461  2.318684 -2.175857
    ##         PC21      PC22       PC23      PC24       PC25      PC26       PC27
    ## 1  1.0821565  1.638081 -3.2145651 -4.758528  1.4895006  0.635360  2.0009133
    ## 2  0.3206081  4.569456 -2.6774190 -3.010654  1.7918365  3.590631 -2.1635115
    ## 3  3.3821109  1.941092  4.1961120 -4.021997 -5.8525362  4.754885  2.6657647
    ## 4 -2.5499057  2.176327  1.5534849 -1.531439 -1.9402259  4.341141 -0.4489018
    ## 5 -5.0932039 -5.987911  2.5989813  4.352300 -0.4705859 -0.934063  5.8693262
    ## 6  4.7052834 -8.561316  0.3865037  4.712929 -4.1930201  1.788674  1.0502503
    ##        PC28       PC29       PC30      PC31      PC32       PC33      PC34
    ## 1 -1.581707 -1.6399420 -0.7706624 -5.984501  0.248612  0.1507356  2.927915
    ## 2 -3.463865 -5.0093831 -5.1212250  2.743881 -3.877156  2.0096261 -3.026366
    ## 3 -1.276755  0.2187924  0.6582345  1.419515  6.519812 -4.2469713  2.050440
    ## 4 -3.382342  2.3276324  1.3256539 -1.157794  4.736545  1.4863124  2.419409
    ## 5  1.730407 -0.9446137 -3.9412674  2.308560 -1.895174  1.0568205 -3.367169
    ## 6  4.143549  0.1490191  0.3017374 -1.830570 -3.446601 -0.8154433 -3.491790
    ##        PC35       PC36      PC37       PC38       PC39       PC40       PC41
    ## 1 -2.841447 -2.2303863 -5.194586  3.6348250  0.8424234  2.8622780   2.628326
    ## 2  1.299273  1.1412407  3.168001  1.6242827 -3.0216367 -2.0602185 -10.803224
    ## 3  3.291428  0.9075436  2.726006 -3.3150204  1.5612708  0.5208121  -1.342924
    ## 4  3.551889  1.5450998 -4.227717 -1.1614388  2.0323782 -1.4266015   4.448948
    ## 5 -2.217989  1.0155489  4.151879 -0.2508915 -0.2837378 -0.7826875   1.481057
    ## 6 -2.603202  0.9258750 -1.139455  0.4752595 -1.8574172  0.4751768   1.001320
    ##          PC42        PC43       PC44       PC45       PC46        PC47
    ## 1 -2.48002257 -0.11045067 -1.7393470 -1.9248111 -0.3191336 -0.79911289
    ## 2 -0.58581460 -1.63028327 -0.4280295 -3.8478125  4.4589670 -3.00720747
    ## 3  0.52444876 -0.02466197 -1.9607043  3.4311289 -2.9119172  5.52716596
    ## 4  0.97885112  3.18849008 -1.7503717  1.8194832 -0.8685634 -3.87708721
    ## 5 -1.13736631 -0.25495689  1.0895862 -0.1080761 -0.2217693  0.87210373
    ## 6  0.09326974 -0.04500355 -0.2011824  0.6357723  0.5372292  0.08497015
    ##         PC48        PC49       PC50       PC51        PC52       PC53
    ## 1  5.7109798  6.94563318  0.6797816 -5.2518982  -0.8784367 -2.7488365
    ## 2 -0.2281665 -6.36445501  1.5811172 -4.0440641  -0.7461741 -1.6924277
    ## 3  0.5651592  0.02482585 -2.2604509  0.1560478   7.6746400 -0.1356940
    ## 4 -2.9680807  0.23359132  1.9370407  5.2416536 -11.2068653  3.7949851
    ## 5 -1.0266180  0.31933152  0.3478020  0.9699128   1.1310521  0.3889535
    ## 6 -1.5429469 -1.62680331 -0.3838217 -0.7817686   0.5375108 -0.2649840
    ##         PC54       PC55       PC56       PC57       PC58        PC59
    ## 1 -4.6071945 -0.3863459  6.3444800  3.5537213  1.1616091 -3.15389603
    ## 2  0.5273151 -1.9474896  0.2501325 -0.6721356 -0.3106257  3.30145895
    ## 3  6.1624506  6.2803873 -1.6187305 -0.9260990 -3.6330121 -0.68476680
    ## 4 -4.2723035 -2.2261186 -5.1536491 -2.6990065  3.0967847 -0.98438295
    ## 5 -0.1725670  0.9312331  1.4664179  0.8019946  0.4362105  0.04677094
    ## 6  1.2178922 -1.7909038 -0.5981598  0.1741735 -0.5310293  0.78874214
    ##          PC60       PC61        PC62       PC63       PC64        PC65
    ## 1 -0.72853948 -0.8179696 -4.00911177 -4.0693248  3.1375273  1.78583138
    ## 2 -2.55751349 -0.2214490 -1.58630140 -0.3162112 -1.3172589  0.01860059
    ## 3  2.84413535 -0.6312307 -0.98974231  4.1579054  1.1305623  0.16868990
    ## 4 -0.45794779  1.3005452  6.16275947  0.9210241 -0.8695925 -1.77768432
    ## 5  0.04059505  0.4436773 -0.33688242 -0.6215107  0.1731113  0.16299754
    ## 6 -0.56294925  0.4402412  0.07419497  0.2639030 -0.2702782  0.72367822
    ##         PC66       PC67       PC68        PC69        PC70       PC71
    ## 1  1.1401989 -1.4302502 -0.3384630 -0.49356834 -0.24770896  0.8230437
    ## 2 -0.7411071  1.5313028 -0.2300210  0.04206224  0.53493672 -0.3846649
    ## 3 -1.7165381 -1.1580465 -0.9399088 -0.55732425 -0.07863755  1.0636665
    ## 4 -0.2969664  0.4104425 -0.3779310  0.15588828 -1.02325043 -1.0341417
    ## 5  0.5819508  0.0898898  1.1664917  0.09933955 -0.19364041 -0.2169556
    ## 6  0.8385376 -0.3760198  0.5991345  0.85131628 -0.26280148 -0.3453531
    ##          PC72         PC73        PC74       PC75        PC76        PC77
    ## 1  0.39714058  0.009503675 -0.68748338 -0.4318030  0.43987525 -1.17577972
    ## 2  0.39616140 -0.721533502 -0.51230940  0.5794666 -0.54750768  0.91586622
    ## 3  1.15209964  0.222064924  2.13236342  0.5072941  0.77433367 -0.03963387
    ## 4 -1.54687391  0.841025150 -0.43746117 -0.4231308 -0.01147578  1.26776979
    ## 5  0.04976158 -0.728094297 -0.06822997 -0.1853361 -0.67847471  0.06380940
    ## 6 -0.53500790 -0.823281722 -0.25782010  0.1481072 -0.42785835 -0.32062580
    ##          PC78        PC79       PC80      PC81      PC82          PC83
    ## 1  1.09035756 -0.21041785  1.3088829 -1.577579 0.7368843  3.255586e-13
    ## 2 -0.44148482  0.12542092  1.0467837 -2.859669 0.6151647  1.641905e-14
    ## 3  0.65142556  0.26660946  0.7064822 -1.556606 0.8014818 -6.832850e-14
    ## 4 -1.34723374 -0.07896840  0.4075884 -1.661268 0.5492144  1.167646e-13
    ## 5 -0.09962754  0.16979451  0.4455882 -1.528856 0.4863900  1.214952e-13
    ## 6 -0.34260833  0.07530536 -0.2144310 -1.901886 0.2858714  4.219935e-15

``` r
fviz_eig(pca, geom = "bar", bar_width = 0.5, 
              addlabels = TRUE, hjust = -0.3,
              barfill = "#BF98A0", barcolor = "#BF98A0") + 
  ggtitle("PCA Variances of Training Data (Filtered, BMIQ Norm, Batch Corrected)") +
  theme_minimal() +
  labs(x = "Principal Components", y = "% of explained variances")
```

![](00_Training_files/figure-markdown_github/pca%20filt_norm_batchcorrected-2.png)

``` r
# rename columns
pca_scores <- pca$x %>% as_tibble() %>% mutate(Sample_Name = trainMeta$Sample_ID) 

# r squared
pca_cor_rs <- lmmatrix(dep = pca_scores[,1:10],
                        ind = as.data.frame(trainMeta) %>% 
                          dplyr::select(GEO_Dataset,
                                        Sex, 
                                        PE,
                                        Class,
                                        GA_RPC,
                                        Predicted_ethnicity_nothresh),
                       metric = 'Rsquared')

pca_plot_rs <- pca_cor_rs %>% as.data.frame() %>% 
  mutate(indep = rownames(pca_cor_rs)) %>%
 pivot_longer(c(-indep), names_to = "PC", values_to = "rsqr" ) %>%
  mutate(
  PC = factor(PC, levels = colnames(pca_cor_rs))) %>% as_tibble()

# Plot PCA 
p1_r <- ggplot(pca_plot_rs, aes(x = PC, y = indep, fill = rsqr)) +
  geom_tile(col = 'lightgrey') + theme_classic() +
  scale_x_discrete(expand = c(0, 0), labels = 1:20) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = '', fill = 'R Squared') +
  scale_fill_gradientn(colors = c("#ff4d00", "#f98010", "#f09623", "#e9b448")) +
  ggtitle("PCA - Training Data (Filtered, BMIQ Normalized, Batch Corrected)") 

p1_r
```

![](00_Training_files/figure-markdown_github/pca%20filt_norm_batchcorrected-3.png)

``` r
as.data.frame(pca_cor_rs) %>%
  filter(PC1 > 0.1) # no variables are no longer strongly associated to PC1 anymore after combat correction (all less than 0.5), PE and class are 0.21
```

    ##             PC1         PC2     PC3        PC4       PC5        PC6        PC7
    ## PE    0.2069356 0.009212334 0.10843 0.03751677 0.1276131 0.09986257 0.01554578
    ## Class 0.2069356 0.009212334 0.10843 0.03751677 0.1276131 0.09986257 0.01554578
    ##               PC8        PC9         PC10
    ## PE    0.002417302 0.02577868 1.682745e-06
    ## Class 0.002417302 0.02577868 1.682745e-06

``` r
# p value
pca_cor_pval <- lmmatrix(dep = pca_scores[,1:10],
                        ind = as.data.frame(trainMeta) %>% 
                          dplyr::select(GEO_Dataset,
                                        Sex, 
                                        PE,
                                        Class,
                                        GA_RPC,
                                        Predicted_ethnicity_nothresh),
                       metric = 'Pvalue')

pca_plot_pval <- pca_cor_pval %>% as.data.frame() %>% 
  mutate(indep = rownames(pca_cor_pval)) %>%
 pivot_longer(c(-indep), names_to = "PC", values_to = "pval" ) %>%
  mutate(
  PC = factor(PC, levels = colnames(pca_cor_pval))) %>% as_tibble()

# Plot PCA 
p1_r <- ggplot(pca_plot_pval, aes(x = PC, y = indep, fill = pval)) +
  geom_tile(col = 'lightgrey') + theme_classic() +
  scale_x_discrete(expand = c(0, 0), labels = 1:20) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = '', fill = 'P Value') +
  scale_fill_gradientn(colors = c("#0047AB", "#6495ED", "#ADD8E6")) +
  ggtitle("PCA - Training Data (Filtered, BMIQ Normalized, Batch Corrected)")

p1_r
```

![](00_Training_files/figure-markdown_github/pca%20filt_norm_batchcorrected-4.png)

``` r
as.data.frame(pca_cor_pval) %>%
  filter(PC1 < 0.05) # PE and Class, but not Dataset.
```

    ##                                       PC1       PC2         PC3        PC4
    ## PE                           1.556865e-05 0.3880486 0.002368086 0.07934194
    ## Class                        1.556865e-05 0.3880486 0.002368086 0.07934194
    ## Predicted_ethnicity_nothresh 3.369244e-02 0.4031463 0.087432410 0.85474352
    ##                                       PC5         PC6         PC7          PC8
    ## PE                           0.0009155405 0.003610802 0.261405331 6.589247e-01
    ## Class                        0.0009155405 0.003610802 0.261405331 6.589247e-01
    ## Predicted_ethnicity_nothresh 0.1355518164 0.307908018 0.002983034 4.612566e-05
    ##                                       PC9        PC10
    ## PE                           1.470615e-01 0.990713721
    ## Class                        1.470615e-01 0.990713721
    ## Predicted_ethnicity_nothresh 3.428688e-09 0.001611903

# Save

``` r
saveRDS(trainComBat, here::here("Code", "2022_GitHub", "01_Data", "train.rds"))
write.csv(trainMeta, here::here("Code", "2022_GitHub", "01_Data", "trainMeta.csv"))
```
