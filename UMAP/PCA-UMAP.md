PCA-UMAP - D. sechellia
================

In this script I will try to detect any signal of population structure
in the sequences of the interesting TEs. To filter the SNPs, I used the
script *mpileup2PCA.py*.

## Filtering strategy

For a position to be selected as “SNP”, it must have at least 2 alleles
which are present in at least *n* samples with a frequency above *m*. In
the script, *n* and *m* are the two parameters to be specified as
arguments: **min-count** and **min-freq**.

The idea behind is to select only alleles which are present in a
relevant fraction of TE insertions in an individual (frequency
threshold), and these insertions with these alleles must be shared among
at least few individuals (to select population specific alleles).

## Metadata

For detecting population structure, we focus on **GDL** samples, cause
old samples without TEs are not interesting in this context and would
only hide potential structure within modern populations.

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.1     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggpubr)
library(patchwork)
library(umap)

theme_set(theme_bw())

sech_more_meta <- read_tsv("/Volumes/Temp1/simulans-old-strains/other-data/sechellia-teisseri/sechellia-moredata-metadata.txt", show_col_types = FALSE) %>% rename(ID = "run_accession") %>% mutate(
    location = case_when(
      startsWith(library_name, "Anro") ~ "Anro",
      startsWith(library_name, "Denis") ~ "Denis",
      startsWith(library_name, "LD") ~ "LD",
      startsWith(library_name, "maria_") ~ "Maria",
      startsWith(library_name, "mariane") ~ "Mariane",
      startsWith(library_name, "PNF") ~ "PNF")) #%>% filter(location %in% c("PNF", "Denis"))
```

## PCA and UMAP functions

``` r
PCA <- function(af, metadata, titlee){
  
  (full_dataset <- inner_join(metadata, af, by="ID") %>% distinct() %>% type_convert())
  pcaable <- full_dataset %>% select_if(~ !all(. == .[1]))
  pca_result <- prcomp(pcaable[, -c(1:3)], center = TRUE, scale = TRUE)
  var_explained <- pca_result$sdev^2/sum(pca_result$sdev^2)
  
  plot1 <- ggplot(data.frame(pca_result$x, ID=full_dataset$ID, location=full_dataset$location), aes(x=PC1,y=PC2, color=location)) + geom_point(size=2) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%"), color="Location") + ggtitle(titlee) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_text(face = "bold"))
  
  plot2 <- ggplot(data.frame(pca_result$x, ID=full_dataset$ID, location=full_dataset$location), aes(x=PC1,y=PC3, color=location)) + geom_point(size=2) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC3: ",round(var_explained[3]*100,1),"%"), color="Location") + ggtitle(titlee) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_text(face = "bold"))
  
  plot3 <- ggplot(data.frame(pca_result$x, ID=full_dataset$ID, location=full_dataset$location), aes(x=PC3,y=PC2, color=location)) + geom_point(size=2) + labs(x=paste0("PC3: ",round(var_explained[3]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%"), color="Location") + ggtitle(titlee) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_text(face = "bold"))
  
  plot4 <- ggplot(data.frame(pca_result$x, ID=full_dataset$ID, location=full_dataset$location), aes(x=PC1,y=PC4, color=location)) + geom_point(size=2) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC4: ",round(var_explained[4]*100,1),"%"), color="Location") + ggtitle(titlee) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_text(face = "bold"))
  
  plot5 <- ggplot(data.frame(pca_result$x, ID=full_dataset$ID, location=full_dataset$location), aes(x=PC4,y=PC2, color=location)) + geom_point(size=2) + labs(x=paste0("PC4: ",round(var_explained[4]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%"), color="Location") + ggtitle(titlee) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_text(face = "bold"))
  
  plot6 <- ggplot(data.frame(pca_result$x, ID=full_dataset$ID, location=full_dataset$location), aes(x=PC3,y=PC4, color=location)) + geom_point(size=2) + labs(x=paste0("PC3: ",round(var_explained[3]*100,1),"%"), y=paste0("PC4: ",round(var_explained[4]*100,1),"%"), color="Location") + ggtitle(titlee) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_text(face = "bold"))

  list(pc1_2 = plot1, pc1_3 = plot2, pc2_3 = plot3, pc1_4 = plot4, pc2_4 = plot5, pc3_4 = plot6)
}

UMAP <- function(af, metadata, titlee){
  
  full_dataset <- inner_join(metadata, af, by="ID") %>% distinct() %>% type_convert()
  pcaable <- full_dataset %>% select_if(~ !all(. == .[1]))
  umappable <- as.matrix(pcaable[, -c(1:3)])
  umap_result <- umap(umappable, n_neighbors = 10, min_dist = 0.3)
  umap <- umap_result$layout %>% as.data.frame() %>% rename(UMAP1="V1",UMAP2="V2")
  
  plot <- umap %>% ggplot(aes(x = UMAP1, y = UMAP2, color = full_dataset$location)) +
  geom_point(size=2)+ labs(x = "UMAP1", y = "UMAP2", title = titlee, color = "Location") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_text(face = "bold"))
  plot
}
```

## Dsech

Parameters used:

- min-freq 0.2
- min-count 2

``` r
shellder <- read_tsv("/Volumes/Temp1/simulans-old-strains/PCA-UMAP/sechellia/Dsec-shellder.PCAable") #%>%
```

    ## Rows: 46 Columns: 18
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): ID
    ## dbl (17): gypsy-29-dsim_442, gypsy-29-dsim_452, gypsy-29-dsim_625, gypsy-29-...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
  #filter(!(ID %in% c("SRR5860639.cleaned", "SRR5860640.cleaned", "SRR5860626.cleaned")))
shellder$ID <- str_remove(shellder$ID, "\\.cleaned$")

(umap_shellder <- UMAP(shellder, sech_more_meta, "shellder"))
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   ID = col_character(),
    ##   library_name = col_character(),
    ##   location = col_character()
    ## )

![](PCA-UMAP_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
pca_shellder <- PCA(shellder, sech_more_meta, "shellder")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   ID = col_character(),
    ##   library_name = col_character(),
    ##   location = col_character()
    ## )

``` r
#pca_shellder$pc1_2
#pca_shellder$pc1_3
#pca_shellder$pc2_3
#pca_shellder$pc1_4
#pca_shellder$pc2_4
#pca_shellder$pc3_4

#ggsave("/Volumes/Temp1/Dmel-stealthTEs/PCA/OPUS/pca.png", pca_OPUS, dpi = 300, width = 8, height = 6)
#ggsave("/Volumes/Temp1/Dmel-stealthTEs/PCA/OPUS/umap.png", umap_OPUS, dpi = 300, width = 8, height = 6)
```

``` r
spoink <- read_tsv("/Volumes/Temp1/simulans-old-strains/PCA-UMAP/sechellia/Dsec-spoink.PCAable") #%>%
```

    ## Rows: 46 Columns: 283
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): ID
    ## dbl (282): gypsy-7-sim1_141, gypsy-7-sim1_333, gypsy-7-sim1_415, gypsy-7-sim...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
  #filter(!(ID %in% c("SRR5860639.cleaned", "SRR5860640.cleaned", "SRR5860626.cleaned")))
spoink$ID <- str_remove(spoink$ID, "\\.cleaned$")

(umap_spoink <- UMAP(spoink, sech_more_meta, "spoink"))
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   ID = col_character(),
    ##   library_name = col_character(),
    ##   location = col_character()
    ## )

![](PCA-UMAP_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
spoink_present <- spoink %>% filter(ID %in% c("SRR5860583","SRR5860625","SRR5860627","SRR5860628","SRR5860629","SRR5860630","SRR5860631","SRR5860632","SRR5860633","SRR5860634","SRR5860643","SRR5860644","SRR5860666","SRR5860667","SRR5860668","SRR5860669","SRR5860670","SRR5860671","SRR5860672","SRR5860673","SRR5860674","SRR5860675"))

pca_spoink <- PCA(spoink_present, sech_more_meta, "spoink")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   ID = col_character(),
    ##   library_name = col_character(),
    ##   location = col_character()
    ## )

``` r
pca_spoink$pc1_2
```

![](PCA-UMAP_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
(umap_spoink <- UMAP(spoink_present, sech_more_meta, "spoink"))
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   ID = col_character(),
    ##   library_name = col_character(),
    ##   location = col_character()
    ## )

![](PCA-UMAP_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
#pca_shellder$pc1_3
#pca_shellder$pc2_3
#pca_shellder$pc1_4
#pca_shellder$pc2_4
#pca_shellder$pc3_4

#ggsave("/Volumes/Temp1/Dmel-stealthTEs/PCA/OPUS/pca.png", pca_OPUS, dpi = 300, width = 8, height = 6)
#ggsave("/Volumes/Temp1/Dmel-stealthTEs/PCA/OPUS/umap.png", umap_OPUS, dpi = 300, width = 8, height = 6)
```

## Private SNPs functions

Function which finds private SNPs in each population. Takes as input a
merged file (PCAable + metadata) and a frequency threshold (es. 0.5). It
counts how many SNPs are present in at least 50% of the individuals in a
population and not present in more than 15 other samples.

``` r
private <- function(data, threshold, percentage_pop_with_snp){
snps <- colnames(data)[3:length(colnames(data))]
populations <- data %>% select(location) %>% distinct() %>% pull()
private_snps <- tibble(location=NULL, SNP=NULL, thr=NULL, perc_shared_in_pop=NULL, number_external=NULL)

for (pop in populations){
  for (snp in snps){
    p <- data %>% select(ID, location, !!snp) %>% filter(location==pop)
    other <- filter(data, location!=pop)
    shared <- sum(p[[snp]] < threshold)
    count_pop <- p %>% summarise(count = n()) %>% pull()
    other_shared <- sum(other[[snp]] < threshold)
    if ((shared > (count_pop/100*percentage_pop_with_snp)) & other_shared<3){
      row <- tibble(location=pop, SNP=snp, thr=threshold, perc_shared_in_pop=shared/count_pop, number_external=other_shared)
      private_snps <- bind_rows(private_snps, row)
    }
  }
}
private_snps
}
```

Function to loop the function “private()” over multiple thresholds (from
0.5 to 0.9). A private SNP with 0.5 frequency is a very strong proof of
a bottleneck during the invasion, while a SNP with only 0.9 is not that
significant.

``` r
different_thresholds <- function(data, percentage_pop_with_snp){
result_tibble <- tibble(location=NULL, SNP=NULL, thr=NULL, perc_shared_in_pop=NULL, number_external=NULL)
for (t in seq(0.5, 0.9, +0.1)){
  tib <- private(data, t, percentage_pop_with_snp)
  result_tibble <- bind_rows(result_tibble, tib)
}
pops <- result_tibble %>% select(location) %>% distinct() %>% pull()
filtered_tibble <- tibble(location=NULL, SNP=NULL, thr=NULL, perc_shared_in_pop=NULL, number_external=NULL)
for (pop in pops){
  only_pop <- filter(result_tibble, location==pop)
  filtered <- only_pop %>% filter(thr==0.5)
for (t in seq(0.6, 0.9, +0.1)){
  th <- filter(only_pop, thr<t) %>% select(SNP) %>% pull()
  non_overlapping <- filter(only_pop, thr==t, !(SNP %in% th))
  filtered <- bind_rows(filtered, non_overlapping)
}
  filtered_tibble <- bind_rows(filtered_tibble, filtered)
}
filtered_tibble
}
```

Function which takes as input the tibble returned from
“different_thresholds” and creates a barplot.

``` r
private_plot <- function(data, titlee){
  plottable <- data %>% group_by(location, thr) %>% summarise(count = n()) %>% mutate(thr=1-thr)
  ggplot(plottable, aes(x = thr, y = count, fill = location)) +
  geom_bar(stat = "identity") +
  labs(x = "Frequency threshold", y = "Private SNPs", fill = "Location") +
  ggtitle(titlee) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ylim(0,25) #+ xlim(0,0.6)
}
```

``` r
meta_shellder <- inner_join(sech_more_meta, shellder, by="ID") %>% filter(location %in% c("Denis","PNF")) %>% filter(!(ID %in% c("SRR5860639", "SRR5860640", "SRR5860626")))
(tib_shellder <- different_thresholds(meta_shellder, 50))
```

    ## # A tibble: 7 × 5
    ##   location SNP                  thr perc_shared_in_pop number_external
    ##   <chr>    <chr>              <dbl>              <dbl>           <int>
    ## 1 PNF      gypsy-29-dsim_4380   0.7              1                   0
    ## 2 PNF      gypsy-29-dsim_3025   0.8              0.571               0
    ## 3 PNF      gypsy-29-dsim_5220   0.8              0.857               0
    ## 4 PNF      gypsy-29-dsim_442    0.9              0.857               0
    ## 5 Denis    gypsy-29-dsim_452    0.8              0.867               1
    ## 6 Denis    gypsy-29-dsim_4488   0.9              0.6                 0
    ## 7 Denis    gypsy-29-dsim_6037   0.9              0.6                 0

``` r
(plot_shellder <- private_plot(tib_shellder, "shellder"))
```

    ## `summarise()` has grouped output by 'location'. You can override using the
    ## `.groups` argument.

![](PCA-UMAP_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
meta_spoink <- inner_join(sech_more_meta, spoink, by="ID") %>% filter(location %in% c("Denis","PNF")) %>% filter(!(ID %in% c("SRR5860639", "SRR5860640", "SRR5860626")))
(tib_spoink <- different_thresholds(meta_spoink, 50))
```

    ## # A tibble: 59 × 5
    ##    location SNP                 thr perc_shared_in_pop number_external
    ##    <chr>    <chr>             <dbl>              <dbl>           <int>
    ##  1 PNF      gypsy-7-sim1_2640   0.5              1                   0
    ##  2 PNF      gypsy-7-sim1_4414   0.5              0.857               0
    ##  3 PNF      gypsy-7-sim1_4629   0.5              0.857               1
    ##  4 PNF      gypsy-7-sim1_4333   0.9              0.714               1
    ##  5 PNF      gypsy-7-sim1_4440   0.9              1                   0
    ##  6 Denis    gypsy-7-sim1_415    0.7              0.8                 0
    ##  7 Denis    gypsy-7-sim1_2897   0.8              0.8                 0
    ##  8 Denis    gypsy-7-sim1_3731   0.8              0.8                 0
    ##  9 Denis    gypsy-7-sim1_333    0.9              0.867               0
    ## 10 Denis    gypsy-7-sim1_1556   0.9              0.667               0
    ## # ℹ 49 more rows

``` r
(plot_spoink <- private_plot(tib_spoink, "spoink"))
```

    ## `summarise()` has grouped output by 'location'. You can override using the
    ## `.groups` argument.

    ## Warning: Removed 1 rows containing missing values (`position_stack()`).

![](PCA-UMAP_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->
