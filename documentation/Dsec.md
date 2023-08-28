Drosophila sechellia - Copynumber analysis
================

``` r
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ggpubr))
theme_set(theme_bw())

knitr::opts_knit$set(root.dir = "/Volumes/Temp1/simulans-old-strains/analysis/plots")
```

``` r
analysis <- function(csv, old, new, year, titlee) {
  
copynumbers <- read_csv(csv, show_col_types = FALSE) %>% filter(Sample!="Sample", Sample%in%c(old, new)) %>% type_convert() %>% mutate(age = ifelse(Sample %in% old, "old", "new")) %>% inner_join(year, by="Sample") %>% arrange(collection_year) %>% mutate(Sample = paste0(Sample, "_", collection_year)) %>% select(-collection_year)

avg_museum_modern <- copynumbers %>% group_by(age, TE) %>% summarise(All_reads=mean(All_reads))
museum <- avg_museum_modern %>% filter(age == "old") %>% ungroup() %>% select(TE, All_reads) %>% rename(old_cn = All_reads)
modern <- avg_museum_modern %>% filter(age == "new") %>% ungroup() %>% select(TE, All_reads) %>% rename(new_cn = All_reads)
fold <- inner_join(museum, modern, by="TE") %>% type_convert() %>% mutate(fold_enrichment = new_cn/old_cn) %>% filter(new_cn > 2)

plot_fold <- ggplot(data = fold, aes(x = TE, y = fold_enrichment)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = fold_enrichment > 2)) +
  scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  xlab("TE") +
  ylab("fold enrichment") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) +
  geom_text(aes(label = ifelse(fold_enrichment > 2, as.character(TE), "")), 
            position = position_dodge(width = 1),
            vjust = -0.5, size = 3, angle=90)+
  ggtitle(titlee)

fold_enriched <- fold %>% filter(fold_enrichment>2) %>% arrange(desc(fold_enrichment))

comparison_foldenriched <- copynumbers %>% ungroup() %>% select(-HQ_reads, -age) %>% pivot_wider(names_from = Sample, values_from = All_reads) %>% filter(TE %in% fold_enriched$TE)# %>% inner_join(fold_enriched, by="fold_enrichment")

list(plot = plot_fold, table = fold_enriched, each_sample = comparison_foldenriched)
}
```

``` r
timeline <- function(csv, age, transposon) {
  single_te <- read_csv(csv, show_col_types = FALSE) %>% filter(Sample!="Sample") %>% type_convert() %>% filter(TE==transposon) %>% inner_join(age, by="Sample") %>% arrange(collection_year)
  
  plot <- ggplot(single_te, aes(x = reorder(Sample, collection_year), y = HQ_reads)) +
    geom_point(size = 1.5) +
    scale_x_discrete(labels = paste0(single_te$collection_year)) +
    labs(y = "copynumber", x = "collection_year") + ggtitle(transposon)
}
```

``` r
samples <- c("SRR14138507", "SRR14138506", "SRR9913024", "SRR5860661", "SRR5860660", "SRR5860658") 
years <- c(2012, 2012, 1980, 2012, 2012, 2012)
age_dsec <- tibble(collection_year = years, Sample = samples)

old <- c("SRR9913024")
new <- c("SRR14138507", "SRR14138506","SRR5860661", "SRR5860660", "SRR5860658")

islands <- c("Anro", "Denis", "Cousin Island", "LD", "Anro", "Anro")
```

## Dmel library (+ Gypsy7 and Gypsy29 from Dsim)

``` r
dmel_TEs <- analysis("/Volumes/Temp1/simulans-old-strains/analysis/csv/Dsec/dmel_TEs/Dsec-dmel_TEs.csv", old, new, age_dsec, "Old vs New")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Sample = col_character(),
    ##   TE = col_character(),
    ##   All_reads = col_double(),
    ##   HQ_reads = col_double()
    ## )

    ## `summarise()` has grouped output by 'age'. You can override using the `.groups`
    ## argument.
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── cols( TE =
    ## col_character() )

``` r
dmel_TEs$plot
```

![](Dsec_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
kable(dmel_TEs$table)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
TE
</th>
<th style="text-align:right;">
old_cn
</th>
<th style="text-align:right;">
new_cn
</th>
<th style="text-align:right;">
fold_enrichment
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
KEPLER
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
2.494
</td>
<td style="text-align:right;">
9.237037
</td>
</tr>
</tbody>
</table>

``` r
kable(dmel_TEs$each_sample)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
TE
</th>
<th style="text-align:right;">
SRR9913024_1980
</th>
<th style="text-align:right;">
SRR14138506_2012
</th>
<th style="text-align:right;">
SRR14138507_2012
</th>
<th style="text-align:right;">
SRR5860658_2012
</th>
<th style="text-align:right;">
SRR5860660_2012
</th>
<th style="text-align:right;">
SRR5860661_2012
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
KEPLER
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
4.19
</td>
<td style="text-align:right;">
4.17
</td>
<td style="text-align:right;">
1.24
</td>
<td style="text-align:right;">
1.59
</td>
<td style="text-align:right;">
1.28
</td>
</tr>
</tbody>
</table>

``` r
dmel_TEs <- analysis("/Volumes/Temp1/simulans-old-strains/analysis/csv/Dsec/dmel_TEs/Dsec-dmel_TEs.csv", new, old, age_dsec, "New vs Old")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Sample = col_character(),
    ##   TE = col_character(),
    ##   All_reads = col_double(),
    ##   HQ_reads = col_double()
    ## )
    ## 
    ## `summarise()` has grouped output by 'age'. You can override using the `.groups` argument.
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   TE = col_character()
    ## )

``` r
dmel_TEs$plot
```

![](Dsec_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
kable(dmel_TEs$table)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
TE
</th>
<th style="text-align:right;">
old_cn
</th>
<th style="text-align:right;">
new_cn
</th>
<th style="text-align:right;">
fold_enrichment
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
gypsy-29-dsim
</td>
<td style="text-align:right;">
0.032
</td>
<td style="text-align:right;">
6.02
</td>
<td style="text-align:right;">
188.125000
</td>
</tr>
<tr>
<td style="text-align:left;">
gypsy-7-dsim
</td>
<td style="text-align:right;">
3.416
</td>
<td style="text-align:right;">
15.95
</td>
<td style="text-align:right;">
4.669204
</td>
</tr>
</tbody>
</table>

``` r
kable(dmel_TEs$each_sample)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
TE
</th>
<th style="text-align:right;">
SRR9913024_1980
</th>
<th style="text-align:right;">
SRR14138506_2012
</th>
<th style="text-align:right;">
SRR14138507_2012
</th>
<th style="text-align:right;">
SRR5860658_2012
</th>
<th style="text-align:right;">
SRR5860660_2012
</th>
<th style="text-align:right;">
SRR5860661_2012
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
gypsy-29-dsim
</td>
<td style="text-align:right;">
6.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
gypsy-7-dsim
</td>
<td style="text-align:right;">
15.95
</td>
<td style="text-align:right;">
5.03
</td>
<td style="text-align:right;">
9.39
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.91
</td>
</tr>
</tbody>
</table>

``` r
candidates <- c("gypsy-29-dsim", "gypsy-7-dsim")

for (t in candidates){
 dmel_TEs_timeline <- timeline("/Volumes/Temp1/simulans-old-strains/analysis/csv/Dsec/dmel_TEs/Dsec-dmel_TEs.csv", age_dsec, t)
  print(dmel_TEs_timeline)
}
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Sample = col_character(),
    ##   TE = col_character(),
    ##   All_reads = col_double(),
    ##   HQ_reads = col_double()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Sample = col_character(),
    ##   TE = col_character(),
    ##   All_reads = col_double(),
    ##   HQ_reads = col_double()
    ## )

![](Dsec_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](Dsec_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

## Clark library

## Dmel library (+ Gypsy7 and Gypsy29 from Dsim)

``` r
clark <- analysis("/Volumes/Temp1/simulans-old-strains/analysis/csv/Dsec/clark/Dsec-clark.csv", old, new, age_dsec, "Old vs New")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Sample = col_character(),
    ##   TE = col_character(),
    ##   All_reads = col_double(),
    ##   HQ_reads = col_double()
    ## )

    ## `summarise()` has grouped output by 'age'. You can override using the `.groups`
    ## argument.
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── cols( TE =
    ## col_character() )

``` r
clark$plot
```

![](Dsec_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
kable(clark$table)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
TE
</th>
<th style="text-align:right;">
old_cn
</th>
<th style="text-align:right;">
new_cn
</th>
<th style="text-align:right;">
fold_enrichment
</th>
</tr>
</thead>
<tbody>
<tr>
</tr>
</tbody>
</table>

``` r
kable(clark$each_sample)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
TE
</th>
<th style="text-align:right;">
SRR9913024_1980
</th>
<th style="text-align:right;">
SRR14138506_2012
</th>
<th style="text-align:right;">
SRR14138507_2012
</th>
<th style="text-align:right;">
SRR5860658_2012
</th>
<th style="text-align:right;">
SRR5860660_2012
</th>
<th style="text-align:right;">
SRR5860661_2012
</th>
</tr>
</thead>
<tbody>
<tr>
</tr>
</tbody>
</table>

``` r
clark <- analysis("/Volumes/Temp1/simulans-old-strains/analysis/csv/Dsec/clark/Dsec-clark.csv", new, old, age_dsec, "New vs Old")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Sample = col_character(),
    ##   TE = col_character(),
    ##   All_reads = col_double(),
    ##   HQ_reads = col_double()
    ## )
    ## 
    ## `summarise()` has grouped output by 'age'. You can override using the `.groups` argument.
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   TE = col_character()
    ## )

``` r
clark$plot
```

![](Dsec_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
kable(clark$table)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
TE
</th>
<th style="text-align:right;">
old_cn
</th>
<th style="text-align:right;">
new_cn
</th>
<th style="text-align:right;">
fold_enrichment
</th>
</tr>
</thead>
<tbody>
<tr>
</tr>
</tbody>
</table>

``` r
kable(clark$each_sample)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
TE
</th>
<th style="text-align:right;">
SRR9913024_1980
</th>
<th style="text-align:right;">
SRR14138506_2012
</th>
<th style="text-align:right;">
SRR14138507_2012
</th>
<th style="text-align:right;">
SRR5860658_2012
</th>
<th style="text-align:right;">
SRR5860660_2012
</th>
<th style="text-align:right;">
SRR5860661_2012
</th>
</tr>
</thead>
<tbody>
<tr>
</tr>
</tbody>
</table>
