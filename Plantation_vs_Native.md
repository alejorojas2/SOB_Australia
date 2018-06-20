---
title: "Plantation vs Native data anlaysis"
author: "Alejandro Rojas"
date: "June 5, 2018"
output:
  html_document: 
    keep_md: TRUE
---


```r
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(vegan)
library(readr)
library(ampvis)

load(file = "SOB_files.rda")
```

## Checking sequencing depth


```r
# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(SOB_data))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

## Standardizing by sequencing depth


```r
#Standardize abundances to the median sequencing depth
total <- median(sample_sums(SOB_data))
standf <- function(x, t=total) round(t * (x/sum(x)))
SOB_data.std <- transform_sample_counts(SOB_data, standf)

#Filter taxa with cutoff 3.0 Coefficient of Variation
SOB_data.stdf <- filter_taxa(SOB_data.std, function(x) sd(x)/mean(x) > 3.0, TRUE)
```

## Bar plots after standardizing


```r
fungi.p <- subset_taxa(SOB_data.stdf, Phylum!="p__unidentified")
fungi.p <- tax_glom(fungi.p, taxrank="Phylum")
fungi.p2 <- subset_samples(fungi.p, site_code != "null")

code_names <- c(INV = "Invasive",
                PL = "Plantation",
                UN = "Native")

bar_phyla_code <- plot_bar(fungi.p2, x = "site_name", fill = "Phylum") 
bar_phyla_code + geom_bar(stat = "identity", color="black", size=0.2,position = "stack") + 
  facet_grid(site_code ~ ., labeller = labeller(site_code = code_names)) +
  scale_fill_brewer(type = "div", palette = "Paired")
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

## Ordination

### Original plot


```r
#Transforming to proportions
SOB_data.prop <- transform_sample_counts(SOB_data.stdf, function(otu) otu/sum(otu))
SOB_data.prop1 <- subset_samples(SOB_data.prop, site_code != "null")

#Palette
library(RColorBrewer)
pal <- colorRampPalette(brewer.pal(12, "Paired"))

#Ordination
ord.SOB <- ordinate(SOB_data.prop1, "NMDS", "bray")
```

```
## Run 0 stress 0.3272439 
## Run 1 stress 0.324991 
## ... New best solution
## ... Procrustes: rmse 0.02641423  max resid 0.2188735 
## Run 2 stress 0.3293062 
## Run 3 stress 0.3243587 
## ... New best solution
## ... Procrustes: rmse 0.02582157  max resid 0.2191589 
## Run 4 stress 0.3236486 
## ... New best solution
## ... Procrustes: rmse 0.02284304  max resid 0.2224103 
## Run 5 stress 0.3253105 
## Run 6 stress 0.3243008 
## Run 7 stress 0.3273157 
## Run 8 stress 0.3272866 
## Run 9 stress 0.3275019 
## Run 10 stress 0.323585 
## ... New best solution
## ... Procrustes: rmse 0.01134536  max resid 0.1814109 
## Run 11 stress 0.3250794 
## Run 12 stress 0.3368188 
## Run 13 stress 0.3235016 
## ... New best solution
## ... Procrustes: rmse 0.023272  max resid 0.2216099 
## Run 14 stress 0.3245983 
## Run 15 stress 0.3234724 
## ... New best solution
## ... Procrustes: rmse 0.01872211  max resid 0.2254138 
## Run 16 stress 0.325958 
## Run 17 stress 0.3308166 
## Run 18 stress 0.3358738 
## Run 19 stress 0.4182711 
## Run 20 stress 0.3250799 
## *** No convergence -- monoMDS stopping criteria:
##      3: no. of iterations >= maxit
##     17: stress ratio > sratmax
```

```r
plot_ordination(SOB_data.prop1, ord.SOB, shape = "site_code", color = "site_name") + 
  geom_point(size = 3) + scale_color_manual(values = pal(42)) + labs(title = "All samples indicated by site name and code")
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


### Focusing on Plantations



```r
#Filtering data
SOB_data.prop.PL <- subset_samples(SOB_data.prop, site_code == "PL")

#Ordination
ord.SOB.PL <- ordinate(SOB_data.prop.PL, "NMDS", "bray")
```

```
## Run 0 stress 0.2891296 
## Run 1 stress 0.2918815 
## Run 2 stress 0.2895129 
## ... Procrustes: rmse 0.01247778  max resid 0.113833 
## Run 3 stress 0.2899605 
## Run 4 stress 0.2907021 
## Run 5 stress 0.2895129 
## ... Procrustes: rmse 0.01229969  max resid 0.1138403 
## Run 6 stress 0.292045 
## Run 7 stress 0.2890911 
## ... New best solution
## ... Procrustes: rmse 0.02137696  max resid 0.1750331 
## Run 8 stress 0.2880312 
## ... New best solution
## ... Procrustes: rmse 0.01740416  max resid 0.1754765 
## Run 9 stress 0.288901 
## Run 10 stress 0.291168 
## Run 11 stress 0.2891288 
## Run 12 stress 0.2900775 
## Run 13 stress 0.2889627 
## Run 14 stress 0.3177512 
## Run 15 stress 0.2891363 
## Run 16 stress 0.2901402 
## Run 17 stress 0.2894934 
## Run 18 stress 0.2881973 
## ... Procrustes: rmse 0.005919558  max resid 0.06167274 
## Run 19 stress 0.2877326 
## ... New best solution
## ... Procrustes: rmse 0.01087331  max resid 0.1147065 
## Run 20 stress 0.2880902 
## ... Procrustes: rmse 0.01030379  max resid 0.1138105 
## *** No convergence -- monoMDS stopping criteria:
##     20: stress ratio > sratmax
```

```r
plot_ordination(SOB_data.prop.PL, ord.SOB.PL, color = "state") + 
  geom_point(size = 3) + scale_color_brewer(type = "div", palette = "Set1",
                                            name = "State",
                                            labels = c("Australian Capital Territory",
                                                       "New South Wales",
                                                       "South Australia",
                                                       "Victoria",
                                                       "Western Australia")) +
  labs(title = "Plantation sites colored by state")
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

### Focusing on Native sites including Nullabor



```r
#Filtering data
SOB_data.prop.UN <- subset_samples(SOB_data.prop, 
                                   site_code == "UN" |
                                    site_code == "null")

#Ordination
ord.SOB.UN <- ordinate(SOB_data.prop.UN, "NMDS", "bray")
```

```
## Run 0 stress 0.295749 
## Run 1 stress 0.2960179 
## ... Procrustes: rmse 0.05379557  max resid 0.1848509 
## Run 2 stress 0.2894091 
## ... New best solution
## ... Procrustes: rmse 0.0649258  max resid 0.2105987 
## Run 3 stress 0.2929966 
## Run 4 stress 0.2900297 
## Run 5 stress 0.291876 
## Run 6 stress 0.2920973 
## Run 7 stress 0.2882637 
## ... New best solution
## ... Procrustes: rmse 0.02047467  max resid 0.09910336 
## Run 8 stress 0.2900333 
## Run 9 stress 0.2957215 
## Run 10 stress 0.3004788 
## Run 11 stress 0.2898515 
## Run 12 stress 0.3051156 
## Run 13 stress 0.2903868 
## Run 14 stress 0.2897819 
## Run 15 stress 0.3011159 
## Run 16 stress 0.2899798 
## Run 17 stress 0.2896673 
## Run 18 stress 0.290406 
## Run 19 stress 0.2898804 
## Run 20 stress 0.2889866 
## *** No convergence -- monoMDS stopping criteria:
##     20: stress ratio > sratmax
```

```r
plot_ordination(SOB_data.prop.UN, ord.SOB.UN, color = "state", shape = "site_code") + 
  geom_point(size = 3) + scale_color_brewer(type = "div", palette = "Set1",
                                            name = "State",
                                            labels = c("Australian Capital Territory",
                                                       "New South Wales",
                                                       "South Australia",
                                                       "Victoria",
                                                       "Western Australia")) +
  scale_shape_discrete(labels = c("Nullabor", "Native"), name = "Site") +
  labs(title = "Native sites, including Nullabor, colored by state")
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-6-1.png)<!-- -->


### Focusing on __only Native sites__



```r
#Filtering data
SOB_data.prop.UN1 <- subset_samples(SOB_data.prop, site_code == "UN")

#Ordination
ord.SOB.UN1 <- ordinate(SOB_data.prop.UN1, "NMDS", "bray")
```

```
## Run 0 stress 0.3250856 
## Run 1 stress 0.3153795 
## ... New best solution
## ... Procrustes: rmse 0.06878835  max resid 0.2661971 
## Run 2 stress 0.3267914 
## Run 3 stress 0.3265395 
## Run 4 stress 0.3293899 
## Run 5 stress 0.3340861 
## Run 6 stress 0.3255799 
## Run 7 stress 0.3145817 
## ... New best solution
## ... Procrustes: rmse 0.03474687  max resid 0.3123634 
## Run 8 stress 0.3150517 
## ... Procrustes: rmse 0.03472002  max resid 0.3151576 
## Run 9 stress 0.31541 
## Run 10 stress 0.315621 
## Run 11 stress 0.3202786 
## Run 12 stress 0.315084 
## Run 13 stress 0.3249198 
## Run 14 stress 0.314252 
## ... New best solution
## ... Procrustes: rmse 0.02134826  max resid 0.09567648 
## Run 15 stress 0.3222605 
## Run 16 stress 0.3149147 
## Run 17 stress 0.3154573 
## Run 18 stress 0.3214494 
## Run 19 stress 0.3289438 
## Run 20 stress 0.315335 
## *** No convergence -- monoMDS stopping criteria:
##      1: no. of iterations >= maxit
##     19: stress ratio > sratmax
```

```r
plot_ordination(SOB_data.prop.UN1, ord.SOB.UN1, color = "state") +
  geom_point(size = 3) + 
  scale_color_brewer(type = "div", palette = "Set1", name = "State",
                     labels = c("Australian Capital Territory",
                                "New South Wales",
                                "South Australia",
                                "Victoria",
                                "Western Australia")) + 
  scale_shape_discrete(labels = c("Nullabor", "Native"), name = "Site") +
  labs(title = "Native sites colored by state")
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

## Heatmap

### Heatmap of plantation and native sites



```r
#Filtering data only to plantation and native sites
SOB_data.stdf.1 <- subset_samples(SOB_data.stdf, 
                                   site_code == "UN" |
                                    site_code == "PL")

#Heatmap
amp_heatmap(data = SOB_data.stdf.1,
            group = c("state", "site_code"),
            tax.show = 50,
            scale.seq = 100,
            plot.text.size = 2,
            tax.aggregate = "Genus",
            tax.add = "Family")
```

```
## Warning: Transformation introduced infinite values in discrete y-axis
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


## Metacoder comparison

### Original Tree



```r
Tree1
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-9-1.png)<!-- -->



### Plantation vs Native



```r
#Tree visual
set.seed(1)
metacoder::heat_tree_matrix(taxa::filter_taxa(obj, taxon_names == "c__Agaricomycetes", subtaxa = TRUE),
                            dataset = "diff_table",
                            node_size = n_obs, 
                            node_label = taxon_names,
                            node_color = log2_median_ratio,
                            node_color_range = c("#a6611a","#dfc27d","#bdbdbd","#80cdc1","#018571"), 
                            node_color_trans = "linear",
                            node_label_max = 120,
                            node_color_interval = c(-1, 1),
                            edge_color_interval = c(-1, 1),
                            node_size_axis_label = "Number of OTUs",
                            node_color_axis_label = "Log2 ratio median proportions",
                            initial_layout = "reingold-tilford", layout = "davidson-harel")
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

## Top 100 by site type

### Heatmap by Plantation, Invasive, and Native


```r
#Heatmap
ht.type <-  amp_heatmap(data = SOB_data.stdf,
            group = "site_code",
            tax.show = 100,
            scale.seq = 100,
            plot.text.size = 2,
            tax.aggregate = "Genus",
            tax.add = "Order")
ht.type + scale_x_discrete(breaks = c("INV", "null", "PL", "UN"),
                           labels = c("Invasive", "Nullabor", "Plantation", "Native"))
```

```
## Warning: Transformation introduced infinite values in discrete y-axis
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


### Rank abudance of top 80 by forest type

__Plantation__


```r
SOB_data.stdf.PL <- subset_samples(SOB_data.stdf, site_code == "PL")
amp_rabund(SOB_data.stdf.PL,
           tax.aggregate = "Genus", 
           tax.add = "Order",
           scale.seq = 10,
           tax.show = 80,
           adjust.zero = 0.1,
           plot.log = TRUE)
```

```
## 80
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

__Native__


```r
SOB_data.stdf.UN <- subset_samples(SOB_data.stdf, site_code == "UN")
amp_rabund(SOB_data.stdf.UN,
           tax.aggregate = "Genus", 
           tax.add = "Order",
           scale.seq = 10,
           tax.show = 80,
           adjust.zero = 0.1,
           plot.log = TRUE)
```

```
## 80
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

__Invasive__


```r
SOB_data.stdf.INV <- subset_samples(SOB_data.stdf, site_code == "INV")
amp_rabund(SOB_data.stdf.INV,
           tax.aggregate = "Genus", 
           tax.add = "Order",
           scale.seq = 10,
           tax.show = 80,
           adjust.zero = 0.1,
           plot.log = TRUE)
```

```
## 80
```

![](Plantation_vs_Native_files/figure-html/unnamed-chunk-14-1.png)<!-- -->



```r
TopNOTUs <- function(sample,N) {
  names(sort(taxa_sums(sample), TRUE)[1:N])
}

PL.sample <- merge_samples(SOB_data.stdf.PL, "site_code")
INV.sample <- merge_samples(SOB_data.stdf.INV, "site_code")
UN.sample <- merge_samples(SOB_data.stdf.UN, "site_code")

top.PL <- TopNOTUs(PL.sample, 100)
top.INV <- TopNOTUs(INV.sample, 100)
top.UN <- TopNOTUs(UN.sample, 100)

PL.100OTUs <- prune_taxa(top.PL, PL.sample) %>% psmelt()
INV.100OTUs <- prune_taxa(top.INV, INV.sample) %>% psmelt()
UN.100OTUs <- prune_taxa(top.UN, UN.sample) %>% psmelt()


write.csv(PL.100OTUs, file = "top_100OTUs_Plantation.csv")
write.csv(INV.100OTUs, file = "top_100OTUs_Invasive.csv")
write.csv(UN.100OTUs, file = "top_100OTUs_Native.csv")
```


