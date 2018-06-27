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
## Run 1 stress 0.3271709 
## ... New best solution
## ... Procrustes: rmse 0.034605  max resid 0.2062421 
## Run 2 stress 0.3231505 
## ... New best solution
## ... Procrustes: rmse 0.0247594  max resid 0.2260914 
## Run 3 stress 0.4182684 
## Run 4 stress 0.3340041 
## Run 5 stress 0.3242311 
## Run 6 stress 0.3235061 
## ... Procrustes: rmse 0.02189246  max resid 0.2233902 
## Run 7 stress 0.3283945 
## Run 8 stress 0.3281971 
## Run 9 stress 0.3253827 
## Run 10 stress 0.3247226 
## Run 11 stress 0.3247521 
## Run 12 stress 0.3363564 
## Run 13 stress 0.3271408 
## Run 14 stress 0.3237851 
## Run 15 stress 0.3265173 
## Run 16 stress 0.3271274 
## Run 17 stress 0.3291448 
## Run 18 stress 0.3257645 
## Run 19 stress 0.3244434 
## Run 20 stress 0.3327246 
## *** No convergence -- monoMDS stopping criteria:
##      2: no. of iterations >= maxit
##     18: stress ratio > sratmax
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
## Run 1 stress 0.2876689 
## ... New best solution
## ... Procrustes: rmse 0.01531088  max resid 0.1597298 
## Run 2 stress 0.2890027 
## Run 3 stress 0.2889662 
## Run 4 stress 0.2891255 
## Run 5 stress 0.2931497 
## Run 6 stress 0.2900311 
## Run 7 stress 0.2891465 
## Run 8 stress 0.2909815 
## Run 9 stress 0.292933 
## Run 10 stress 0.2880395 
## ... Procrustes: rmse 0.009450927  max resid 0.1145649 
## Run 11 stress 0.2906672 
## Run 12 stress 0.2888722 
## Run 13 stress 0.2882818 
## Run 14 stress 0.2919494 
## Run 15 stress 0.2910357 
## Run 16 stress 0.2904799 
## Run 17 stress 0.2882752 
## Run 18 stress 0.2928396 
## Run 19 stress 0.288037 
## ... Procrustes: rmse 0.009386043  max resid 0.1141034 
## Run 20 stress 0.2892791 
## *** No convergence -- monoMDS stopping criteria:
##      1: no. of iterations >= maxit
##     19: stress ratio > sratmax
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
## Run 1 stress 0.2884708 
## ... New best solution
## ... Procrustes: rmse 0.0533707  max resid 0.2266338 
## Run 2 stress 0.2955162 
## Run 3 stress 0.2954637 
## Run 4 stress 0.2931377 
## Run 5 stress 0.2981559 
## Run 6 stress 0.2882965 
## ... New best solution
## ... Procrustes: rmse 0.01458933  max resid 0.1022946 
## Run 7 stress 0.2996505 
## Run 8 stress 0.2926715 
## Run 9 stress 0.3012843 
## Run 10 stress 0.2885983 
## ... Procrustes: rmse 0.02435067  max resid 0.2509616 
## Run 11 stress 0.2923056 
## Run 12 stress 0.2886835 
## ... Procrustes: rmse 0.02704382  max resid 0.2445538 
## Run 13 stress 0.2898205 
## Run 14 stress 0.2971962 
## Run 15 stress 0.2952565 
## Run 16 stress 0.2901546 
## Run 17 stress 0.2931713 
## Run 18 stress 0.2883606 
## ... Procrustes: rmse 0.003211131  max resid 0.02041259 
## Run 19 stress 0.2895781 
## Run 20 stress 0.2889773 
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
## Run 1 stress 0.3204581 
## ... New best solution
## ... Procrustes: rmse 0.05730455  max resid 0.2222996 
## Run 2 stress 0.3186165 
## ... New best solution
## ... Procrustes: rmse 0.04564094  max resid 0.2216435 
## Run 3 stress 0.3238501 
## Run 4 stress 0.318702 
## ... Procrustes: rmse 0.05246132  max resid 0.3038673 
## Run 5 stress 0.3169293 
## ... New best solution
## ... Procrustes: rmse 0.03768056  max resid 0.226009 
## Run 6 stress 0.3278974 
## Run 7 stress 0.3151927 
## ... New best solution
## ... Procrustes: rmse 0.02801938  max resid 0.1707524 
## Run 8 stress 0.3157709 
## Run 9 stress 0.3182548 
## Run 10 stress 0.3313282 
## Run 11 stress 0.3237317 
## Run 12 stress 0.3254074 
## Run 13 stress 0.3182641 
## Run 14 stress 0.3148973 
## ... New best solution
## ... Procrustes: rmse 0.03175343  max resid 0.3041805 
## Run 15 stress 0.3258368 
## Run 16 stress 0.315013 
## ... Procrustes: rmse 0.03862069  max resid 0.3016616 
## Run 17 stress 0.3161033 
## Run 18 stress 0.3149625 
## ... Procrustes: rmse 0.005119172  max resid 0.03398091 
## Run 19 stress 0.3151797 
## ... Procrustes: rmse 0.01933077  max resid 0.1415934 
## Run 20 stress 0.3264354 
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


