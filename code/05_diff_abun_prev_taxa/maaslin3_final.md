# Differential abundance of bacterial taxa at various taxonomic levels
Sarah Tanja
2026-03-16

- [<span class="toc-section-number">1</span> Background](#background)
- [<span class="toc-section-number">2</span> Install
  MaAsLin3](#install-maaslin3)
- [<span class="toc-section-number">3</span> Load
  Libraries](#load-libraries)
- [<span class="toc-section-number">4</span> Setup](#setup)
  - [<span class="toc-section-number">4.1</span> Set custom ggplot
    theme](#set-custom-ggplot-theme)
  - [<span class="toc-section-number">4.2</span> Set
    colorschemes](#set-colorschemes)
  - [<span class="toc-section-number">4.3</span> Define file
    paths](#define-file-paths)
  - [<span class="toc-section-number">4.4</span> Load
    metadata](#load-metadata)
  - [<span class="toc-section-number">4.5</span> Put read depth per
    sample into the
    metadata](#put-read-depth-per-sample-into-the-metadata)
- [<span class="toc-section-number">5</span> Load Data](#load-data)
  - [<span class="toc-section-number">5.1</span> Feature tables
    L4-7](#feature-tables-l4-7)
- [<span class="toc-section-number">6</span> Run
  MaAslin3](#run-maaslin3)
- [<span class="toc-section-number">7</span> Pivot data longer for
  plotting](#pivot-data-longer-for-plotting)
- [<span class="toc-section-number">8</span> Abundance](#abundance)
  - [<span class="toc-section-number">8.1</span> Species level
    (L7)](#species-level-l7)
  - [<span class="toc-section-number">8.2</span> Family
    (L5)](#family-l5)
  - [<span class="toc-section-number">8.3</span> Order (L4)](#order-l4)
- [<span class="toc-section-number">9</span> Prevalence](#prevalence)
  - [<span class="toc-section-number">9.1</span> Species
    (L7)](#species-l7)
  - [<span class="toc-section-number">9.2</span> Family
    (L5)](#family-l5-1)
  - [<span class="toc-section-number">9.3</span> Order
    (L4)](#order-l4-1)

# Background

H<sub>o</sub> : “PVC Leachate does not alter microbiome trajectory over
developmental time”

H<sub>a</sub> : “PVC Leachate alters microbiome trajectory over
developmental time”

Microbiome Multivariable Associations with Linear Models (MaAsLin)

[MaAslin3](https://huttenhower.sph.harvard.edu/MaAsLin3) improves on
MaAslin2 by accounting for compositionality and accommodates
cross-sectional studies (that’s ours!)

First read and get familiar with the: [MaAslin3 package
README](https://github.com/biobakery/maaslin3) & the [MaAslin3
tutorial](https://github.com/biobakery/biobakery/wiki/maaslin3) from the
Stanford University Huttenhower Lab.

Here we aim to test changes in relative abundance (how many) and
prevalence (presence/absence) of specific taxa due to leachate, and the
interaction of leachate exposure over time, while controlling for spawn
night and read depth.

Our PERMANOVA detected that spawn night was a significant batch effect
and introduced nuisance variation into our data. To account for this
`(1 | spawn_night)` adds a random intercept for each unique spawn night
— that is, each spawn night gets its own baseline shift in abundance,
but the effects of `hpf` and `leachate` is assumed to be the same across
nights.

[MaAsLin3 Wiki](https://github.com/biobakery/biobakery/wiki/maaslin3)
Any significant abundance associations with a categorical variable
should usually have at least 10 observations in each category.
Significant prevalence associations with categorical variables should
also have at least 10 samples in which the feature was present and at
least 10 samples in which it was absent for each significant category.
Significant abundance associations with continuous metadata should be
checked visually for influential outliers.

> [!IMPORTANT]
>
> There are also a few rules of thumb to keep in mind:
>
> - Models should ideally have about **10 times as many samples** (all
>   samples for logistic fits, non-zero samples for linear fits) **as
>   covariate terms** (**all continuous variables plus all categorical
>   variable levels**).
>
> <!-- -->
>
>     We have 63 samples... so the maximum number of terms we should use is 6 
>
> - Significant associations for MaAsLin 3 are results with no model
>   fitting errors, and joint q-value less than 0.1
>
> - Significant abundance associations with continuous metadata should
>   be checked visually for influential outliers.

# Install MaAsLin3

``` r
#library("devtools")
#install_github("biobakery/maaslin3")
```

# Load Libraries

``` r
library(maaslin3)
library(qiime2R)
library(tidyverse)
```

> MaAsLin3 requires two input files, one for taxonomic or functional
> feature abundances, and one for sample metadata.

1.  Data (or features) file

- This file is tab-delimited.
- Formatted with features as columns and samples as rows.
- The transpose of this format is also okay.
- Possible features in this file include **microbes**, genes, pathways,
  etc.

2.  Metadata file

- Formatted with metadata as columns and samples as rownames
- Is a data.frame object
- Includes per sample read counts

# Setup

## Set custom ggplot theme

``` r
library(ggsidekick) # theme by sean anderson
theme_sleek_axe <- function() {
  theme_sleek() +
    theme(
      panel.border = element_rect(
        color = "grey70",
        fill = NA,
        linewidth = 0.5
      ),
      axis.line.x  = element_line(color = "grey70"),
      axis.line.y  = element_line(color = "grey70"),
      
      plot.title = element_text(
        size  = 10,
        color = "grey40",
        hjust = 0.5,
        margin = margin(b = 20)
      ),
      
      strip.text = element_text(size = 8, color = "grey40"),
  
      
      plot.subtitle = element_text(
        size   = 8,
        color  = "grey40",
        hjust  = 0.5,
        margin = margin(b = 40)
      ),
      
      axis.title = element_text(size = 8, color = "grey40"),
      axis.text = element_text(size = 7, color = "grey40"),
      legend.title = element_text(size = 8, color = "grey40"),
      legend.text = element_text(size = 7, color = "grey40")
    )
}
```

## Set colorschemes

``` r
leachate.colors <- c(control = "#AEF1FF", 
                     low     = "#BBC7FF",
                     mid     = "#7D8BFF", 
                     high    = "#592F7D")

stage.5.colors <- c(egg           = "#FFE362",
                    cleavage      = "#EBA600", 
                    morula        = "#E6AA83",
                    prawnchip     = "#D9685B", 
                    earlygastrula = "#A2223C")

stage.3.colors <- c(cleavage      = "#EBA600", 
                    prawnchip     = "#D9685B", 
                    earlygastrula = "#A2223C")

status.colors <- c(typical   = "#75C165", 
                   uncertain = "#E3FAA5", 
                   malformed = "#8B0069")

night.colors <- c(July_6th = "#E3FAA5", 
                  July_7th = "#578B21", 
                  July_8th = "#1E2440")

night.colors.alt <- c(July_6th = "#21918C", 
                      July_7th = "#201158", 
                      July_8th = "#1E2440")
```

## Define file paths

``` r
feature_table_path <- "../../salipante/241121_StonyCoral/270x200/"
metadata_path <- "../../metadata/meta.csv"
output_path <- "../../output/maaslin/"
fig_path <- "../../output/figs"
```

## Load metadata

``` r
# Load metadata
metadata <- read_csv(metadata_path)

# set factors
metadata <- metadata %>% 
  mutate(
    collection_date = as.Date(collection_date, format = "%d-%b-%Y"),
    stage    = factor(stage,    levels = c("cleavage", "prawnchip", "earlygastrula"), ordered = TRUE),
    leachate = factor(leachate, levels = c("control", "low", "mid", "high"),        ordered = TRUE),
    spawn_night = factor(
      collection_date,
      levels  = as.Date(c("06-Jul-2024", "07-Jul-2024", "08-Jul-2024"), format = "%d-%b-%Y"),
      labels  = c("July 6th", "July 7th", "July 8th"),
      ordered = TRUE
    )
  )
```

## Put read depth per sample into the metadata

> Because MaAsLin 3 identifies prevalence (presence/absence)
> associations, sample read depth (number of reads) should be included
> as a covariate if available. Deeper sequencing will likely increase
> feature detection in a way that could spuriously correlate with
> metadata of interest when read depth is not included in the model.

``` r
readepth <- read_csv("../../salipante/Sarah_StonyCoral/241121_StonyCoral_readcounts.csv")
```

``` r
metadata <- metadata %>% 
  left_join(readepth, by = c("sample_id" = "sample"))
```

``` r
# make samples rownames
meta <- metadata %>% 
  column_to_rownames(var = "sample_id") %>% 
  as.data.frame()

# View metadata structure
str(meta)
```

    'data.frame':   63 obs. of  10 variables:
     $ collection_date: Date, format: "2024-07-08" "2024-07-08" ...
     $ parents        : num  101112 101112 101112 101112 101112 ...
     $ group          : chr  "C14" "C4" "C9" "H14" ...
     $ hpf            : num  14 4 9 14 4 9 14 4 9 14 ...
     $ stage          : Ord.factor w/ 3 levels "cleavage"<"prawnchip"<..: 3 1 2 3 1 2 3 1 2 3 ...
     $ leachate       : Ord.factor w/ 4 levels "control"<"low"<..: 1 1 1 4 4 4 2 2 2 3 ...
     $ leachate_mgL   : num  0 0 0 1 1 1 0.01 0.01 0.01 0.1 ...
     $ spawn_night    : Ord.factor w/ 3 levels "July 6th"<"July 7th"<..: 3 3 3 3 3 3 3 3 3 3 ...
     $ reads          : num  1315572 1691582 1663949 1339947 2341461 ...
     $ % of reads     : chr  "0.99%" "1.27%" "1.25%" "1.01%" ...

# Load Data

## Feature tables L4-7

``` r
# L4 Order
## Load feature table from QIIME2 artifact
ft_l4 <- read_qza(file.path(feature_table_path, "collapsed-l4.qza"))$data

## Convert to data.frame
ft_l4 <- ft_l4 %>% 
  as.data.frame()

# L5 Family
## Load feature table from QIIME2 artifact
ft_l5 <- read_qza(file.path(feature_table_path, "collapsed-l5.qza"))$data

## Convert to data.frame
ft_l5 <- ft_l5 %>% 
  as.data.frame()

# L6 Genus
## Load feature table from QIIME2 artifact
ft_l6 <- read_qza(file.path(feature_table_path, "collapsed-l6.qza"))$data

## Convert to data.frame
ft_l6 <- ft_l6 %>% 
  as.data.frame()

# L7 Species
## Load feature table from QIIME2 artifact
ft_l7 <- read_qza(file.path(feature_table_path, "collapsed-l7.qza"))$data

## Convert to data.frame
ft_l7 <- ft_l7 %>% 
  as.data.frame()
```

# Run MaAslin3

> What if we look only for the main effect of leachate and the
> interaction between leachate and time, but exclude the coefficient for
> hpf itself? That would let us focus on features that show significant
> differences due either to leachate or the interaction between leachate
> and time, without requiring a significant linear trend with time
> itself. That might best capture features that respond to leachate in a
> non-linear way across development, while filtering out the features
> that are only responding to time.

**Model formula:**
`Abundance|Prevelance ~ leachate + leachate:hpf + reads + (1|spawn_night)`

**This model tests:**

- Main effect of `leachate` (categorical; control, low, mid, high):
  Differences in the response between the levels of the `leachate`
  treatment, averaged over all `hpf`values. We choose to represent
  leachate categorically instead of a continuous variable here because
  we do not assume that response to increasing leachate levels is
  linear, and we want to capture any non-linear patterns of response
  across the different leachate treatments.

- Interaction between `leachate` and `hpf`: Whether the slope of the
  response across hpf differs depending on the leachate category. This
  tests if `leachate` modifies how the response changes continuously
  over time; e.g., some leachate levels may accelerate or slow the
  effect of `hpf`.

- Effect of `reads` (continuous covariate): Controls for the influence
  of read depth on the response, adjusting other effects accordingly.

- Random intercept for `spawn_night`: Accounts for variability and
  non-independence among samples collected on the same spawn night by
  allowing different baseline response levels for each spawn night.

Interpretation-wise, if the interaction term is significant, it
indicates that the “trajectory” or slope of response change with `hpf`
differs by `leachate` group. If the interaction is not significant but
`leachate` is, then `leachate` affects overall response level but not
the shape of the response over time.

Graphically, imagine separate regression lines (response vs. `hpf`) for
each `leachate` group. The interaction tests whether these lines have
different slopes. The main effects test differences in intercepts and
general trend.

This model lets you hone in on how `leachate` influences response
patterns over continuous developmental time (`hpf`), adjusting for
sequencing depth (`reads`) and including random intercepts for
`spawn_night`.

``` r
set.seed(03132026)
maaslin_fit <- maaslin3(
  input_data = ft_l7,
  input_metadata = meta,
  output = file.path(output_path, "Interaction/L7_leachatehpf_cat.con"),
  formula = ~ leachate + leachate:hpf + reads + (1|spawn_night),
  verbosity = 'ERROR',
  cores = 3, 
  plot_associations = FALSE
)
```

``` r
L7_sig <- read_tsv(file.path(output_path, "Interaction/L7_leachatehpf_cat.con/significant_results.tsv"))

L7_normalized<- read_tsv(file.path(output_path, "Interaction/L7_leachatehpf_cat.con/features/data_norm.tsv"))

L7_transformed <- read_tsv(file.path(output_path, "Interaction/L7_leachatehpf_cat.con/features/data_transformed.tsv"))

L7_filtered <- read_tsv(file.path(output_path, "Interaction/L7_leachatehpf_cat.con/features/filtered_data.tsv"))
```

``` r
L7_sig <- L7_sig %>% 
  mutate(
    CI_lwr = coef - 1.98 * stderr,
    CI_upr = coef + 1.98 * stderr
  ) %>% 
  filter(is.na(error)) %>% 
  filter(metadata != "reads")

L7_feat <- unique(L7_sig$feature)
```

There are 57 unique features

``` r
L7_sig %>% 
  group_by(model) %>%
  summarise(n = n_distinct(feature))
```

    # A tibble: 2 × 2
      model          n
      <chr>      <int>
    1 abundance     57
    2 prevalence    25

``` r
# Extract only unique features to a character vector for prevalence and abundance models 
feat_prev <- L7_sig %>% 
  filter(model == "prevalence") %>% 
  distinct(feature) %>% 
  pull(feature)

feat_abun <- L7_sig %>% 
  filter(model == "abundance") %>% 
  distinct(feature) %>% 
  pull(feature)

# Which overlap?
feat_overlap <- intersect(feat_prev, feat_abun)

# How many unique features are in each model?
# Only in prevalence model
paste0("There are ", length(feat_prev), " total and ", length(setdiff(feat_prev, feat_abun)), " unique significant features in prevalence model")
```

    [1] "There are 25 total and 0 unique significant features in prevalence model"

``` r
# Only in abundance model
paste0("There are ", length(feat_abun), " total and ", length(setdiff(feat_abun, feat_prev)), " unique significant features in abundance model" )
```

    [1] "There are 57 total and 32 unique significant features in abundance model"

``` r
# In both models
paste0("There are ", length(feat_overlap), " significant features in both")
```

    [1] "There are 25 significant features in both"

> [!IMPORTANT]
>
> ALL taxa that are significant in the prevalence model are also
> significant in the abundance model, but there are some taxa that are
> only significant in the abundance model and not the prevalence model.
> This suggests that some taxa show significant changes in relative
> abundance across leachate and time, but do not show significant
> changes in presence/absence. In other words, these taxa may be
> consistently present across samples but vary in how abundant they are,
> rather than being present in some samples and absent in others.
> Meanwhile, all taxa that show significant changes in presence/absence
> also show significant changes in relative abundance, which makes sense
> because if a taxon is present in some samples and absent in others, it
> is likely to also show differences in how abundant it is across those
> samples.

> In abundance models, a one-unit change in the metadatum variable
> corresponds to a 2coef fold change in the relative abundance of the
> feature.

> In prevalence models, a one-unit change in the metadatum variable
> corresponds to a coef change in the log-odds of a feature being
> present.

# Pivot data longer for plotting

##### Make function to extract taxa levels

``` r
extract_tax_level <- function(x, rank) {
  # rank like "d", "p", "c", "o", "f", "g", "s"
  str_extract(x, paste0("(?<=", rank, "__)[^;]+"))
}
```

example to use funtion like: df2 \<- df %\>% mutate( domain =
extract_tax_level(feature, “d”), phylum = extract_tax_level(feature,
“p”), class = extract_tax_level(feature, “c”), order =
extract_tax_level(feature, “o”), family = extract_tax_level(feature,
“f”), genus = extract_tax_level(feature, “g”), species =
extract_tax_level(feature, “s”) )

``` r
L7_sig_long <- L7_normalized %>% 
  as_tibble() %>% 
  rename(sample_id = feature) %>% 
  pivot_longer(
    cols = 2:last_col(),
    names_to = "feature",
    values_to = "norm_abundance"
  ) %>% 
  filter(feature %in% L7_feat) %>%
  replace_na(list(norm_abundance = 0)) %>%
  left_join(metadata, by = "sample_id") %>% 
  mutate(
    order = extract_tax_level(feature, "o"),
    family = extract_tax_level(feature, "f"), 
    genus = extract_tax_level(feature, "g"),
    species = extract_tax_level(feature, "s"), 
    genus_species = paste0(genus, " ", species)
  ) %>% 
  mutate(presence = ifelse(norm_abundance > 0, 1, 0))

L7_feat_all <- unique(L7_sig_long$genus_species)
```

# Abundance

## Species level (L7)

``` r
## Relative abundance (normalized read counts)
L7_sig_long %>% 
ggplot(., aes(x = factor(hpf), y = norm_abundance, color = leachate, fill = leachate)) +
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.8, position = position_dodge(width = 0.8)) +
  facet_wrap(~genus_species, scales = "free_y") +
  scale_fill_manual(values = leachate.colors) +
  scale_color_manual(values = leachate.colors) +
  labs(
      title = "Relative Abundance of leachate responsive bacterial Species in 14 hours of embryonic development",
      x = "Hours Post Fertilization",
      y = "Relative Abundance (normalized read counts)"
    ) +
  theme_sleek_axe()
```

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-18-1.png)

``` r
L7_sig_long %>%
  group_split(genus_species) %>%
  walk(function(df) {

    p <- ggplot(df, aes(x = factor(hpf), y = norm_abundance,
                        color = leachate, fill = leachate)) +
      geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA,
                   position = position_dodge(width = 0.8)) +
      geom_jitter(alpha = 0.8,
                  position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = leachate.colors) +
      scale_color_manual(values = leachate.colors) +
      labs(
        title = paste("Relative abundance of", unique(df$genus_species)),
        x = "Hours Post Fertilization",
        y = "Relative Abundance (normalized read counts)"
      ) +
      theme_sleek_axe()

    ggsave(
      filename = file.path(fig_path, paste0("species_abundance_boxplots/L7_", unique(df$genus_species), ".png")),
      plot = p,
      width = 6,
      height = 4,
      dpi = 900
    )
    
    print(p)

  })
```

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-1.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-2.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-3.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-4.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-5.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-6.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-7.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-8.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-9.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-10.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-11.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-12.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-13.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-14.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-15.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-16.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-17.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-18.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-19.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-20.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-21.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-22.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-23.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-24.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-25.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-26.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-27.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-28.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-29.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-30.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-31.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-32.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-33.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-34.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-35.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-36.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-37.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-38.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-39.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-40.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-41.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-42.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-43.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-44.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-45.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-46.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-47.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-48.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-49.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-50.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-51.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-52.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-53.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-19-54.png)

## Family (L5)

``` r
L7_sig_long %>%
  group_split(family) %>%
  walk(function(df) {

    p <- ggplot(df, aes(x = factor(hpf), y = norm_abundance,
                        color = leachate, fill = leachate)) +
      geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA,
                   position = position_dodge(width = 0.8)) +
      geom_jitter(alpha = 0.8,
                  position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = leachate.colors) +
      scale_color_manual(values = leachate.colors) +
      labs(
        title = paste("Relative abundance of", unique(df$family)),
        x = "Hours Post Fertilization",
        y = "Relative Abundance (normalized read counts)"
      ) +
      theme_sleek_axe()

    ggsave(
      filename = file.path(fig_path, paste0("family_abundance_boxplots/L7_", unique(df$family), ".png")),
      plot = p,
      width = 6,
      height = 4,
      dpi = 900
    )
    
    print(p)

  })
```

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-1.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-2.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-3.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-4.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-5.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-6.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-7.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-8.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-9.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-10.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-11.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-12.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-13.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-14.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-15.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-16.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-17.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-18.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-19.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-20.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-21.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-22.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-20-23.png)

``` r
## Relative abundance (normalized read counts)
L7_sig_long %>% 
ggplot(., aes(x = factor(hpf), y = norm_abundance, color = leachate, fill = leachate)) +
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.8, position = position_dodge(width = 0.8)) +
  facet_wrap(~family, scales = "free_y") +
  scale_fill_manual(values = leachate.colors) +
  scale_color_manual(values = leachate.colors) +
  labs(
      title = "Relative Abundance of leachate responsive bacterial Famiy in 14 hours of embryonic development",
      x = "Hours Post Fertilization",
      y = "Relative Abundance (normalized read counts)"
    ) +
  theme_sleek_axe()
```

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-21-1.png)

``` r
ggsave(
  filename = file.path(fig_path, "family_abundance.png"),
  width = 12,
  height = 12,
  dpi = 900
)
```

## Order (L4)

``` r
L7_sig_long %>%
  group_split(order) %>%
  walk(function(df) {

    p <- ggplot(df, aes(x = factor(hpf), y = norm_abundance,
                        color = leachate, fill = leachate)) +
      geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA,
                   position = position_dodge(width = 0.8)) +
      geom_jitter(alpha = 0.8,
                  position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = leachate.colors) +
      scale_color_manual(values = leachate.colors) +
      labs(
        title = paste("Relative abundance of", unique(df$order)),
        y = "Relative Abundance (normalized read counts)",
      ) +
      theme_sleek_axe()

    ggsave(
      filename = file.path(fig_path, paste0("order_abundance_boxplots/L7_", unique(df$order), ".png")),
      plot = p,
      width = 6,
      height = 4,
      dpi = 900
    )
    
    print(p)

  })
```

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-1.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-2.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-3.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-4.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-5.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-6.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-7.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-8.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-9.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-10.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-11.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-12.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-13.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-14.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-15.png)

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-23-16.png)

``` r
## Relative abundance (normalized read counts)
L7_sig_long %>% 
ggplot(., aes(x = factor(hpf), y = norm_abundance, color = leachate, fill = leachate)) +
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.8, position = position_dodge(width = 0.8)) +
  facet_wrap(~order, scales = "free_y") +
  scale_fill_manual(values = leachate.colors) +
  scale_color_manual(values = leachate.colors) +
  labs(
      title = "Relative Abundance of leachate responsive bacterial Order in 14 hours of embryonic development",
      x = "Hours Post Fertilization",
      y = "Relative Abundance (normalized read counts)"
    ) +
  theme_sleek_axe()
```

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-24-1.png)

``` r
ggsave(
  filename = file.path(fig_path, "order_abundance.png"),
  width = 12,
  height = 12,
  dpi = 900
)
```

# Prevalence

## Species (L7)

``` r
L7_sig_long %>% filter(feature %in% feat_prev) %>%
ggplot(.,
       aes(x = leachate,
           y = as.numeric(presence),
           color = leachate,
           fill = leachate)) +
    
  geom_smooth(
    aes(group = 1),
    method = "glm",
    method.args = list(family = "binomial"),
    se = FALSE,
    color = "grey40",
    fill = "grey70",
    linewidth = 0.8
  ) +
  geom_rug(
    aes(x = leachate_mgL),
    sides = "bt",
    color = "grey70",
    alpha = 0.5
  ) +
  geom_dotplot(
    binaxis = "y",
    stackdir = "centerwhole",
    alpha = 0.5,
    shape = 1,
    dotsize = 1.5
  ) +
    # mean points
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4
  ) +

  facet_grid(genus_species~hpf) +
  coord_fixed(ratio = 2.8) +
  scale_fill_manual(values = leachate.colors) +
  scale_color_manual(values = leachate.colors) +
  labs(
    title = "Prevalence of leachate responsive bacteria Species across leachate levels and developmental stages",
    x = "Leachate Level",
    y = "Prevalence (proportion of samples with bacteria Species present)"
  ) +
  theme_sleek_axe()+ 
  theme(
  strip.text.y = element_text(size = 8, color = "grey40", angle = 0, hjust = 0)
)
```

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-26-1.png)

``` r
ggsave(
  filename = file.path(fig_path, "species_prevalence.png"),
  width = 12,
  height = 12,
  dpi = 900
)
```

## Family (L5)

``` r
L7_sig_long %>% filter(feature %in% feat_prev) %>%
ggplot(.,
       aes(x = leachate,
           y = as.numeric(presence),
           color = leachate,
           fill = leachate)) +
    
  geom_smooth(
    aes(group = 1),
    method = "glm",
    method.args = list(family = "binomial"),
    se = FALSE,
    color = "grey40",
    fill = "grey70",
    linewidth = 0.8
  ) +
  geom_rug(
    aes(x = leachate_mgL),
    sides = "bt",
    color = "grey70",
    alpha = 0.5
  ) +
  geom_dotplot(
    binaxis = "y",
    stackdir = "centerwhole",
    alpha = 0.5,
    shape = 1,
    dotsize = 1.5
  ) +
    # mean points
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4
  ) +

  facet_grid(family~hpf) +
  coord_fixed(ratio = 2.8) +
  scale_fill_manual(values = leachate.colors) +
  scale_color_manual(values = leachate.colors) +
  labs(
    title = "Prevalence of leachate responsive bacteria Family across leachate levels and developmental stages",
    x = "Leachate Level",
    y = "Prevalence (proportion of samples with bacteria Species present)"
  ) +
  theme_sleek_axe()+ 
  theme(
  strip.text.y = element_text(size = 8, color = "grey40", angle = 0, hjust = 0)
)
```

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-28-1.png)

``` r
ggsave(
  filename = file.path(fig_path, "family_prevalence.png"),
  width = 14,
  height = 30,
  dpi = 900)
```

## Order (L4)

``` r
L7_sig_long %>% filter(feature %in% feat_prev) %>%
ggplot(.,
       aes(x = leachate,
           y = as.numeric(presence),
           color = leachate,
           fill = leachate)) +
    
  geom_smooth(
    aes(group = 1),
    method = "glm",
    method.args = list(family = "binomial"),
    se = FALSE,
    color = "grey40",
    fill = "grey70",
    linewidth = 0.8
  ) +
  geom_rug(
    aes(x = leachate_mgL),
    sides = "bt",
    color = "grey70",
    alpha = 0.5
  ) +
  geom_dotplot(
    binaxis = "y",
    stackdir = "centerwhole",
    alpha = 0.5,
    shape = 1,
    dotsize = 1.5
  ) +
    # mean points
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4
  ) +

  facet_grid(order~hpf) +
  coord_fixed(ratio = 2.8) +
  scale_fill_manual(values = leachate.colors) +
  scale_color_manual(values = leachate.colors) +
  labs(
    title = "Prevalence of leachate responsive bacteria Order across leachate levels and developmental stages",
    x = "Leachate Level",
    y = "Prevalence (proportion of samples with bacteria Species present)"
  ) +
  theme_sleek_axe()+ 
  theme(
  strip.text.y = element_text(size = 8, color = "grey40", angle = 0, hjust = 0)
)
```

![](maaslin3_final_files/figure-commonmark/unnamed-chunk-30-1.png)

``` r
ggsave(
  filename = file.path(fig_path, "order_prevalence.png"),
  width = 12,
  height = 18,
  dpi = 900
)
```

``` r
sessionInfo()
```

    R version 4.5.1 (2025-06-13 ucrt)
    Platform: x86_64-w64-mingw32/x64
    Running under: Windows 11 x64 (build 26200)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=English_United States.utf8 
    [2] LC_CTYPE=English_United States.utf8   
    [3] LC_MONETARY=English_United States.utf8
    [4] LC_NUMERIC=C                          
    [5] LC_TIME=English_United States.utf8    

    time zone: America/Los_Angeles
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] ggsidekick_0.0.3 lubridate_1.9.4  forcats_1.0.1    stringr_1.6.0   
     [5] dplyr_1.1.4      purrr_1.2.1      readr_2.1.6      tidyr_1.3.2     
     [9] tibble_3.3.1     ggplot2_4.0.1    tidyverse_2.0.0  qiime2R_0.99.6  
    [13] maaslin3_1.0.2  

    loaded via a namespace (and not attached):
      [1] RColorBrewer_1.1-3              rstudioapi_0.18.0              
      [3] jsonlite_2.0.0                  magrittr_2.0.4                 
      [5] farver_2.1.2                    rmarkdown_2.30                 
      [7] fs_1.6.6                        ragg_1.5.0                     
      [9] vctrs_0.6.5                     multtest_2.64.0                
     [11] base64enc_0.1-3                 htmltools_0.5.9                
     [13] S4Arrays_1.8.1                  truncnorm_1.0-9                
     [15] Rhdf5lib_1.30.0                 SparseArray_1.8.1              
     [17] Formula_1.2-5                   rhdf5_2.52.1                   
     [19] htmlwidgets_1.6.4               plyr_1.8.9                     
     [21] igraph_2.2.1                    lifecycle_1.0.5                
     [23] iterators_1.0.14                pkgconfig_2.0.3                
     [25] Matrix_1.7-4                    R6_2.6.1                       
     [27] fastmap_1.2.0                   GenomeInfoDbData_1.2.14        
     [29] MatrixGenerics_1.20.0           digest_0.6.39                  
     [31] colorspace_2.1-2                S4Vectors_0.46.0               
     [33] textshaping_1.0.4               Hmisc_5.2-5                    
     [35] GenomicRanges_1.60.0            vegan_2.7-2                    
     [37] labeling_0.4.3                  timechange_0.3.0               
     [39] httr_1.4.7                      TreeSummarizedExperiment_2.16.1
     [41] abind_1.4-8                     mgcv_1.9-4                     
     [43] compiler_4.5.1                  bit64_4.6.0-1                  
     [45] withr_3.0.2                     htmlTable_2.4.3                
     [47] S7_0.2.1                        backports_1.5.0                
     [49] BiocParallel_1.42.1             MASS_7.3-65                    
     [51] rappdirs_0.3.3                  DelayedArray_0.34.1            
     [53] biomformat_1.36.0               permute_0.9-8                  
     [55] optparse_1.7.5                  tools_4.5.1                    
     [57] foreign_0.8-90                  otel_0.2.0                     
     [59] ape_5.8-1                       nnet_7.3-20                    
     [61] glue_1.8.0                      nlme_3.1-168                   
     [63] rhdf5filters_1.20.0             grid_4.5.1                     
     [65] checkmate_2.3.3                 cluster_2.1.8.1                
     [67] reshape2_1.4.5                  ade4_1.7-23                    
     [69] generics_0.1.4                  gtable_0.3.6                   
     [71] tzdb_0.5.0                      data.table_1.18.0              
     [73] hms_1.1.4                       utf8_1.2.6                     
     [75] XVector_0.48.0                  BiocGenerics_0.54.1            
     [77] foreach_1.5.2                   pillar_1.11.1                  
     [79] yulab.utils_0.2.3               vroom_1.6.7                    
     [81] splines_4.5.1                   getopt_1.20.4                  
     [83] treeio_1.32.0                   lattice_0.22-7                 
     [85] survival_3.8-3                  bit_4.6.0                      
     [87] tidyselect_1.2.1                SingleCellExperiment_1.30.1    
     [89] Biostrings_2.76.0               knitr_1.51                     
     [91] gridExtra_2.3                   phyloseq_1.52.0                
     [93] IRanges_2.42.0                  SummarizedExperiment_1.38.1    
     [95] zCompositions_1.5.0-5           stats4_4.5.1                   
     [97] xfun_0.54                       Biobase_2.68.0                 
     [99] matrixStats_1.5.0               DT_0.34.0                      
    [101] stringi_1.8.7                   UCSC.utils_1.4.0               
    [103] lazyeval_0.2.2                  yaml_2.3.12                    
    [105] evaluate_1.0.5                  codetools_0.2-20               
    [107] cli_3.6.5                       rpart_4.1.24                   
    [109] systemfonts_1.3.1               dichromat_2.0-0.1              
    [111] Rcpp_1.1.1                      GenomeInfoDb_1.44.3            
    [113] parallel_4.5.1                  tidytree_0.4.7                 
    [115] scales_1.4.0                    crayon_1.5.3                   
    [117] rlang_1.1.6                    
