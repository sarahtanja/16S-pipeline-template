# Exploring results from mim_c
Sarah Tanja

- [Install packages & Load libraries](#install-packages--load-libraries)
- [Import data](#import-data)
- [Filter to bacteria that are significantly different in leachate
  treatments](#filter-to-bacteria-that-are-significantly-different-in-leachate-treatments)
- [Collapse duplicates](#collapse-duplicates)
- [Top 10](#top-10)
  - [Arcobacteraceae (family)](#arcobacteraceae-family)
  - [Dolosigranulum (genus)](#dolosigranulum-genus)
  - [Lactobacillus.\_\_ (genus or sp.?) & Lactobacillus
    (genus)](#lactobacillus__-genus-or-sp--lactobacillus-genus)
  - [Microbacteriaceae (family)](#microbacteriaceae-family)
  - [Alteromonadaceae uncultured
    (genus)](#alteromonadaceae-uncultured-genus)
  - [Oleiphilus_sp. (species) & Oleiphilus
    (genus)](#oleiphilus_sp-species--oleiphilus-genus)
  - [Atopostipes.s\_\_uncultured_bacterium
    (species)](#atopostipess__uncultured_bacterium-species)
  - [Cognatishimia.s\_\_uncultured_bacterium
    (species)](#cognatishimias__uncultured_bacterium-species)

The microbiome analysis was completed by Adam Waalkes of mim_c from the
Salipante Lab. The contents of his analysis are in the
`Sarah_StonyCoral` directory.

> A couple of things to note about the shannon diversity analysis we
> did. First, after we filtered out the Uassigned, mitochodria and
> chloroplast the diversity went up in all samples (average from 4.6 to
> 5.7) and we lost statistical significance in timeframe where we had
> had it in the previous analysis. The leachate 0 to 1 was still not
> statistically significant. Also, three of the samples gave errors due
> to a bug in QIIME2 calculating their shannon score so they were not in
> this analysis. Those three samples were in the rest of the analysis. -
> email correspondance with A.W.

# Install packages & Load libraries

``` r
library(tidyverse)
```

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ✔ purrr     1.0.2     
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

Let’s jump straight to significant results from the Maaslin2 significant
differential features

# Import data

``` r
significant <- read_tsv(file = "../Sarah_StonyCoral/Level7_filtered_organism_output/significant_results.tsv")
```

    Rows: 732 Columns: 9
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (3): feature, metadata, value
    dbl (6): coef, stderr, N, N.not.0, pval, qval

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

There were 732 significantly different bacteria across leachate level
(`leachate`) and developmental time in hours post fertilization
(`timepoint`)

# Filter to bacteria that are significantly different in leachate treatments

``` r
sig_leach <- significant %>% 
  filter(metadata == "leachate")
```

125 bacteria were significantly different in leachate treatments. It
seems there are many duplicates.. how do we consolidate? What are these
bacteria, and what are their functions?

Also… it seems we’ve lost all granularity on when the shifts happen and
at which concentrations. I’ll need to circle back to Adam to clarify.

# Collapse duplicates

A feature that is showing up as significant on multiple levels has the
same `coef`, `stderr`, `pval` and `qval` values.

The following code will: - Groups rows that have the same coef, stderr,
pval, and qval - Within each group, selects the feature with the longest
character length using nchar() and which.max()

The longest character length of each group is the deepest taxonomic
level

``` r
collapsed_sig_leach <- sig_leach %>% 
    group_by(coef, stderr, pval, qval) %>%
  summarise(
    feature = feature[which.max(nchar(feature))],
    .groups = "drop"
  )
```

Now we’re down to 91 unique features …

# Top 10

Let’s just quickly explore the top 10 most significant ones

``` r
top10 <- collapsed_sig_leach %>% 
  arrange(desc(abs(coef))) %>% # sort by effect size (magnitude)
  slice_head(n = 10) # keep top 10     

print(top10$feature)
```

     [1] "d__Bacteria.p__Campilobacterota.c__Campylobacteria.o__Campylobacterales.f__Arcobacteraceae.__.__"                                    
     [2] "d__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Carnobacteriaceae.g__Dolosigranulum.__"                                   
     [3] "d__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Lactobacillaceae.g__Lactobacillus.__"                                     
     [4] "d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Micrococcales.f__Microbacteriaceae"                                             
     [5] "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae.g__uncultured"                           
     [6] "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Oceanospirillales.f__Oleiphilaceae.g__Oleiphilus.s__Oleiphilus_sp."          
     [7] "d__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Lactobacillaceae.g__Lactobacillus"                                        
     [8] "d__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Carnobacteriaceae.g__Atopostipes.s__uncultured_bacterium"                 
     [9] "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Oceanospirillales.f__Oleiphilaceae.g__Oleiphilus"                            
    [10] "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhodobacterales.f__Rhodobacteraceae.g__Cognatishimia.s__uncultured_bacterium"

## Arcobacteraceae (family)

> The phylogenomic tree of *Arcobacteraceae* is divided into three large
> clades, among which members of clades A and B are almost all from
> terrestrial environments, while those of clade C are widely
> distributed in various marine habitats in addition to some terrestrial
> origins. All clades harbor genes putatively involved in chitin
> degradation, sulfide oxidation, hydrogen oxidation, thiosulfate
> oxidation, denitrification, dissimilatory nitrate reduction to
> ammonium, microaerophilic respiration, and metal (iron/manganese)
> reduction. Additionally, in clade C, more unique pathways were
> retrieved, including thiosulfate disproportionation, ethanol
> fermentation, methane oxidation, fatty acid oxidation, cobalamin
> synthesis, and dissimilatory reductions of sulfate, perchlorate, and
> arsenate. Within this clade, two mixotrophic Candidatus genera
> represented by UBA6211 and CAIJNA01 harbor genes putatively involved
> in the reverse tricarboxylic acid pathway for carbon fixation.
> Moreover, the metatranscriptomic data in deep-sea *in situ*
> incubations indicated that the latter genus is a mixotroph that
> conducts carbon fixation by coupling sulfur oxidation and
> denitrification and metabolizing organic matter. Furthermore, global
> metatranscriptomic data confirmed the ubiquitous distribution and
> global relevance of *Arcobacteraceae* in the expression of those
> corresponding genes across all oceanic regions and depths.
>
> - (Li et al. 2024)

## Dolosigranulum (genus)

## Lactobacillus.\_\_ (genus or sp.?) & Lactobacillus (genus)

## Microbacteriaceae (family)

## Alteromonadaceae uncultured (genus)

## Oleiphilus_sp. (species) & Oleiphilus (genus)

## Atopostipes.s\_\_uncultured_bacterium (species)

## Cognatishimia.s\_\_uncultured_bacterium (species)

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Li2024-yr" class="csl-entry">

Li, Jianyang, Shizheng Xiang, Yufei Li, Ruolin Cheng, Qiliang Lai,
Liping Wang, Guizhen Li, Chunming Dong, and Zongze Shao. 2024.
“Arcobacteraceae Are Ubiquitous Mixotrophic Bacteria Playing Important
Roles in Carbon, Nitrogen, and Sulfur Cycling in Global Oceans.”
*mSystems* 9 (July): e0051324.
<https://doi.org/10.1128/msystems.00513-24>.

</div>

</div>
