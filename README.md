# coral-embryo-microbiome

A single stressor experiment exploring the microbial community during embryonic development of *Montipora capitata* corals when exposed to increasing levels of PVC leachate.

This analysis uses:
- QIIME2 verson 2023.9 and version 2023.5
- [qiime2 q2-picrust2-2023.2_0 plugin](https://github.com/picrust/q2-picrust2)
- Maaslin2

## Analysis Documentation

### Beta Diversity and PERMANOVA
- **[PERMANOVA Quick Reference](code/PERMANOVA_QuickRef.md)** - One-page summary of PERMANOVA for microbiome data
- **[PERMANOVA Comprehensive Guide](code/PERMANOVA_README.md)** - Complete theoretical explanation and best practices
- **[PERMANOVA Analysis Code](code/permanova_analysis.qmd)** - Working examples with all distance metrics
- **[Core Metrics Visualization](code/core_metrics_viz.qmd)** - Beta diversity PCoA plots
- **[Compositional Analysis](code/compositional_analysis.qmd)** - CLR transformation and Aitchison distance

### Multi-Omics Integration
- **[Mantel Test Quick Reference](code/MANTEL_QuickRef.md)** - One-page summary of Mantel tests for multi-omics correlation
- **[Mantel Test Guide](code/MANTEL_TEST_GUIDE.md)** - Comprehensive guide for correlation analysis between microbiome and RNA-seq data
- **[Mantel Test Analysis](code/mantel_test_analysis.qmd)** - Executable workflow for microbiome-transcriptome correlation
- **[Mantel Test Example](code/EXAMPLE_mantel_test.md)** - Step-by-step example using actual file paths from this repository

### Other Analyses
- **[Taxa Barplot Guide](code/README_taxa_barplot.md)** - Taxonomic composition visualization
- **[MaAsLin2 Guide](code/MAASLIN2_README.md)** - Differential abundance analysis
- **[Example Usage](code/EXAMPLE_usage.md)** - General workflow examples
