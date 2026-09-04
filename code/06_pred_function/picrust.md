# Running PICRUSt2
Sarah Tanja

- [Goals](#goals)
- [Setup](#setup)
- [Load packages](#load-packages)
- [Install miniconda](#install-miniconda)
- [Install QIIME2 as a conda
  environment](#install-qiime2-as-a-conda-environment)
  - [2023.5](#20235)
  - [2023.9](#20239)
  - [2024.4](#20244)
- [Activate qiime2 environment in R](#activate-qiime2-environment-in-r)
- [Install picrust2](#install-picrust2)
- [Install q2-picrust2 plugin](#install-q2-picrust2-plugin)
- [PICRUSt2 qiime pipeline](#picrust2-qiime-pipeline)
  - [Feature table (BIOM)](#feature-table-biom)
  - [Sequence file](#sequence-file)
  - [Run
    `qiime picrust2 full-pipeline`](#run-qiime-picrust2-full-pipeline)
    - [Install SEQTK](#install-seqtk)
    - [Reverse compliment](#reverse-compliment)
    - [Run picrust with reverse
      compliment](#run-picrust-with-reverse-compliment)
  - [Lower –p-min-align](#lower-p-min-align)
    - [Troubleshoot error](#troubleshoot-error)
  - [Set `LD_LIBRARY_PATH` when running
    QIIME](#set-ld_library_path-when-running-qiime)
- [Run with `mp` & `epa-ng`](#run-with-mp--epa-ng)

# Goals

> Yes it is reasonable to look at picrust output and see if there are
> functional shifts, but it is not clear how adequate picrust libraries
> are to allow metabolic inference of these particular populations.  I
> would recommend that the actual analysis of metabolic shifts be
> conducted in Maaslin, as before, to account for covariates.
>
> it is possible to give sarah the input files she would need for the
> analysis. \[see
> [https://forum.qiime2.org/t/picrust2-input-files-for-bacteria/27968](https://urldefense.com/v3/__https://forum.qiime2.org/t/picrust2-input-files-for-bacteria/27968__;!!K-Hz7m0Vt54!iZIXAiB2i2QtliEgrFMztEKZHQZuAILNAY54WW_klnDOn0rCNia6tnE-GrFq0efQ8IeJkkAwEJWtIM0$)\].
> For the record, you could isolate the ASV table from DADA2 in qiime
> (these are stored in rep-seqs-dada2.qza, according to
> [https://forum.qiime2.org/t/how-to-generate-an-asv-table-from-dada2-output/9973](https://urldefense.com/v3/__https://forum.qiime2.org/t/how-to-generate-an-asv-table-from-dada2-output/9973__;!!K-Hz7m0Vt54!iZIXAiB2i2QtliEgrFMztEKZHQZuAILNAY54WW_klnDOn0rCNia6tnE-GrFq0efQ8IeJkkAwpWyxhvY$)),
> and convert it to biom format
> \[[https://forum.qiime2.org/t/how-to-generate-an-asv-table-from-dada2-output/9973/2](https://urldefense.com/v3/__https://forum.qiime2.org/t/how-to-generate-an-asv-table-from-dada2-output/9973/2__;!!K-Hz7m0Vt54!iZIXAiB2i2QtliEgrFMztEKZHQZuAILNAY54WW_klnDOn0rCNia6tnE-GrFq0efQ8IeJkkAwc-L0gtk$)\]. 
>
> having said this it may just be easier to use the picrust2 plugin for
> qiime and hand over those outputs, since we already have the inputs
> ready in qiime
>
> [https://github.com/picrust/picrust2/wiki/q2-picrust2-Tutorial](https://urldefense.com/v3/__https://github.com/picrust/picrust2/wiki/q2-picrust2-Tutorial__;!!K-Hz7m0Vt54!iZIXAiB2i2QtliEgrFMztEKZHQZuAILNAY54WW_klnDOn0rCNia6tnE-GrFq0efQ8IeJkkAwVcNVfvQ$)
>
> [https://github.com/picrust/picrust2/wiki/Manually-install-QIIME-2-plugin](https://urldefense.com/v3/__https://github.com/picrust/picrust2/wiki/Manually-install-QIIME-2-plugin__;!!K-Hz7m0Vt54!iZIXAiB2i2QtliEgrFMztEKZHQZuAILNAY54WW_klnDOn0rCNia6tnE-GrFq0efQ8IeJkkAwAif4ef0$)  
>   
> id say generate those outputs and then run them through maaslin to
> look for significant feature abundances. – S.S.

# Setup

# Load packages

``` r
library(knitr)
library(dplyr)
library(reticulate)
```

# Install miniconda

Installation instruction for miniconda are
[here](https://www.anaconda.com/docs/getting-started/miniconda/install#linux) -
create a new directory named “miniconda3” in your home directory. -
download the Linux Miniconda installation script for your chosen chip
architecture and save the script as miniconda.sh in the miniconda3
directory. - run the miniconda.sh installation script in silent mode
using bash. - remove the miniconda.sh installation script file after
installation is complete.

``` bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```

After installiong Miniconda make sure you’re running the latest version
of `conda`.

``` bash
conda update conda
```

# Install QIIME2 as a conda environment

Here I’m installing qiime2 version 2023.9 because it’s the same one used
by the Salipante Lab …but it is using python 3.8 and picrust needs
python 3.9 so there is a incompatible package conflict

## 2023.5

``` bash
wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-linux-conda.yml
conda env create -n qiime2-2023.5 --file qiime2-2023.5-py38-linux-conda.yml
rm qiime2-2023.5-py38-linux-conda.yml
```

## 2023.9

> [!NOTE]
>
> Use this version if possible because it’s the same one the Salipante
> Lab is using!

``` bash
conda env create \
  --name qiime2-amplicon-2023.9 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2023.9/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
```

## 2024.4

``` bash
conda env create \
  --name qiime2-amplicon-2024.2 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2024.2/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
```

# Activate qiime2 environment in R

``` r
use_condaenv(condaenv = "qiime2-amplicon-2023.9", conda = "/home/shared/8TB_HDD_02/stanja/miniconda3/condabin/conda")

# Check successful env loading
py_config()
```

    python:         /home/shared/8TB_HDD_02/stanja/miniconda3/envs/qiime2-amplicon-2023.9/bin/python
    libpython:      /home/shared/8TB_HDD_02/stanja/miniconda3/envs/qiime2-amplicon-2023.9/lib/libpython3.8.so
    pythonhome:     /home/shared/8TB_HDD_02/stanja/miniconda3/envs/qiime2-amplicon-2023.9:/home/shared/8TB_HDD_02/stanja/miniconda3/envs/qiime2-amplicon-2023.9
    version:        3.8.15 | packaged by conda-forge | (default, Nov 22 2022, 08:46:39)  [GCC 10.4.0]
    numpy:          /home/shared/8TB_HDD_02/stanja/miniconda3/envs/qiime2-amplicon-2023.9/lib/python3.8/site-packages/numpy
    numpy_version:  1.24.4

    NOTE: Python version was forced by use_python() function

If this is successful, the first line of output should show that the
Python environment being used is the one in your conda environment path.

Your conda and qiime environments should be listed in your \$PATH:

``` bash
echo $PATH
```

# Install picrust2

This must be done from within the active qiime2 conda environment!

`qiime2` depends on `python 3.8` but `picrust2` uses `python 3.9` … This
results in package conflicts!!! The qiime2 library plugin for
q2-picrust2 warms about this and suggests installing the plugin manually
to avoid the “mysterious incompatibility errors”

> If you run into incompatibility errors … try [installing the plugin
> manually](https://github.com/picrust/picrust2/wiki/Manually-install-QIIME-2-plugin)
> instead.

Run the following code in the Terminal and type ‘y’ when prompted to
first install picrust2.

    conda install -c bioconda -c conda-forge picrust2

# Install q2-picrust2 plugin

Download, unzip, and install the q2-picrust2 qiime2 plugin!

1\. Get and unzip the plugin

``` bash
cd /home/shared/8TB_HDD_02/stanja
wget https://github.com/picrust/q2-picrust2/archive/refs/tags/2024.5_1.tar.gz
tar -xzvf 2024.5_1.tar.gz
```

2.  Move into the plugin directory and `pip install` the `setup.py`
    script.

``` bash
cd /home/shared/8TB_HDD_02/stanja/q2-picrust2-2024.5_1
pip install -e .
```

3.  Refresh the qiime cache

``` bash
qiime dev refresh-cache
```

4.  Test that the plugin is working

``` bash
qiime picrust2 --help
```

# PICRUSt2 qiime pipeline

``` bash
qiime picrust2 full-pipeline --help
```

> The required inputs are –i-table and –i-seq, which need to correspond
> to QIIME 2 artifacts of types FeatureTable\[Frequency\] and
> FeatureData\[Sequence\], respectively. The Feature Table needs to
> contain the abundances of ASVs (i.e. a BIOM table) and the sequence
> file needs to be a FASTA file containing the sequences for each ASV.

> These files correspond to the ASV count table, ASV sequences, and
> metadata for the samples.

## Feature table (BIOM)

located at path
`../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_featureTable_filtered.qza`

``` bash
qiime tools peek ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_featureTable_filtered.qza
```

    UUID:        3cc95575-462c-4b0f-a4b0-b83a259906d9
    Type:        FeatureTable[Frequency]
    Data format: BIOMV210DirFmt

## Sequence file

located at path
`../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_representative-sequences.qza`

``` bash
qiime tools peek ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_representative-sequences.qza
```

    UUID:        25349e49-dda5-4954-8d2e-46bdcf0ef058
    Type:        FeatureData[Sequence]
    Data format: DNASequencesDirectoryFormat

## Run `qiime picrust2 full-pipeline`

``` bash
qiime picrust2 full-pipeline \
   --i-table ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_featureTable_filtered.qza \
   --i-seq ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_representative-sequences.qza \
   --output-dir ../output/picrust2 \
   --p-placement-tool sepp \
   --p-threads 4 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose
```

> [!WARNING]
>
> Standard error of the above failed command: Warning - 1169 input
> sequences aligned poorly to reference sequences (–min_align option
> specified a minimum proportion of 0.8 aligning to reference
> sequences).

See [issue#122](https://github.com/picrust/picrust2/issues/122)

> There’s a good chance your input sequences are on the negative strand
> and aren’t being aligned properly. In case some one may need this, use
> seqtk to get reverse complementary seqs, then using
> picrust2_pipeline.py. seqtk seq -r otus.fa.format.fasta \> otus_rc.fa

### Install SEQTK

Run in terminal and proceed with ‘y’

    conda install -c bioconda seqtk

### Reverse compliment

Representative sequences (.qza → .fasta)

``` bash
qiime tools export --help
```

``` bash
qiime tools export \
  --input-path ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_representative-sequences.qza \
  --output-path ../output/
```

Reverse compliment the sequences with `seqtk`

``` bash
seqtk seq -r ../output/dna-sequences.fasta > ../output/rep_seqs_rc.fasta
```

``` bash
qiime tools list-types
```

``` bash
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ../output/rep_seqs_rc.fasta \
  --output-path ../output/rep_seqs_rc.qza \
  --input-format DNAFASTAFormat
```

``` bash
qiime tools peek ../output/rep_seqs_rc.qza
```

### Run picrust with reverse compliment

``` bash
qiime picrust2 full-pipeline \
   --i-table ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_featureTable_filtered.qza \
   --i-seq ../output/rep_seqs_rc.qza \
   --output-dir ../output/picrust2 \
   --p-placement-tool sepp \
   --p-threads 4 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose
```

> [!WARNING]
>
>     Standard error of the above failed command:
>     Stopping - all 7182 input sequences aligned poorly to reference sequences (--min_align option specified a minimum proportion of 0.8 aligning to reference sequences).
>
> OK this likely indicates that input sequences were originally on the
> positive strand, and the seqtk reverse compliment we made is the
> negative strand.. the number of input sequences that failed to align
> increased from 1169/7182 (16%) to 7182/7182 (100%)

> can also consider lowering the --min_align flag so reads have a better
> chance of matching poorly characterized taxa.

–p-min-align PROPORTION Range(0.0, 1.0) Proportion of the total length
of an input query sequence that must align with reference sequences. Any
sequences with lengths below this value after making an alignment with
reference sequences will be excluded from the placement and all
subsequent steps. \[default: 0.8\]

## Lower –p-min-align

``` bash
qiime picrust2 full-pipeline \
   --i-table ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_featureTable_filtered.qza \
   --i-seq ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_representative-sequences.qza \
   --output-dir ../output/picrust2 \
   --p-placement-tool sepp \
   --p-threads 4 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose \
   --p-min-align 0.1
```

### Troubleshoot error

> [!WARNING]
>
> Error running this command: mv
> /tmp/tmpulv7acg9/picrust2_out/intermediate/place_seqs/sepp_out/output_placement.newick
> /tmp/tmpulv7acg9/picrust2_out/out.tre
>
> Standard error of the above failed command: mv: relocation error:
> /lib/x86_64-linux-gnu/libacl.so.1: symbol getxattr version ATTR_1.0
> not defined in file libattr.so.1 with link time reference

The error indicates there are multiple or corrupted versions between
`libacl.so.1` and `libattr.so.1` on your system. To get around this we
can set the library path when running qiime commands.

## Set `LD_LIBRARY_PATH` when running QIIME

``` bash
LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH \
qiime picrust2 full-pipeline \
   --i-table ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_featureTable_filtered.qza \
   --i-seq ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_representative-sequences.qza \
   --output-dir ../output/picrust2 \
   --p-placement-tool sepp \
   --p-threads 4 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose \
   --p-min-align 0.1
```

> > [!WARNING]
> >
> > Warning - 661 input sequences aligned poorly to reference sequences
> > (–min_align option specified a minimum proportion of 0.1 aligning to
> > reference sequences). These input sequences will not be placed and
> > will be excluded from downstream steps. 323 of 6521 ASVs were above
> > the max NSTI cut-off of 2.0 and were removed from the downstream
> > analyses.

The output artifacts of this command are the red boxes in the flowchart
[here](https://huttenhower.sph.harvard.edu/picrust/#:~:text=PICRUSt%202.0%20Flowchart).

These output files in (q2-picrust2_output) are:

- ec_metagenome.qza : EC metagenome predictions (rows are EC numbers and
  columns are samples).

- ko_metagenome.qza : KO metagenome predictions (rows are KOs and
  columns are samples).

- pathway_abundance.qza : MetaCyc pathway abundance predictions (rows
  are pathways and columns are samples).

The artifacts are all of type FeatureTable\[Frequency\], which means
they can be used with QIIME 2 plugins that process and analyze these
datatypes.

# Run with `mp` & `epa-ng`

> if users have at least 16 GB of RAM they should place reads with
> `epa-ng`

    --p-placement-tool TEXT Choices('epa-ng', 'sepp')
                           Placement tool to use when placing sequences into
                           reference tree. EPA-ng is the default, but SEPP is less
                           memory intensive.                   [default: 'epa-ng']

> **We recommend that in practice users use the `mp` method**

    --p-hsp-method TEXT Choices('mp', 'emp_prob', 'pic', 'scp',
        'subtree_average') Which hidden-state prediction method to use.
                                                                   [default: 'mp']

``` bash
LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH \
qiime picrust2 full-pipeline \
   --i-table ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_featureTable_filtered.qza \
   --i-seq ../salipante/241121_StonyCoral/270x200/250414_StonyCoral_270x200_representative-sequences.qza \
   --output-dir ../output/picrust2_epa-ng \
   --p-placement-tool epa-ng \
   --p-threads 20 \
   --p-hsp-method mp \
   --p-max-nsti 2 \
   --verbose \
   --p-min-align 0.1
```

> [!WARNING]
>
> Warning - 661 input sequences aligned poorly to reference sequences
> (–min_align option specified a minimum proportion of 0.1 aligning to
> reference sequences). These input sequences will not be placed and
> will be excluded from downstream steps. 287 of 6521 ASVs were above
> the max NSTI cut-off of 2.0 and were removed from the downstream
> analyses.
