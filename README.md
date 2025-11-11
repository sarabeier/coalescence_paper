Coalescence describes the mixing of entire communities, a process frequently encountered by naturally occurring microbial communities. In order to assess how organisms with differnt life histories affect the outcome of coalescence, we experimentally mixed communities from two aquatic sites with contrasting levels of eutrophication and disturbance histories and inspected compositional and functional community parameters after 5 days of incubation.

This repository holds all scripts used to generate the data and figures for the paper:

**Dominance of opportunistic generalist species after aquatic microbial community coalescence**
accepted in [Ecosphere](https://esajournals.onlinelibrary.wiley.com/journal/21508925)

## Table of contents
1. [script01_map.envdata.pdf](#script01_map.envdata.pdf)
2. [script02_community.dada.R](#script02_community.dada.R)
3. [script03_genomic_traits_estimation.md](#script03_genomic_traits_estimation.md)
4. [script04_community.patterns.pdf](#script04_community.patterns.pdf)
5. [script05_anova_stats.pdf](#script05_anova_stats.pdf)
6. [Stables.xlsx](#Stables.xlsx)

## script01_map.envdata.pdf
This script is written in R and contains the code to create Figure 1 of the manuscript (map of sample sites and boxplot of time series chlorophyll a and salinity data from the sites SOLA and Canet).

## script02_community.dada.R
This script is written in R and includes the code used to process metabarcoding sequence data via the [DADA2](https://benjjneb.github.io/dada2/index.html) pipline.

## script03_genomic_traits_estimation.md
This script contains code written in R and Bash. It utilizes genomic trait values from previously sequenced genomes ([Beier et al. 2022](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2022.985216/full)) to predict corresponding traits for closely related ASVs in the sequenced communities using the [PICRUSt2](https://huttenhower.sph.harvard.edu/picrust/) software. Predicted genomic traits—including genome size, the fraction of transcription factors per genome, the number of 16S rRNA gene copies, and maximal growth rates estimated via codon usage bias.

## script04_community.patterns.pdf
This script is written in R and includes the code to obtain and visualize community related patterns, such as community weightyed means for the estimated genomioc traits, or community structural statistics (PERMANOVA) and patterns based on community composition. This script further contains the code for Figure 3 (PCA displaying community characteristics of parent communitiues from SOLA and the Canet Lagoon), Fig. 4 (growth curves and barplots displaying community composition at the order level) and fig. 5 (PCoA biplots).

## script05_anova_stats.pdf
This script is written in R and contains code for two way and reapeated measurement ANOVAS of parameters displayed in the boxplot graphics of Figure 6.

## Stables.xlsx
This file includes contextual data published as supplementary material and is used as input for the supplied code.

