---
title: "Fundulus Differential DNA Methylation Analysis"
author: "Mathia Colwell"
output: html_document
---
This document outlines our bioinformatics pipeline to call differential DNA Methylation from .BAM files generated from our long-read Oxford Nanopore Technology sequencing data and downstream analysis.


##### Call Modified bases

CpG methylation and hydroxymethyaltion global frequencies were estimated using Modkit v0.5.1. Modkits pileup command was run on the sorted BAM files to generate a per-site CpG summary.

### Load Libraries
```{r setup, message=FALSE, warning=FALSE}
library(methylKit)
library(GenomicRanges)
library(tidyverse)
```

### Set Up Sample Data
```{r}
sample_info <- tibble(
  sample_id = c(
    "PAH-sensitive_+PAH",
    "PAH-sensitive_control",
    "PAH-tolerant_+PAH",
    "PAH-tolerant_control"
  ),
  treatment = c(1, 0, 1, 0),
  file = c(
    "sensitive_PAH.m.bed",
    "sensitive_control.m.bed",
    "tolerant_PAH.m.bed",
    "tolerant_control.m.bed"
  )
)
```

### Import Methylation Data
```{r}
modkit_cols <- list(
  fraction = FALSE,
  chr.col = 1,
  start.col = 2,
  end.col = 3,
  coverage.col = 5,
  strand.col = 6,
  freqC.col = 11
)

meth_obj <- methRead(
  location = sample_info$file,
  sample.id = sample_info$sample_id,
  assembly = "Fundulus_heteroclitus",
  treatment = sample_info$treatment,
  context = "CpG",
  pipeline = modkit_cols,
  mincov = 10
)
```

## Tiling and Differnetial Methylation 
```{r}
tiles_united <- unite(tiles)

dmr_results <- calculateDiffMeth(tiles_united)

dmrs <- getMethylDiff(
  dmr_results,
  difference = 25,
  qvalue = 0.01
)

dmrs_df <- getData(dmrs)
```

### Summary statistics of DMRs 
```{r}
dmrs_df %>%
  summarize(
    total_dmrs = n(),
    hypermethylated = sum(meth.diff > 0),
    hypomethylated = sum(meth.diff < 0),
    mean_abs_diff = mean(abs(meth.diff))
  )
```

