---
title: "Variant Analysis in Fundulus heteroclitus"
author: "Mathia Colwell"
output: html_document
---

## Variant calling
Single-nucleotide polymorphisms and structural variants were idenfieid from long read Oxford Nanopore
sequencing data using Clair3.  Larger insertions, deletions, and other structural variants 
were idnetified using Sniffles. Variants were called against the Kings
Creek *Fundulus heteroclitus* genome assembly (accession no. processing).


Variant call files (VCFs) were filtered for a quality score of ≥ 30, read depth of ≥ 10, 
and biallelic SNPs. Thresholds were chosen to reduce false positives.

## Load Libraries
```{r setup, message=FALSE, warning=FALSE}
library(tidyverse)
library(VariantAnnotation)
library(GenomicRanges)
```

## Sample of VCF file format:
```{r}
variants <- tibble(
  chr = sample(paste0("chr", 1:20), 5000, replace = TRUE),
  pos = sample(1:5e6, 5000),
  ref = sample(c("A", "C", "G", "T"), 5000, replace = TRUE),
  alt = sample(c("A", "C", "G", "T"), 5000, replace = TRUE),
  qual = runif(5000, 20, 60),
  depth = sample(10:100, 5000, replace = TRUE)
) %>%
  filter(ref != alt)
```

## Variant filtering 
```{r}
variants_filt <- variants %>%
  filter(
    qual >= 30,
    depth >= 10
  ) %>%
  mutate(
    variant_type = case_when(
      caller == "Clair3" & nchar(ref) == 1 & nchar(alt) == 1 ~ "SNV",
      caller == "Sniffles" ~ "Structural variant",
      TRUE ~ "Other"
    ),
    substitution = ifelse(
      variant_type == "SNV",
      paste0(ref, ">", alt),
      NA
    )
  )
```

## Identifiying C>T Mutations
```{r}
# Read SNPs
snps <- read_tsv("KC_snps.vcf.tsv")

# Identify C>T mutations
ct_snps <- snps %>%
  filter(REF == "C" & ALT == "T")
```
