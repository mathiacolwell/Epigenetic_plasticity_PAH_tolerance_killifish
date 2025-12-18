---
title: "Workflow for Fundulus Multi-omics"
author: "Mathia Colwell, PhD"
date: "2025-12-09"
output: html_document
---


## Description 
How we handled our long-read sequencing data to assess:
1. Align reads to a non-model organsim genome
*Call modified bases
2. Identify differentially methylated regions (DMRs)
3. Call SNPs, SV, and C>T mutations
* Explore C>T mutations in DMRs and CpG sites
4. Measure Differentially Expressed Genes (DEGs) using short-read illumina sequencing
* Align reads to a non-model organism genome
* Parse parent transcripts from novel transcripts
5. Identify cis-nats in the AHR gene-battery
6. Overlap DEGs with DMRs and Regions of Interest

## Align reads to the Kings Creek Assembly

Long-read unmapped sequences from samples (KC21Clean example) we re-aligned to the Kings Creek fundulus heteroclitus reference genome using Dorado v1.2.0. Dorado alignmnents were piped directly into SAMtools v1.19 and sorted by genomic coordinates. Resulting BAM files were indexed.

```{bash}
# Align and sort 
dorado aligner \
  kings_creek_genome_25.fa \
  KC21Clean.unmapped.bam \
  | samtools sort \
      -@ 32 \
      -o KC21Clean.mapped.sorted.bam
      
# Index files 
samtools index KC21Clean.mapped.sorted.bam 
```

##### Call Modified bases

CpG methylation and hydroxymethyaltion global frequencies were estimated using Modkit v0.5.1. Modkits pileup command was run on the sorted BAM files to generate a per-site CpG summary.

```{bash}
modkit pileup \
  -- ref kings_creek_genome_25.fa \
  --cpg \
  KC21Clean.mapped.sorted.bam \
  KC21Clean.all.mods.bed
```

The resulting BED files contain the total number of calls and counts of the canonical cytosine, 5-methylcytosine (5mC) and 5-hydroxymethylcytosine (5hmC). For detailed column specs, refer to Modkit GitHub.

```{bash}
for bed in *.bed; do
    awk -v file="$bed" '
        $4 == "m" {
            can  += $13;   # canonical (unmodified C)
            mod  += $12;   # methylated (5mC)
            oth  += $14;   # hydroxymethylated (5hmC or other)
            valid+= $10;   # total observations
        }
        END {
            printf "%s\tCpG canonical %.6f\tCpG methyl %.6f\tCpG hydroxy %.6f\n",
                file, (can/valid), (mod/valid), (oth/valid)
        }
    ' "$bed"
done >> methylation-summary.txt
```

Next, we filtered the modkit output to contain only methylated CpGs

```{bash}
# 5-methylcytosine
awk '$4=="m"' KCClean.all.mods.bed > KCClean.m.bed

# 5-hydroxymethylatosine
awk '$4=="h"' KCClean.all.mods.bed > KCClean.h.bed
```

## Call DMRs using a 25% Cutoff

```{r}
# Load Libraries
library(methylKit)
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(genomation)
library(readr)
library(rtracklayer)

# Load files 
file.list <- list(
"KC21ERSE.m.bed",
"KC21Clean.m.bed",
"RP21ERSE.m.bed",
"RP21Clean.m.bed"
)

# Set Modkit Column Names
modkit_cols <- list(
fraction = FALSE,
chr.col = 1,
start.col = 2,
end.col = 3,
coverage.col = 5,
strand.col = 6,
freqC.col = 11
)

# Read the files with no minimum coverage
myobj_lowCov <- methRead(file.list,
                         sample.id = list("PAH-sensitive_+PAH", 
                                          "PAH-sensitive_control", 
                                          "PAH-tolerant_+PAH", 
                                          "PAH-tolerant_control"),
                         assembly = "fundulus",
                         treatment = c(1, 0, 1, 0),
                         context = "CpG",
                         pipeline = modkit_cols
)

# Use 1000 bp tiled windows 
tiles <- tileMethylCounts(
  myobj_lowCov,
  win.size  = 1000,
  step.size = 1000,
  cov.bases = 10
)

# Run on all comparisons we are interested in 
process_comparison <- function(tiles, sample1, sample2) {

  # Subset methylation tiles for pairwise comparison
  tiles_subset <- reorganize(
    tiles,
    sample.ids = c(sample1, sample2),
    treatment  = c(1, 0)
  )

  # Unite tiles into methylBase object
  meth.tiles <- methylKit::unite(tiles_subset)

  # Differential methylation
  myDiff <- calculateDiffMeth(meth.tiles)

  # Filter: ≥ 25% difference, q ≤ 0.01
  dmr <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01)

  # Save DMR table
  out_file <- file.path(
    "/path",
    paste0(sample1, "vs", sample2, "_DMRs.txt")
  )
  write.table(dmr, file = out_file, sep = "\t",
              row.names = FALSE, quote = FALSE)

  list(dmrs = dmr)
}

comparisons <- list(
  c("PAH-sensitive_+PAH", "PAH-sensitive_control"),
  c("PAH-tolerant_control", "PAH-sensitive_control"),
  c("PAH-tolerant_+PAH", "PAH-sensitive_+PAH"),
  c("PAH-tolerant_+PAH", "PAH-tolerant_control")
)

esults <- lapply(comparisons, function(pair) {

  result <- process_comparison(tiles, pair[1], pair[2])

  cat("\n----", pair[1], "vs", pair[2], "----\n")
  cat("Total DMRs:", nrow(result$dmrs), "\n")
  cat("Hypermethylated:", sum(result$dmrs$meth.diff > 0), "\n")
  cat("Hypomethylated:", sum(result$dmrs$meth.diff < 0), "\n")
  cat("Mean |delta methylation|:",
      mean(abs(result$dmrs$meth.diff)), "\n")

  return(result)
})

names(results) <- sapply(comparisons, paste, collapse = "vs")
```

Create an upset plot of the comparisons

```{r}
library(UpSetR)
library(dplyr)

# Read DMR data for each comparison
read_dmr_data <- function(filename) {
  read.table(filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

dmr_files <- list(
  "/PAH-tolerant_+PAHvsPAH-sensitive_+PAH_DMRs.txt",
   "/PAH-tolerant_controlvsPAH-sensitive_control_DMRs.txt",
   "/PAH-tolerant_+PAHvsPAH-tolerant_control_DMRs.txt", 
  "/PAH-sensitive_+PAHvsPAH-sensitive_control_DMRs.txt"
)

dmr_data <- lapply(dmr_files, read_dmr_data)

# Create a list of DMR locations for each comparison
dmr_locations <- lapply(dmr_data, function(x) paste(x$chr, x$start, x$end, sep = "_"))

# Create new name list for UpSetR
names(dmr_locations) <- c(
  "PAH tolerant +PAH vs. PAH-sensitive +PAH",
  "PAH-tolerant control vs. PAH-sensitive control",
  "PAH-tolerant +PAH vs. PAH-tolerant control",
  "PAH-sensitive +PAH vs. PAH-sensitive control"
)

# Desired order
set_order <- c(
  "PAH tolerant +PAH vs. PAH-sensitive +PAH",
  "PAH-tolerant control vs. PAH-sensitive control",
  "PAH-tolerant +PAH vs. PAH-tolerant control",
  "PAH-sensitive +PAH vs. PAH-sensitive control"
)

# Create the UpSet plot
upset(fromList(dmr_locations), 
      nsets = 4, 
      order.by = "freq",
      sets = set_order,  
      keep.order = TRUE,
      main.bar.color = "#0000FF", 
      sets.bar.color = "#FF0000",
      matrix.color = "#008000",
      text.scale = c(1.2, 0.8, 0.8, 0.8, 1.2, 1),
      point.size = 2,
      line.size = 0.5,
      mainbar.y.label = "DMR Intersections",
      sets.x.label = "DMRs per Comparison",
      mb.ratio = c(0.55, 0.45))
# Save the plot
ggsave("/DMR_UpSet_plot.svg",
       width = 8,
       height = 9)
```

Add gene names to the DMRs, then identify closest genes and distance away from closest gene. This pipeline is specific to the **Kings Creek Fundulus heteroclitus** GTF, not the 4.1 Missouri Reference Genome.

```{r}
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(readr)

# Specific to our GTF
process_liftoff_annotations <- function(gtf_path) {
  # Import GTF file
  gtf <- import(gtf_path, format = "gtf")
  
  # Extract only transcript features 
  transcripts <- gtf[gtf$type == "transcript", ]
  
  # Create standardized gene annotations
  gene_annotations <- data.frame(
    chr = as.character(seqnames(transcripts)),
    start = start(transcripts),
    end = end(transcripts),
    strand = as.character(strand(transcripts)),
    gene_name = tolower(transcripts$gene_name),  
    stringsAsFactors = FALSE
  ) %>% 
    distinct(gene_name, .keep_all = TRUE)
  
  # Convert to GRanges
  GRanges(
    seqnames = gene_annotations$chr,
    ranges = IRanges(start = gene_annotations$start, 
                    end = gene_annotations$end),
    strand = gene_annotations$strand,
    gene_name = gene_annotations$gene_name
  )
}

# Annotate the DMRs with the gene name only 
annotate_dmrs_with_genename <- function(dmr_file, gene_annotations) {
  # Read DMR data
  dmrs <- read_tsv(dmr_file, show_col_types = FALSE) %>%
    mutate(chr = as.character(chr))
  
  # Convert to GRanges
  dmr_gr <- makeGRangesFromDataFrame(dmrs, 
                                   keep.extra.columns = TRUE,
                                   seqnames.field = "chr",
                                   start.field = "start",
                                   end.field = "end")
  
  # Find overlapping genes
  overlaps <- findOverlaps(dmr_gr, gene_annotations)
  dmrs$gene_name <- NA_character_
  dmrs$gene_name[queryHits(overlaps)] <- gene_annotations$gene_name[subjectHits(overlaps)]
  
  # For non-overlapping DMRs, find nearest gene
  non_overlapping <- which(is.na(dmrs$gene_name))
  if(length(non_overlapping) > 0) {
    nearest <- distanceToNearest(dmr_gr[non_overlapping], gene_annotations)
    dmrs$nearest_gene <- NA_character_
    dmrs$distance_to_gene <- NA_integer_
    dmrs$nearest_gene[non_overlapping] <- gene_annotations$gene_name[subjectHits(nearest)]
    dmrs$distance_to_gene[non_overlapping] <- mcols(nearest)$distance
  }
  
  # Fill nearest_gene with gene_name when gene_name exists but nearest_gene is NA
  dmrs <- dmrs %>%
    mutate(nearest_gene = ifelse(!is.na(gene_name) & is.na(nearest_gene), 
                                gene_name, 
                                nearest_gene),
           distance_to_gene = ifelse(!is.na(gene_name) & is.na(distance_to_gene),
                                   0,
                                   distance_to_gene))
  
  # Select and order output columns
  dmrs %>%
    select(chr, start, end, strand, pvalue, qvalue, meth.diff,
           gene_name, nearest_gene, distance_to_gene)
}

# Process annotations
gene_annotations <- process_liftoff_annotations(
  "/path/liftoff.gtf.gz"
)

# List of DMR files
dmr_files <- c(
  "PAH-sensitive_+PAHvsPAH-sensitive_control_DMRs.txt",
  "PAH-tolerant_+PAHvsPAH-sensitive_+PAH_DMRs.txt",
  "PAH-tolerant_controlvsPAH-sensitive_control_DMRs.txt",
  "PAH-tolerant_+PAHvsPAH-tolerant_control_DMRs.txt"
)

# Process all DMR files
for(file in dmr_files) {
  annotated_dmrs <- annotate_dmrs_with_genename(file, gene_annotations)
  
  output_file <- gsub("_DMRs.txt", "_annotated_genename.tsv", file)
  write_tsv(annotated_dmrs, output_file)
  message("Saved ", output_file, " with ", 
          sum(!is.na(annotated_dmrs$gene_name)), 
          " direct gene overlaps")
}
```

We were interested in detecting region specific methylation. The region input files were .bed files with the gene coordinates.

```{bash}

# Set the path to the  ROI file
ROI_FILE="/path/gene_coordinates.tsv"

# Create an array of .bed files
BED_FILES=("/path/"/*.m.bed)

# Create output file and write header
OUTPUT_FILE="methylation_results.txt"
echo -e "Gene\tChromosome\tStart\tEnd\t$(basename -a "${BED_FILES[@]}" | sed 's/\.m\.bed$/_Methylation/' | tr '\n' '\t')$(basename -a "${BED_FILES[@]}" | sed 's/\.m\.bed$/_CpG_Count/' | tr '\n' '\t' | sed 's/\t$//')" > "$OUTPUT_FILE"

# Loop through each line in the ROI file
while IFS=$'\t' read -r chr start end gene; do
    echo "Processing region: $gene ($chr:$start-$end)"
    results=("$gene" "$chr" "$start" "$end")

    # Process each .bed file
    for bed_file in "${BED_FILES[@]}"; do
        echo "Analyzing $(basename "$bed_file")..."
        result=$(awk -v chr="$chr" -v start="$start" -v end="$end" '
            $1 == chr && $2 >= start && $3 <= end && $4 == "m" {
                mod += $12
                valid += $10
                count++
            }
            END {
                if (valid > 0) {
                    printf "%.6f\t%d", mod/valid, count
                } else {
                    print "NA\t0"
                }
            }
        ' "$bed_file")
        results+=($result)
    done

    # Write results to output file
    echo -e "$(IFS=$'\t'; echo "${results[*]}")" >> "$OUTPUT_FILE"
done < "$ROI_FILE"
```

Below is the DMR methyaltion specific region. Example below is voltage gated K channels (vkg for short)

```{bash}
#!/bin/bash

#Pull DMR methylation based off gene name

RPCvKCC="/DMRs_PAH-tolerant_controlvsPAH-sensitive_control_annotated_genename.no.header.tsv"
RPEvRPC="/DMRs_PAH-tolerant_+PAHvsPAH-tolerant_control_annotated_genename.no.header.tsv"
RPEvKCE="/DMRs_PAH-tolerant_+PAHvsPAH-sensitive_+PAH_annotated_genename.no.header.tsv"
KCEvKCC="/DMRs_PAH-sensitive_+PAHvsPAH-sensitive_control_annotated_genename.no.header.tsv"

ROI="/path/ROI_coordinates.txt"

mkdir -p vgk_intersects

# Function to run intersection
extract_dmrs () {
    local file=$1
    local label=$2

    bedtools intersect \
        -a "$ROI" \
        -b "$file" \
        -wa -wb > "vgk_intersects/${label}.raw.tsv"

    # Extract useful columns:
    awk -v L="$label" '
        {
          roi_chr=$1; roi_start=$2; roi_end=$3; gene=$4;
          dmr_chr=$5; dmr_start=$6; dmr_end=$7;
          meth=$11; nearest=$13;

          print gene, roi_chr, roi_start, roi_end, dmr_chr, dmr_start, dmr_end, meth, nearest, L;
        }
    ' OFS="\t" "vgk_intersects/${label}.raw.tsv" \
      > "vgk_intersects/${label}.clean.tsv"
}

extract_dmrs "$RPCvKCC" "RPCvKCC"
extract_dmrs "$RPEvRPC" "RPEvRPC"
extract_dmrs "$RPEvKCE" "RPEvKCE"
extract_dmrs "$KCEvKCC" "KCEvKCC"
```

#### Create Methylartist Graphs

```{bash}
#!/bin/bash

# Set variables
ROI_FILE="/region.txt"
REF_GENOME="/Kings_creek_genome_25.fa"
GTF_FILE="/liftoff3.gtf.gz"
LOCUS_PLOT_DIR="/output_directory/"
RESULTS_FILE="/Methylartist_loop_methylation_results.txt"

# Create header for results file
echo -e "Gene\tChromosome\tStart\tEnd\tAverage_Methylation_Sensitive\tAverage_Methylation_Tolerant" > "$RESULTS_FILE"

# Process each region in the ROI file
while IFS=$'\t' read -r chr start end gene; do
    # Skip header or invalid lines
    if [[ "$chr" == "chr" || -z "$chr" || -z "$start" || -z "$end" || -z "$gene" ]]; then
        echo "$(date): Skipping invalid entry: $chr $start $end $gene"
        continue
    fi

    # Adjust start and end positions to include 5000 bps flanks, this will help make the graph look nice. 
    flanked_start=$((start - 5000))
    flanked_end=$((end + 5000))

    # Ensure flanked_start is not less than 1
    if (( flanked_start < 1 )); then
        flanked_start=1
    fi

    # Define the interval with flanks and the highlight region (the actual ROI)
    INTERVAL="${chr}:${flanked_start}-${flanked_end}"
    HIGHLIGHT="${chr}:${start}-${end}"
    echo "$(date): Processing region: $gene ($INTERVAL) with highlight: $HIGHLIGHT"

    # Generate locus plot for Pair 1 (PAH-sensitive)
    LOCUS_OUTPUT1="$LOCUS_PLOT_DIR/${gene}_sensitive_locus.svg"
    echo "$(date): Generating locus plot for $gene (PAH-sensitive)..."
    methylartist locus \
        -b KC21ERSE_tiag.q10.bam,KC21Clean_tiag.q10.bam \
        -i "$INTERVAL" \
        -l "$HIGHLIGHT" \
        --ref "$REF_GENOME" \
        --gtf "$GTF_FILE" \
        --motif CG \
        --primary_only \
        --labelgenes \
        --hidelegend \
        --mods m \
        --svg \
        -o "$LOCUS_OUTPUT1"

    # Check if locus plot generation was successful; failure to generate plot is usually caused by a problem in the GTF
    if [ ! -f "$LOCUS_OUTPUT1" ]; then
        echo "$(date): Error: Failed to generate locus plot for $gene (PAH-sensitive)"
        continue
    fi

    # Generate locus plot for Pair 2 (PAH-tolerant; I changed the colors for this after the .bam files for easy comparison between the tolerant and sensitive samples)
    LOCUS_OUTPUT2="$LOCUS_PLOT_DIR/${gene}_tolerant_locus.svg"
    echo "$(date): Generating locus plot for $gene (PAH-tolerant)..."
    methylartist locus \
        -b RP21ERSE_tiag.q10.bam:#d62728,RP21Clean_tia.q10.bam:#2ca02c \
        -i "$INTERVAL" \
        -l "$HIGHLIGHT" \
        --ref "$REF_GENOME" \
        --gtf "$GTF_FILE" \
        --motif CG \
        --labelgenes \
        --hidelegend \
        --primary_only \
        --mods m \
        --svg \
        -o "$LOCUS_OUTPUT2"

    # Check if locus plot generation was successful
    if [ ! -f "$LOCUS_OUTPUT2" ]; then
        echo "$(date): Error: Failed to generate locus plot for $gene (PAH-tolerant)"
        continue
    fi

    # Calculate average methylation for the region
    echo "Calculating average methylation for $gene..."
    AVG_METH_SENSITIVE=$(methylartist segmeth \
        -b KC21ERSE_tiag.q10.sorted.bam,KC21Clean_tiag.q10.sorted.bam \
        -i "$INTERVAL" \
        --ref "$REF_GENOME" \
        --motif CG \
        --primary_only \
        --mods m \
        --output /dev/stdout | \
        awk 'NR>1 {sum+=$5; count++} END {print (count>0) ? sum/count : "NA"}')

    AVG_METH_TOLERANT=$(methylartist segmeth \
        -b RP21ERSE.modmapped.q10.sorted.bam,RP21Clean_tiag.q10.sorted.bam \
        -i "$INTERVAL" \
        --ref "$REF_GENOME" \
        --motif CG \
        --primary_only \
        --mods m \
        --output /dev/stdout | \
        awk 'NR>1 {sum+=$5; count++} END {print (count>0) ? sum/count : "NA"}')

    # Results
    echo -e "$gene\t$chr\t$start\t$end\t$AVG_METH_SENSITIVE\t$AVG_METH_TOLERANT" >> "$RESULTS_FILE"

    echo "Completed processing $gene"
done < "$ROI_FILE"
```
Create a heatmap using region of interest DNA Methyaltion values

```{r}
# Read in methylation data
methylation_data <- read_delim(
  "/gene_methylation.txt",
  delim = "\t",
  col_names = TRUE
)

# Rename the columns (samples)
colnames(methylation_data)[5:8] <- c(
  "PAH-Sensitive Control",
  "PAH-Sensitive PAH Exposed",
  "PAH-Tolerant Control",
  "PAH-Tolerant Exposed"
)

# Long format
methylation_long <- methylation_data %>%
  pivot_longer(
    cols = 5:8,
    names_to = "sample",
    values_to = "methylation"
  ) %>%
  mutate(methylation = as.numeric(methylation))

# Average by gene + sample
methylation_aggregated <- methylation_long %>%
  group_by(Gene, sample) %>%
  summarize(mean_methylation = mean(methylation, na.rm = TRUE), .groups = "drop")

# Wide matrix
methylation_matrix <- methylation_aggregated %>%
  pivot_wider(
    names_from = sample,
    values_from = mean_methylation,
    values_fill = list(mean_methylation = 0)
  ) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Reorder samples
sample_order <- c(
  "PAH-Sensitive Control",
  "PAH-Sensitive PAH Exposed",
  "PAH-Tolerant Control",
  "PAH-Tolerant Exposed"
)
methylation_matrix <- methylation_matrix[, sample_order]

# Sort genes 
gene_order <- order(
  methylation_matrix[, "PAH-Sensitive Control"],
  decreasing = TRUE
)
methylation_sorted <- methylation_matrix[gene_order, ]

# Compute relative z-scores 
control_mean <- mean(
  methylation_sorted[, "PAH-Sensitive Control"],
  na.rm = TRUE
)
control_sd <- sd(
  methylation_sorted[, "PAH-Sensitive Control"],
  na.rm = TRUE
)

relative_z <- (methylation_sorted - control_mean) / control_sd

relative_z[relative_z >  3] <-  3
relative_z[relative_z < -3] <- -3

metab_wrapped_labels <- c(
  "PAH-Sensitive\nControl",
  "PAH-Sensitive\nPAH Exposed",
  "PAH-Tolerant\nControl",
  "PAH-Tolerant\nExposed"
)

heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(256)
heatmap_breaks <- seq(-3, 3, length.out = 257)

heatmap <- pheatmap(
  metab_relative_z,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = metab_heatmap_colors,
  breaks = metab_heatmap_breaks,
  fontsize = 16,
  fontsize_row = 8,
  fontsize_col = 12,
  angle_col = 0,
  labels_col = metab_wrapped_labels,
  family = "Arial",
  main = "Gene Methylation",
  na_col = "white",
  border_color = "grey"
)

# Save as SVG
save_pheatmap_svg(
  heatmap,
  "/methylation.heatmap.scaled.svg"
)
```

## Call Differentially Expressed Genes

Type of Sequencing: Illumina Short Read

```{bash}
#!/bin/bash

# Reference files
GENOME_FASTA="/Kings_Creek_genome_25.fa"
ANNOTATION_GTF="/liftoff3.gtf"
RAW_DATA_DIR="/directory_to_raw_RNA_seq_data/"
OUTPUT_DIR="/star_results"
STAR_INDEX_DIR="/star_results/star_index"


THREADS=32
SAMTOOLS_THREADS=32

# Create output directories
mkdir -p ${OUTPUT_DIR}/{alignments,assemblies,counts_}

# Generate genome index 
if [ ! -d "${STAR_INDEX_DIR}/SAindex" ]; then
  echo "Generating STAR index..."
  mkdir -p ${STAR_INDEX_DIR}

  STAR --runThreadN ${THREADS} \
    --runMode genomeGenerate \
    --genomeDir ${STAR_INDEX_DIR} \
    --genomeFastaFiles ${GENOME_FASTA} \
    --sjdbGTFfile ${ANNOTATION_GTF} \
    --sjdbOverhang 149 \
    --genomeSAindexNbases 13 \
    --genomeChrBinNbits 12 \
    --limitGenomeGenerateRAM 30000000000 \
    --sjdbGTFfeatureExon exon \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFtagExonParentGene gene_id
fi

#Temp dir setup
mkdir -p /media/tia/athena_storage/tmp
export TMPDIR=$(mktemp -d -p /tmp star_temp.XXXXXX)
chmod 700 "$TMPDIR"

ulimit -n 10000

# Alignment loop
for R1 in ${RAW_DATA_DIR}/*_R1_001.fastq.gz; do
  sample=$(basename ${R1} _R1_001.fastq.gz)
  R2=${RAW_DATA_DIR}/${sample}_R2_001.fastq.gz

  echo "Processing ${sample}"
 
  # Post-Alignment
if [[ -s "${OUTPUT_DIR}/alignments/${sample}_Aligned.sortedByCoord.out.bam" ]]; then
  mv "${OUTPUT_DIR}/alignments/${sample}_Aligned.sortedByCoord.out.bam" \
     "${OUTPUT_DIR}/alignments/${sample}.bam"
  samtools index -@ ${SAMTOOLS_THREADS} "${OUTPUT_DIR}/alignments/${sample}_v3.bam"
  samtools quickcheck "${OUTPUT_DIR}/alignments/${sample}_v3.bam" || exit 1
  grep 'strandInfo' "${OUTPUT_DIR}/alignments/${sample}_Log.final.out" >> "${OUTPUT_DIR}/strand_report.txt"
else
  echo "ERROR: STAR failed for ${sample}"
  exit 1
fi

  # Transcript assembly
  stringtie \
    "${OUTPUT_DIR}/alignments/${sample}.bam" \
     -G $GTF_REFERENCE \
      -e -B -eB -p 32 \
      -o ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.gtf
      -p ${THREADS}

done

# Cleanup the temp
rm -rf "$TMPDIR"
```

Separate cis-nats and parent transcripts from main file


**Based on transcript naming patterns used in StringTie/STAR output**

A cis-NAT is defined as:
* transcript ID begins with "STRG."
* gene_name = transcript_id
The Parent transcripts are all remaining features.

# Filter for Novel Transcripts

```{r}
# Load data
abundance <- read.delim("/KCClean_output.tsv", header = TRUE, sep = "\t")

# Calculate transcript length
abundance$Length <- abundance$End - abundance$Start + 1

# Print starting count
cat("Total input transcripts:", nrow(abundance), "\n")

# Filter length >= 200
filtered_length <- subset(abundance, Length >= 200)

# Filter TPM >= 1
filtered_tpm <- subset(filtered_length, TPM >= 1)

# Filter Coverage >= 10
filtered_data <- subset(filtered_tpm, Coverage >= 10)

# Identify novel transcripts (Gene ID starts with "STRG.")
novel_transcripts <- subset(filtered_data, grepl("^STRG\\.", Gene.ID))

# Output summary
cat("Transcripts after length >= 200 filter:", nrow(filtered_length), "\n")
cat("Transcripts after TPM >= 1 filter:", nrow(filtered_tpm), "\n")
cat("Transcripts after Coverage >= 10 filter:", nrow(filtered_data), "\n")
cat("Remaining novel transcripts (STRG.):", nrow(novel_transcripts), "\n")

# Save filtered tables
write.table(filtered_data, file = "/filtered_abundance.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(novel_transcripts, file = "/novel_transcripts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

```

#### Variant Calling

```{bash}
#!/bin/bash

REFERENCE_GENOME="/Kings_Creek_genome_25.fa"
OUTPUT_DIR="/all_bcfs"
MODEL_PATH="/bin/models/r941_prom_sup_g5014"
THREADS=32


# Location of my BAM files
BAM_FILES=(
    "KC21Clean.sorted.bam"
    "KC21ERSE.sorted.bam"
    "RP21Clean.sorted.bam"
    "RP21ERSE.sorted.bam"
    )

# Function to process single BAM file
process_bam(){
    local bam_file=$1
    local sample_name=$(basename "$bam_file" .bam)
    local sample_output_dir="${OUTPUT_DIR}/${sample_name}"

    echo "Processing $bam_file..."
    
    mkdir -p "$sample_output_dir"

# Use the absolute path to run clair3
    /aboslute_path/run_clair3.sh --bam_fn="$bam_file" \
    --ref_fn="$REFERENCE_GENOME" \
    --output="$sample_output_dir" \
    --threads="$THREADS" \
    --platform="ont" \
    --model_path="$MODEL_PATH" \
    --include_all_ctgs \
    --use_whatshap_for_final_output_phasing

    echo "Finished processing $sample_name"
    }

# To process each BAM file:
for bam_file in "${BAM_FILES[@]}"; do
    process_bam "$bam_file"
    done

    echo "All samples are processed!"
```

Remove SNP overlap of KCClean with other samples

```{bash}
#!/bin/bash

# Set the path to your VCF files
KC21CLEAN_VCF="/KC21Clean.merged.filtered.vcf.gz"
INPUT_DIR="/input_dir"
OUTPUT_DIR="/output_dir"

# Ensure the KC21CLEAN VCF file is indexed
bcftools index -f $KC21CLEAN_VCF

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Process all files ending with filtered.vcf.gz
for INPUT_VCF in ${INPUT_DIR}/*.merged.filtered.vcf.gz
do
    # Skip KC21Clean.merged.filtered.vcf.gz
    if [[ $INPUT_VCF == *"KC21Clean.merged.filtered.vcf.gz" ]]; then
        continue
    fi

    # Extract the base name of the input file
    BASE_NAME=$(basename $INPUT_VCF .merged.filtered.vcf.gz)

    # Set the output VCF name
    OUTPUT_VCF="${OUTPUT_DIR}/${BASE_NAME}_non_overlapping.vcf.gz"

    echo "Processing $INPUT_VCF..."

    # Ensure the input VCF file is indexed
    bcftools index -f $INPUT_VCF

    # Create a temporary directory
    TEMP_DIR=$(mktemp -d)

    # Use bcftools isec to find non-overlapping sites
    bcftools isec -C -c none -p $TEMP_DIR $INPUT_VCF $KC21CLEAN_VCF

    # Check if the isec command was successful
    if [ -f "$TEMP_DIR/0000.vcf" ]; then
        # Compress and index the resulting VCF
        bgzip $TEMP_DIR/0000.vcf
        bcftools index $TEMP_DIR/0000.vcf.gz

        # Move the final VCF to the desired output location
        mv $TEMP_DIR/0000.vcf.gz $OUTPUT_VCF
        mv $TEMP_DIR/0000.vcf.gz.csi ${OUTPUT_VCF}.csi

        echo "Non-overlapping sites for $BASE_NAME have been written to $OUTPUT_VCF"
    else
        echo "Error: bcftools isec failed for $INPUT_VCF"
    fi

    # Clean up temporary directory
    rm -rf $TEMP_DIR
done

echo "All files have been processed."
```

Count mutations

```{bash}
#!/bin/bash

# Set the variables

INPUT_DIR="/input_dir"
OUTPUT_DIR="/output_dir"
SUMMARY_FILE="summary_mutation_counts.txt"
THREADS=32

mkdir -p "OUTPUT_DIR"

# Function to process a single BCF file
process_vcf() {
	local vcf_file=$1
	local base_name=$(basename "$vcf_file" overlapping.vcf.gz)
	local output_file="OUTPUT_DIR/${base_name}_mutation_counts.txt"

	echo "Processing $vcf_file..."


	# bcftools to extract SNPs, count mutations, and calculate percent
	bcftools view -v snps "$vcf_file" | \
	bcftools query -f '%REF\t%ALT\n' | \
	awk '
    BEGIN {
        OFS="\t"
        print "Mutation_Type", "Count", "Percentage"
    }
    {
        mutation = $1 ">" $2
        count[mutation]++
        total++
    }
    END {
        for (mut in count) {
            percentage = (count[mut] / total) * 100
            printf "%s\t%d\t%.2f%%\n", mut, count[mut], percentage
        }
        print "Total_SNPs", total, "100%"
    }
    ' | sort -k2,2nr > "$output_file"

    echo "Results saved to $output_file"
}

export -f process_vcf

# Find all BCF files in the input directory
VCF_FILES=("$INPUT_DIR"/*overlapping.vcf.gz)

# Check if BCF files are found
if [ ${#VCF_FILES[@]} -eq 0 ]; then
    echo "No BCF files found in $INPUT_DIR"
    exit 1
fi

# Use GNU Parallel to process BCF files in parallel
parallel -j "$THREADS" process_vcf ::: "${VCF_FILES[@]}"

# Combine all individual results into a summary file
echo "Combining results..."
echo -e "File\tMutation_Type\tCount\tPercentage" > "$OUTPUT_DIR/$SUMMARY_FILE"
for vcf_file in "${VCF_FILES[@]}"; do
    base_name=$(basename "$Vcf_file" overlapping.vcf.gz)
    result_file="$OUTPUT_DIR/${base_name}_mutation_counts.txt"
    if [ -f "$result_file" ]; then
        awk -v file="$vcf_file" 'NR>1 {print file "\t" $0}' "$result_file" >> "$OUTPUT_DIR/$SUMMARY_FILE"
    fi
done

echo "Summary saved to $OUTPUT_DIR/$SUMMARY_FILE"
```

Find all C\>T Mutations in CpG Sites

```{bash}
#!/bin/bash

for vcf_file in *non_overlapping.vcf.gz; do
  base_name=$(basename "$vcf_file" non_overlapping.vcf.gz)

  # Count total lines and C>T mutations in CpG sites
  total_lines=$(zcat "$vcf_file" | wc -l)
  ct_cpg_mutations=$(zcat "$vcf_file" | awk '!/^#/ && ($4 == "C" && $5 == "T") && substr($0, $2, 2) ~ /CG/' | wc -l)

  echo "Total lines in $vcf_file: $total_lines"
  echo "C>T mutations in CpG sites found: $ct_cpg_mutations"

  # Filter C>T mutations in CpG sites and create a new VCF file (without header)
  zcat "$vcf_file" | awk '!/^#/ && ($4 == "C" && $5 == "T") && substr($0, $2, 2) ~ /CG/' > "${base_name}CT_CpG.vcf"

  # Create a BED file with genomic coordinates and sort it
  awk -v OFS='\t' '{print $1, $2-1, $2}' "${base_name}CT_CpG.vcf" | sort -k1,1 -k2,2n > "${base_name}CT_CpG.bed"

  # Check the size of the resulting BED file
  bed_size=$(wc -l < "${base_name}CT_CpG.bed")
  echo "Lines in ${base_name}CT_CpG.bed: $bed_size"

  echo "Processed $vcf_file"
done

echo "All files processed successfully"
```

Call Variants

```{bash}
#### Generate VCF and call SVs using Sniffles 
Sniffles2, introduced recently, is a redesign of the Sniffles pipeline for higher accuracy and time efficiency, and it also extends to population-scale SV calling. In the clustering step, Sniffles2 implements a three-phase clustering. It first clusters the raw SV signals by their SV type and genome start location, then corrects alignment errors in highly repetitive regions, and finally splits the preliminary clusters to represent different supported SV lengths. 

#Use KCClean as reference file 
samtools consensus KCClean.bam > KCClean_consensus.fasta

#Prepare and align  with other samples
samtools bam2fq sample.bam > sample.fastq

#Run the minimap again and prepare for Sniffles
minimap2 -ax map-ont -y KCClean_consensus.fasta sample.fastq > sample_aligned_to_KCClean.sam

samtools view -bS sample_aligned_to_KCClean.sam | samtools sort -o sample_aligned_to_KCClean.bam

samtools index sample_aligned_to_KCClean.bam

# Sort BAM file
samtools sort -o aligned_reads.sorted.bam aligned_reads.bam

# Index sorted BAM file
samtools index aligned_reads.sorted.bam

#Run command
sniffles --reference reference.fasta --input sample1.bam --vcf sample.vcf 
```

#### Specific Overlaps

##### DMRs vs. DEGs
```{r}
library(dplyr)
# Start with KCE v KCC for example
KCEvKCC_deg <- read.table("/DEGs_KCEvKCC.tsv", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
KCEvKCC_dmr <- read.table("/DMRs_PAH-sensitive_+PAHvsPAH-sensitive_control.tsv", header = TRUE, sep ='\t', stringsAsFactors = FALSE)

# merge based on 'gene_name'
 merged <- inner_join(KCEvKCC_deg, KCEvKCC_dmr, by = "gene_name")
 
dmr_summary <- merged %>%
  group_by(gene_name) %>%  
  summarise(
    chr = dplyr::first(chr),  
    start = min(start),
    end = max(end),
    LOGF2_mean = mean(as.numeric(LOGF2), na.rm = TRUE),
    meth_diff_mean = mean(as.numeric(meth.diff), na.rm = TRUE),
    dmr_pvalues = mean(as.numeric(pvalue), na.rm = TRUE),  
    dmr_qvalues = mean(as.numeric(qvalue), na.rm = TRUE),  
    .groups = "drop"
  ) %>%
  # Rename/reorder columns to match your desired output
  select(chr, start, end, gene_name, LOGF2_mean, meth_diff_mean, dmr_pvalues, dmr_qvalues)

# Write output
write.table(dmr_summary, "/KCEvKCC_DEG_DMR_overlap.tsv", sep='\t', quote= FALSE)
```

##### DMRs vs. SNPs

##### DMRs vs. Variants 
```{r}
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(purrr)

dmr_files <- c(
  "/DMRs_PAH-sensitive_+PAHvsPAH-sensitive_control_annotated_genename.tsv",
  "/DMRs_PAH-tolerant_+PAHvsPAH-sensitive_+PAH_annotated_genename.tsv",
  "/DMRs_PAH-tolerant_+PAHvsPAH-tolerant_control_annotated_genename.tsv",
  "/DMRs_PAH-tolerant_controlvsPAH-sensitive_control_annotated_genename.tsv"
)

sv_files <- c(
  "/Sniffles_KCClean.bed",
  "/Sniffles_KCERSE.bed",
  "/Sniffles_RPClean.bed",
  "/Sniffles_RPERSE.bed"
)

# Names for output
dmr_names <- gsub("\\.tsv$", "", basename(dmr_files))
sv_names  <- gsub("\\.bed$", "", basename(sv_files))

load_dmr <- function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t")
  
  GRanges(
    seqnames = df[[1]],
    ranges   = IRanges(start = df[[2]], end = df[[3]]),
    meth_diff = df[[7]],
    gene_name = df[[8]]
  )
}

load_sv <- function(file) {
  df <- read.delim(file, header = FALSE, sep = "\t")
  
  GRanges(
    seqnames = df[[1]],
    ranges   = IRanges(start = df[[2]], end = df[[3]]),
    svtype   = df[[4]]
  )
}

# Load all files
dmr_gr <- setNames(lapply(dmr_files, load_dmr), dmr_names)
sv_gr  <- setNames(lapply(sv_files,  load_sv), sv_names)

# Total Overlaps: Every SV/DMR overlap event, including multiple SVs in one DMR and multiple DMRs in one SV
# Unique SVs: Number of SVs that hit a DMR at least once. No duplicate events are counted (if one SV overlaps 10 DMRs, it still counts as one unique hit)
# Unique DMRs: Number of DMRs that were hit by a SV at least once. No duplicate events are counted (if 5 SVs hit one DMR, it still counts as one unique hit)
overlap_summary <- function(dmr_name, dmr, sv_name, sv) {
  hits <- findOverlaps(sv, dmr)
  
  data.frame(
    DMR_file       = dmr_name,
    SV_file        = sv_name,
    total_overlaps = length(hits),
    unique_SVs     = length(unique(queryHits(hits))),
    unique_DMRs    = length(unique(subjectHits(hits)))
  )
}

# Cycles through DMR name x4 and SV file x4
results <- bind_rows(purrr::map2(
  rep(names(dmr_gr), each = length(sv_gr)),
  rep(names(sv_gr), times = length(dmr_gr)),
  ~ overlap_summary(.x, dmr_gr[[.x]], .y, sv_gr[[.y]])
))


print(results)
write.csv(results, "/DMR_SV_overlap_summary.csv", row.names = FALSE)
```
##### SNPs vs. CpG Sites 




##### Gene body methylation vs. TPM values 
```{r}
library(ggplot2)
library(dplyr)
library(ggpubr)

df <- read.delim(
  "/methvtpm.txt",
  header = TRUE, sep = "\t"
)

# Convert % to numeric

meth_cols <- grep("Meth$", colnames(df))
df[meth_cols] <- lapply(df[meth_cols], function(x) as.numeric(gsub("%", "", x)))

# Define Meth/TPM pairs

comparisons <- list(
  "PAH Sensitive Control" = c("PAH_Sensitive_Control_Meth", "PAH_Sensitive_Control_TPM"),
  "PAH Sensitive + PAH"   = c("PAH_Sensitive_.PAH_Meth",   "PAH_Sensitive_.PAH_TPM"),
  "PAH Tolerant Control"  = c("PAH_Tolerant_Control_Meth",  "PAH_Tolerant_Control_TPM"),
  "PAH Tolerant + PAH"    = c("PAH_Tolerant_.PAH_Meth",    "PAH_Tolerant_.PAH_TPM")
)

outdir <- "/output_dir"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Compute global Y-axis limits (log TPM)

all_tpm_cols <- unlist(lapply(comparisons, function(x) x[2]))
global_y <- log2(unlist(df[all_tpm_cols]) + 1)

ymin <- floor(min(global_y, na.rm = TRUE))
ymax <- ceiling(max(global_y, na.rm = TRUE))

for (cond in names(comparisons)) {

  meth_col <- comparisons[[cond]][1]
  tpm_col  <- comparisons[[cond]][2]

  temp <- df %>%
    select(Gene,
           Methylation = all_of(meth_col),
           Expression  = all_of(tpm_col)) %>%
    mutate(logExpression = log2(Expression + 1))

  p <- ggplot(temp, aes(x = Methylation, y = logExpression)) +
    geom_point(color = "blue", alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    stat_cor(aes(label = ..rr.label..),
             size = 5,
             label.x = 5,
             label.y = ymax - 0.5) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_y_continuous(
      limits = c(ymin, ymax),
      breaks = seq(ymin, ymax, by = 1)
    ) +
    theme_classic(base_size = 14) +
    labs(
      x = "Gene Methylation (%)",
      y = "Gene Expression (TPM)",
      title = paste0(cond, "\nGene Methylation vs Expression\n Genes")
    ) +
    theme(plot.title = element_text(hjust = 0.5))

  # SAVE SVG
  outfile_svg <- paste0(outdir, "/", gsub(" ", "_", cond), "_Meth_vs_Exp_logTPM.svg")
  ggsave(outfile_svg, p, width = 6, height = 5, device = "svg")

  print(p)
}
```
