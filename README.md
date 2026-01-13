# smallRNA_seq_pipeline

Pipeline for reanalysis of public small RNA-seq datasets, with an emphasis on tRNA fragment profiling.

## Overview

This pipeline runs end-to-end automatically. The only required user input is a list of SRA run accessions (SRR codes).

Once started, the pipeline builds reference annotations and alignment indexes, downloads data, trims adapters, performs sequential alignments, classifies tRNA fragments, merges counts, and outputs analysis-ready tables.

The pipeline was developed and tested on **human** brain small RNA-seq dataset GSE48552 and should be applicable to other small RNA-seq datasets.

> **Note:**
> - Reference files included in this repository are currently configured for **human RNA**.
> - The pipeline itself is **species-agnostic** and will run correctly with references from other organisms, provided that appropriate reference **FASTA (`.fa`) files** are supplied and the species identifier in the scripts is updated (e.g., changing `human` → `mouse`).
> - Alignment indexes are generated automatically during the `make_reference` step.
> - Enabling the `--qc` flag runs QC on the **output of each pipeline step** (e.g., raw reads after download, trimmed reads after adapter removal, and files generated after each alignment step). When analyzing a new dataset, it is recommended to first run the pipeline on a small set of test samples with `--qc`, then proceed with the full dataset or adjust parameters as needed.
> - Each pipeline step can also be run independently by executing the corresponding script in `scripts/`, provided its required inputs are available. For example:
>   ```bash
>   bash scripts/1_download.sh [--qc]
>   ```

## Alignment strategy

Following adapter trimming and contaminant removal (UniVec), reads are sequentially aligned to:
- ribosomal RNA (rRNA)
- transfer RNA (tRNA), followed by tRNA fragment (tRF) classification
- microRNA, snoRNA, and scaRNA
- piRNA
- protein-coding mRNA (longest transcripts only), long noncoding RNA (lncRNA), and miRNA hairpin references

Reference annotations are derived from RNAcentral (rRNA), GtRNAdb (nuclear tRNA), miRBase (miRNAs), piRBase (piRNAs), and Ensembl (mt-tRNA, snRNA, snoRNA, other ncRNAs, mRNA, and lncRNA).

UniVec filtering is performed using Bowtie (v1) only. All other alignment steps use a two-stage Bowtie-based strategy, with Bowtie (v1) as the primary aligner and Bowtie2 used to rescue reads that fail Bowtie1 alignment due to short overhangs or end mismatches.

Reads are assigned using a best-alignment strategy. When a read maps equally well to multiple features within the same alignment step, **counts are assigned fractionally** to avoid overcounting, and output count values may therefore be non-integer.

Using this multi-step alignment workflow, approximately **>90% of reads** are successfully mapped across samples.

## tRNA fragment classification

tRNA fragments are classified following the tRAX framework (PMID: 40444975) based on read position relative to mature tRNA boundaries:
- reads within 10 nucleotides of the 5′ end are classified as 5′ tRFs
- reads within 10 nucleotides of the 3′ end are classified as 3′ tRFs
- reads mapping to internal regions and not ovlerapping either boundaries are classified as internal tRFs
- reads spanning both boundaries are classified as full-length tRNAs (mostly not detected due to the short read length of <50 nt of small RNA-seq libraries)

## Adapter trimming

The trimming step is configured for TruSeq small RNA libraries. Only reads 15 nt or longer are retained.

If a different library construction was used, adapter sequences and trimming parameters must be modified accordingly.

## Usage

Edit the run list at `scripts/runlist.txt`.

The current run list contains two test runs:
- SRR1103937
- SRR1103939

Run from the repository root: `bash scripts/pipeline.sh [--qc]`.

## Output

Final results are written to `analysis/`. 

Running the pipeline with the default run list will produce results in the existing `analysis/` folder.

Key output files include:
- rna_counts.csv
- mapped_reads_by_source.csv
- percentage_mapped_by_source.csv
 
The primary output `rna_counts.csv` is a raw count table generated directly from the sequential alignment workflow. This table includes all detected RNA species and tRNA-derived fragments prior to any downstream filtering, normalization, or statistical analysis.

Additional filtering (e.g., low-count filtering, biotype selection, sample exclusion) should be applied separately depending on the analysis context.



## Dependencies

Tested on Linux (Ubuntu / WSL).

Required software:
- bash (GNU bash 5.2.21)
- standard Unix utilities (awk, sed, grep, sort, cut, gzip, coreutils)
- wget
- SRA Toolkit (fasterq-dump 3.0.3)
- cutadapt
- Bowtie v1.3.1
- Bowtie2 v2.5.2
- samtools v1.19.2
- R v4.3.3

Required R packages:
- tidyverse
- data.table
