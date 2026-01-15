# smallRNA_seq_pipeline

Pipeline for reanalysis of public small RNA-seq datasets, with an emphasis on tRNA fragment profiling.

## Overview

This pipeline runs end-to-end automatically. The only required user input is a list of SRA run accessions (SRR codes) in runlist.txt in input folder. The metadata.csv file in the input folder is only required for the sample analysis step.

Once started, the pipeline builds reference annotations and alignment indexes, downloads data, trims adapters, performs sequential alignments, classifies tRNA fragments, merges counts, and outputs analysis-ready tables.

The pipeline was developed and tested on **human** brain small RNA-seq dataset GSE48552 and should be applicable to other small RNA-seq datasets.

> **Note:**
> - Reference files included in this repository are currently configured for **human RNA**.
> - The pipeline itself is **species-agnostic** and will run correctly with references from other organisms, provided that appropriate reference **FASTA (`.fa`) files** are supplied and the species identifier in the scripts is updated (e.g., changing `human` → `mouse`).
> - Alignment indexes are generated automatically from the provided `.fa` files during the `00_make_reference.sh` step.
> - The adapter trimming step (`02_trim.sh`) is currently configured for **TruSeq small RNA libraries** and retains reads **≥15 nt**. This step may need to be adjusted depending on how the library was constructed (e.g., adapter sequences, read length distribution) and the overall quality of the reads.
> - Enabling the `--qc` flag runs **FastQC and MultiQC** on the output of each pipeline step (e.g., raw reads after download, trimmed reads after adapter removal, and files generated after each alignment step). When analyzing a new dataset, it is recommended to first run the pipeline on a small set of test samples with `--qc`, then proceed with the full dataset or adjust parameters as needed.
> - Each pipeline step can also be run independently by executing the corresponding script in `scripts/`, provided its required inputs are available. For example:
>   ```bash
>   bash scripts/1_download.sh [--qc]
>   ```

## Alignment strategy

Following adapter trimming and contaminant removal (UniVec), reads are sequentially aligned to:
- ribosomal RNA (rRNA)
- transfer RNA (tRNA), followed by tRNA fragment (tRF) classification
- miRNA, snoRNA, and other small RNAs
- piRNA
- protein-coding mRNA (longest transcripts only), long noncoding RNA (lncRNA), and miRNA hairpin references

Reference annotations are derived from RNAcentral (rRNA), GtRNAdb (nuclear tRNA), miRBase (miRNAs), piRBase (piRNAs), and Ensembl (mt-tRNA, snRNA, snoRNA, other ncRNAs, mRNA, and lncRNA).

UniVec filtering is performed using Bowtie (v1) only. All other alignment steps use a two-stage Bowtie-based strategy, with Bowtie (v1) as the primary aligner and Bowtie2 used to rescue reads that fail Bowtie1 alignment due to short overhangs or end mismatches.

Reads are assigned using a best-alignment strategy. When a read maps equally well to multiple features within the same alignment step, **counts are assigned fractionally** to avoid overcounting, and output count values may therefore be non-integer.

Using this multi-step alignment workflow, approximately **>90% of reads** are successfully mapped across samples.

## tRNA fragment classification

tRNA fragments are classified following the tRAX framework (PMID: 40444975) based on read position relative to mature tRNA boundaries, including the 3′ CCA tail.
- reads within 10 nucleotides of the 5′ end are classified as 5′ tRFs
- reads within 10 nucleotides of the 3′ end are classified as 3′ tRFs
- reads mapping to internal regions and not overlapping either boundaries are classified as internal tRFs
- reads spanning both boundaries are classified as full-length tRNAs (mostly not detected due to the short read length of <50 nt of small RNA-seq libraries)

## Usage

Edit the run list at `input/runlist.txt`.

The current run list contains 12 SRR codes for samples from GSE48552.

Run from the repository root: `bash scripts/pipeline.sh [--qc]`.

To run in docker: 

```bash
docker run --rm -it
-v "$PWD/input:/opt/smallRNA_seq_pipeline/input"
-v "$PWD/analysis:/opt/smallRNA_seq_pipeline/analysis"
smallrna-pipeline
bash scripts/pipeline.sh
```

## Output

Final results are written to `analysis/`. 

Running the pipeline with the default run list will produce results in the existing `analysis/` folder.

Key output files include:
- rna_counts.csv
- mapped_reads_by_biotype_2.csv
- percentage_mapped_by_biotype2.csv
 
The primary output `rna_counts.csv` is a raw count table generated directly from the sequential alignment workflow. This table includes all detected RNA species and tRNA-derived fragments prior to any downstream filtering, normalization, or statistical analysis.

Additional filtering (e.g., low-count filtering, biotype selection, sample exclusion) should be applied depending on the analysis context.

#  Sample Analysis Step

This script (`15_analysis.R`) performs a global differential expression analysis using DESeq2 after collapsing tRF isodecoders. Differential expression is computed across all retained RNA features using sample metadata provided in `input/metadata.csv`.

**Figure 1h** displays only nuclear and mitochondrial tRFs extracted from the global DESeq2 results

**Figure 1i** shows a principal component analysis (PCA) of tRNA and mitochondrial tRNA expression.

This script is provided as a reproducible example of the analysis workflow used to generate key figures in the manuscript.



## Dependencies

Tested on Linux (Ubuntu / WSL).

Required software:
- bash (GNU bash 5.2.21)
- standard Unix utilities (awk, sed, grep, sort, cut, gzip, coreutils)
- wget
- SRA Toolkit (fasterq-dump 3.0.3)
- cutadapt
- FastQC
- MultiQC
- Bowtie v1.3.1
- Bowtie2 v2.5.2
- samtools v1.19.2
- R v4.3.3

Required R packages:
- tidyverse
- data.table

