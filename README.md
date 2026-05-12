# smallRNA_seq_pipeline

Pipeline for reanalysis of public **single-end (unpaired) small RNA-seq datasets**, with an emphasis on tRNA fragment profiling.

## Overview

This pipeline runs end-to-end automatically. The only required user input is a list of SRA run accessions (SRR codes) in `input/runlist.txt`. The `input/metadata.csv` file is only required for the sample analysis step.

Once started, the pipeline builds reference annotations and alignment indexes, downloads data, trims adapters, performs sequential alignments, classifies tRNA fragments, merges counts, and outputs analysis-ready tables.

The pipeline was developed and tested on the following publicly available datasets:

| GEO | Species | Tissue | Library Kit | Platform | Read length | PMID |
|-----|---------|--------|-------------|----------|-------------|------|
| GSE48552 | Human | Prefrontal cortex | TruSeq Small RNA (Illumina) | HiSeq 2000 | 50 nt | 24014289 |
| GSE111279 | Mouse (C57BL/6) | Brain hemisphere (ages 2–30 mo) | TruSeq Small RNA (Illumina) | HiSeq 2500 | 50 nt | — |
| GSE217458 | Mouse (C57BL/6JN) | 16 organs (ages 1–27 mo, 10 timepoints) | MGIEasy Small RNA Library Prep (MGI SP-960) | BGISEQ-500RS | 50 nt | 37106037 |
| GSE282205 | Mouse (C57BL/6JN) | 15 brain regions (ages 3–28 mo) | MGIEasy Small RNA Library Prep (MGI SP-960) | BGISEQ-500 | 50 nt | 40382330 |
| E-MTAB-12731 | Human | Frontal lobe, post-mortem (FTD + controls) | NEXTflex Small RNA-seq v3 (Bioo Scientific) | NextSeq 550 | 50 nt | 38040703 |
| GSE131695 | Human | CSF (24 h post-TBI vs. control) | NEBNext Small RNA Library Prep (NEB) | NextSeq 550 | 75 nt | 32094409 |

The pipeline is **species-agnostic** and should be applicable to other single-end small RNA-seq datasets with appropriate modifications, as outlined below.

> **Notes:**
>
> - Reference files included in this repository support both **human** and **mouse**. The species is selected at runtime via `--species mouse` or `--species human`.
> - Alignment indexes are generated automatically from the provided `.fa` files during `00_make_reference.sh`.
> - The adapter trimming configuration in `02_trim.sh` is currently set for **TruSeq Small RNA libraries** (GSE48552, GSE111279) and retains reads **≥15 nt**. Other library kits require updating the adapter sequences in `02_trim.sh` before running. Read length thresholds may also need adjustment.
> - Enabling `--qc` runs **FastQC and MultiQC** on the output of each pipeline step. When analyzing a new dataset, it is recommended to first run on a small test set with `--qc`, then proceed with the full dataset.
> - Each pipeline step can be run independently, provided its required inputs are available. For example:
>
>   ```bash
>   bash alignment_scripts/04_align_rrna.sh --species mouse [--qc]
>   ```

## Alignment strategy

Following adapter trimming and contaminant removal (UniVec), reads are sequentially aligned to:

1. Ribosomal RNA (rRNA)
2. Transfer RNA (tRNA), followed by tRNA fragment (tRF) classification
3. miRNA, snoRNA, and other small RNAs
4. piRNA
5. Protein-coding mRNA (longest transcripts only), lncRNA, and miRNA hairpin references

Reference annotations are derived from RNAcentral (rRNA), GtRNAdb (nuclear tRNA), miRBase (miRNAs), piRBase (piRNAs), and Ensembl (mt-tRNA, snRNA, snoRNA, other ncRNAs, mRNA, and lncRNA).

UniVec filtering uses Bowtie (v1) only. All other alignment steps use a two-stage strategy: Bowtie (v1) as the primary aligner, followed by Bowtie2 to rescue reads that fail Bowtie1 alignment due to short overhangs or end mismatches. The tRNA rescue step uses relaxed mismatch/gap penalties to recover reads carrying tRNA modifications (m1A, m3C, pseudouridine, inosine) that cause RT misincorporations.

When a read maps equally well to multiple features within the same step, **counts are assigned fractionally** (NH tag). Output count values may therefore be non-integer.

Using this workflow, an average of **90–95%** of reads were successfully mapped across samples; mapping rates were lower in CSF samples (GSE131695; 50-75%)

## tRNA fragment classification

tRNA fragments are classified following the tRAX framework (PMID: 40444975, 39952700) based on read position relative to mature tRNA boundaries, including the 3′ CCA tail:

| Class | Rule |
|-------|------|
| 5′-tRF | read start within 10 nt of the tRNA 5′ end |
| 3′-tRF | read end within 10 nt of the tRNA 3′ end |
| internal-tRF | read does not overlap either boundary |
| FL-tRNA | read spans both boundaries (rare at <50 nt read length) |

## Usage

> **Important:** Reference FASTA files are stored using [Git LFS](https://git-lfs.com/). After cloning, run `git lfs pull` to download them before running the pipeline.

**1. Edit the run list:**

```
input/runlist.txt   # one SRR accession per line
```

**2. Run the full pipeline:**

```bash
bash alignment_scripts/pipeline.sh --species mouse [--qc]
bash alignment_scripts/pipeline.sh --species human [--qc]
```

The current `runlist.txt` contains SRR codes for the GSE111279 mouse dataset.

## Output

Final results are written to `analysis/data/`.

| File | Description |
|------|-------------|
| `rna_counts.csv` | Raw count table (all RNA species and tRFs) |
| `alignment_stats.csv` | Per-sample mapped reads by RNA biotype |
| `alignment_pct.csv` | Per-sample mapping percentages by RNA biotype |
| `pipeline_runtime.csv` | Per-step runtime log |

`rna_counts.csv` is generated directly from the alignment workflow, prior to any filtering, normalization, or statistical analysis. Downstream filtering (low-count removal, biotype selection, sample exclusion) should be applied based on the analysis context.

## Dependencies

Tested on **Linux (Ubuntu / WSL2)**. Install all dependencies at once using [pixi](https://pixi.sh):

```bash
pixi install
pixi shell   # activates the environment
```

The `pixi.toml` in the repository root pins all required tools:

| Tool | Purpose |
|------|---------|
| bowtie ≥1.3.1 | Primary aligner (rRNA, tRNA, small RNA, piRNA, longRNA) |
| bowtie2 ≥2.5.1 | Rescue aligner |
| samtools ≥1.17 | BAM processing |
| subread ≥2.0.6 | featureCounts quantification |
| cutadapt ≥4.4 | Adapter trimming |
| fastqc ≥0.12.1 | Per-sample QC |
| multiqc ≥1.19 | Aggregated QC reports |
| wget ≥1.20 | SRA FASTQ download |
| R ≥4.3 + Bioconductor/CRAN packages | tRF classification, count merging, analysis |

Standard Unix utilities (awk, gzip, sort, etc.) are assumed to be available in the environment.

## License

This project is licensed under the MIT License.
