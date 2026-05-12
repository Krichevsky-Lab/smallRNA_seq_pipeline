#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# classify_tRF.R
#
# Classify tRNA-derived fragments (tRFs) from BAM alignments and quantify
# fragment abundance using fractional assignment (NH tag).
#
# Inputs:
#   - BAM files aligned to tRNA reference
#   - SAF annotation defining tRNA loci
#   - tRNA annotation TSV
#
# Outputs:
#   - Per-sample tRF count tables
#   - Combined annotated tRF count table (CSV)
#
# Notes:
#   - Fragment classes: FL-tRNA, 5'-tRF, 3'-tRF, internal-tRF
#   - mt-tRNA gene IDs are normalized to Ensembl IDs before annotation join
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(tools)
  library(tidyr)
})

# ----------------------------------------------------------------------
# Parse arguments
# ----------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: classify_tRF.R <bam_dir> <saf_file> <ann_file> <out_csv> <per_sample_dir>")
}

bam_dir    <- args[1]
saf_file   <- args[2]
ann_file   <- args[3]
out_csv    <- args[4]
per_sample <- args[5]
yield_size <- 1e6

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(per_sample, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------
# Load SAF annotation
# ----------------------------------------------------------------------
saf <- read_tsv(saf_file, show_col_types = FALSE)

trnas <- GRanges(
  seqnames = saf$Chr,
  ranges   = IRanges(saf$Start, saf$End),
  strand   = rep("*", nrow(saf))
)
mcols(trnas)$gene_id <- saf$GeneID

# ----------------------------------------------------------------------
# Classifier (10 nt rule)
# ----------------------------------------------------------------------
classify_tRF <- function(read_start, read_end, trna_start, trna_end) {
  len <- trna_end - trna_start + 1L
  rel_start <- read_start - trna_start
  rel_end   <- read_end   - trna_start

  near_5 <- rel_start <= 10L
  near_3 <- rel_end   >= (len - 10L)

  if (near_5 && near_3) "FL-tRNA"
  else if (near_5) "5'-tRF"
  else if (near_3) "3'-tRF"
  else "internal-tRF"
}

# ----------------------------------------------------------------------
# Identify BAM files
# ----------------------------------------------------------------------
bam_files <- list.files(bam_dir, pattern = "_trna.*\\.bam$", full.names = TRUE)
if (length(bam_files) == 0)
  stop("No BAM files found in: ", bam_dir)

sample_names <- sub("_trna.*$", "", file_path_sans_ext(basename(bam_files)))
names(bam_files) <- sample_names

# ----------------------------------------------------------------------
# Process one BAM
# ----------------------------------------------------------------------
process_one_bam <- function(bam_path, sample_id) {

  out_file <- file.path(per_sample, paste0(sample_id, "_tRF.csv"))
  if (file.exists(out_file)) {
    message("[SKIP]  Skipping sample ", sample_id)
    return(read_csv(out_file, show_col_types = FALSE))
  }

  bf <- BamFile(bam_path, yieldSize = yield_size)
  param <- ScanBamParam(
    what = c("qname", "rname", "pos", "qwidth"),
    tag  = c("NH"),
    flag = scanBamFlag(isUnmappedQuery = FALSE)
  )

  open(bf); on.exit(close(bf), add = TRUE)
  message("[INFO]  Processing sample: ", sample_id)

  chunks <- list()

  repeat {
    x <- scanBam(bf, param = param)[[1]]
    if (length(x$qname) == 0) break

    keep <- !is.na(x$pos) & !is.na(x$qwidth) & !is.na(x$rname)
    if (!any(keep)) next

    pos     <- x$pos[keep]
    qwidth  <- x$qwidth[keep]
    nh_tag  <- x$tag$NH[keep]
    nh_tag[is.na(nh_tag) | nh_tag <= 0] <- 1L

    read_start <- pos
    read_end   <- pos + qwidth - 1L

    gr_reads <- GRanges(
      seqnames = Rle(x$rname[keep]),
      ranges   = IRanges(read_start, read_end),
      strand   = Rle("*")
    )
    mcols(gr_reads)$nh <- nh_tag

    hits <- findOverlaps(gr_reads, trnas, ignore.strand = TRUE)
    if (length(hits) == 0) next

    qh <- queryHits(hits)
    sh <- subjectHits(hits)

    df <- tibble(
      trna_id    = mcols(trnas)$gene_id[sh],
      nh         = mcols(gr_reads)$nh[qh],
      read_start = start(gr_reads)[qh],
      read_end   = end(gr_reads)[qh],
      trna_start = start(trnas)[sh],
      trna_end   = end(trnas)[sh]
    ) %>%
      mutate(
        fragment = mapply(classify_tRF, read_start, read_end, trna_start, trna_end),
        count    = 1 / nh
      ) %>%
      group_by(fragment, trna_id) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      mutate(
        genename = paste0(fragment, "-", trna_id),
        geneid   = trna_id,
        biotype_2 = fragment,
        !!sample_id := count
      ) %>%
      select(genename, geneid, biotype_2, !!sample_id)

    chunks[[length(chunks) + 1]] <- df
  }

  df_final <- if (length(chunks) == 0) {
    tibble(genename=character(), geneid=character(), biotype_2=character(),
           !!sample_id := numeric())
  } else {
    bind_rows(chunks) %>%
      group_by(genename, geneid, biotype_2) %>%
      summarise(!!sample_id := sum(.data[[sample_id]]), .groups="drop")
  }

  write_csv(df_final, out_file)
  message("[DONE]  Saved per-sample tRFs -> ", out_file)
  df_final
}

# ----------------------------------------------------------------------
# Process all samples
# ----------------------------------------------------------------------
trf_list <- map2(bam_files, names(bam_files), process_one_bam)

final_df <- reduce(trf_list, full_join,
                   by = c("genename","geneid","biotype_2")) %>%
  mutate(across(-c(genename, geneid, biotype_2),
                ~ replace(., is.na(.), 0)))

# ----------------------------------------------------------------------
# Build unified annotation (SAF + annotation TSV)
# ----------------------------------------------------------------------
ann_file_df <- read_tsv(ann_file, show_col_types = FALSE)

# SAF-based annotation (always valid)
ann_saf <- saf %>%
  transmute(
    geneid      = GeneID,
    biotype_1   = "tRNA",
    chromosome  = Chr,
    start       = Start,
    end         = End,
    strand      = Strand,
    length      = End - Start + 1,
    source      = "saf"
  )

# Annotation TSV (preferred when available)
ann_file_df <- ann_file_df %>%
  transmute(
    geneid      = genename,   # IMPORTANT: match mt-tRNA-Ala, not ENSG ID
    biotype_1,
    chromosome,
    start,
    end,
    strand,
    length,
    source
  )

# Merge SAF + annotation (annotation wins)
unified_ann <- ann_saf %>%
  full_join(ann_file_df, by = "geneid", suffix = c("_saf","_ann")) %>%
  mutate(
    biotype_1  = coalesce(biotype_1_ann, biotype_1_saf, "tRNA"),
    chromosome = coalesce(chromosome_ann, chromosome_saf),
    start      = coalesce(start_ann, start_saf),
    end        = coalesce(end_ann, end_saf),
    strand     = coalesce(strand_ann, strand_saf, "*"),
    length     = coalesce(length_ann, length_saf),
    source     = coalesce(source_ann, source_saf)
  ) %>%
  select(geneid, biotype_1, chromosome, start, end, strand, length, source)

# Join annotation to tRF table
final_df <- final_df %>%
  left_join(unified_ann, by = "geneid")


required <- c("genename","geneid","biotype_1","biotype_2",
              "chromosome","start","end","strand","length","source")
sample_cols <- setdiff(colnames(final_df), required)

out_df <- final_df %>%
  select(genename, geneid,
         biotype_1, biotype_2,
         chromosome, start, end, strand, length, source,
         all_of(sample_cols))

names(out_df) <- sub("_trimmed$", "", names(out_df))

write_csv(out_df, out_csv)
message("[DONE]  Saved final combined tRF table -> ", out_csv)
message("[INFO]  Samples included: ", length(sample_cols))
message("[INFO]  Total tRF entries: ", nrow(out_df))
