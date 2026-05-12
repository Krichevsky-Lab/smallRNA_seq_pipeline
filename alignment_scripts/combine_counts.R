#!/usr/bin/env Rscript
# combine_counts.R — merge per-sample featureCounts files and join annotation
#
# Usage:
#   Rscript alignment_scripts/combine_counts.R \
#     <counts_dir> <file_pattern> <suffix_to_strip> <annot_tsv> <output_csv>
#
# Arguments:
#   counts_dir     directory containing *_counts.txt files
#   file_pattern   regex to match count files (e.g. "_rrna_bowtie1_counts\\.txt$")
#   suffix_strip   literal suffix removed from filename to get sample name
#   annot_tsv      annotation TSV produced by 00_make_reference.sh
#   output_csv     path to write the merged CSV

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: combine_counts.R <counts_dir> <pattern> <suffix> <annot> <output>")
}

counts_dir  <- args[1]
pattern     <- args[2]
suffix      <- args[3]
annot_file  <- args[4]
output_csv  <- args[5]

files <- list.files(counts_dir, pattern = pattern, full.names = TRUE)
if (length(files) == 0) {
  message("[INFO]  No files matching '", pattern, "' in ", counts_dir, " — skipping.")
  quit(status = 0)
}

ann <- read_tsv(annot_file, show_col_types = FALSE)

df_list <- lapply(files, function(f) {
  sample <- sub(suffix, "", basename(f))
  sample <- sub("_trimmed$", "", sample)
  dat <- read_tsv(f, comment = "#", show_col_types = FALSE)
  count_col <- names(dat)[ncol(dat)]
  dat |>
    rename(genename = Geneid) |>
    select(genename, Length, !!sym(count_col)) |>
    rename(fc_length = Length, !!sample := !!sym(count_col))
})

fc <- Reduce(function(x, y) full_join(x, y, by = c("genename", "fc_length")), df_list)

meta_cols <- c("genename", "geneid", "biotype_1", "biotype_2",
               "chromosome", "start", "end", "strand", "length", "source")

merged <- fc |>
  left_join(ann, by = "genename") |>
  mutate(length = coalesce(length, fc_length)) |>
  select(all_of(meta_cols), everything(), -fc_length)

write_csv(merged, output_csv)
message("[DONE]  ", output_csv, " (", nrow(merged), " genes, ",
        ncol(merged) - length(meta_cols), " samples)")
