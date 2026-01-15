#!/usr/bin/env Rscript

################################################################################
# 14_merge.R
#
# Step 14: Merge RNA count tables and compute mapped reads / percentages
#
# Outputs:
#   analysis/rna_counts.csv
#   analysis/mapped_reads_by_biotype_2.csv
#   analysis/percentage_mapped_by_biotype_2.csv
#
################################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(fs)
  library(tidyr)
})

# ================================
# Paths
# ================================
ROOT_DIR     <- getwd()
COUNTS_DIR   <- file.path(ROOT_DIR, "counts")
ANALYSIS_DIR <- file.path(ROOT_DIR, "analysis")
LOG_FILE     <- file.path(ROOT_DIR, "log/rrna_bowtie1/read_counts.csv")

dir_create(ANALYSIS_DIR)

OUT_COUNTS   <- file.path(ANALYSIS_DIR, "rna_counts.csv")
OUT_READS    <- file.path(ANALYSIS_DIR, "mapped_reads_by_biotype_2.csv")
OUT_PCT      <- file.path(ANALYSIS_DIR, "percentage_mapped_by_biotype_2.csv")

message("📂 Pipeline root: ", ROOT_DIR)

# ================================
# Safety checks
# ================================
if (!dir_exists(COUNTS_DIR)) stop("❌ Missing counts/")
if (!file.exists(LOG_FILE)) stop("❌ Missing read_counts.csv")

# ================================
# Constants
# ================================
EXCLUDE_RE <- "/trna_bowtie[12]/"

meta_cols <- c(
  "genename","geneid","biotype_1","biotype_2",
  "chromosome","start","end","strand","length","source"
)

# ================================
# Discover count tables
# ================================
bowtie1_files <- dir_ls(COUNTS_DIR, recurse = TRUE, regexp = "bowtie1.*\\.csv$")
bowtie2_files <- dir_ls(COUNTS_DIR, recurse = TRUE, regexp = "bowtie2.*\\.csv$")

bowtie1_files <- bowtie1_files[!grepl(EXCLUDE_RE, bowtie1_files)]
bowtie2_files <- bowtie2_files[!grepl(EXCLUDE_RE, bowtie2_files)]

bases <- intersect(
  str_replace(basename(bowtie1_files), "_bowtie1\\.csv$", ""),
  str_replace(basename(bowtie2_files), "_bowtie2\\.csv$", "")
)

# ================================
# Merge helper
# ================================
merge_two <- function(f1, f2) {
  message("  🔗 ", f1, " + ", f2)

  d1 <- read_csv(f1, show_col_types = FALSE)
  d2 <- read_csv(f2, show_col_types = FALSE)

  sample_cols <- union(
    setdiff(names(d1), meta_cols),
    setdiff(names(d2), meta_cols)
  )

  merged <- full_join(d1, d2, by = "genename", suffix = c(".f1",".f2"))

  for (m in setdiff(meta_cols, "genename")) {
    c1 <- paste0(m, ".f1"); c2 <- paste0(m, ".f2")
    if (c1 %in% names(merged) && c2 %in% names(merged)) {
      merged[[m]] <- ifelse(!is.na(merged[[c1]]), merged[[c1]], merged[[c2]])
      merged[[c1]] <- merged[[c2]] <- NULL
    }
  }

  for (s in sample_cols) {
    c1 <- paste0(s, ".f1"); c2 <- paste0(s, ".f2")
    v1 <- if (c1 %in% names(merged)) ifelse(is.na(merged[[c1]]),0,merged[[c1]]) else 0
    v2 <- if (c2 %in% names(merged)) ifelse(is.na(merged[[c2]]),0,merged[[c2]]) else 0
    merged[[s]] <- v1 + v2
    if (c1 %in% names(merged)) merged[[c1]] <- NULL
    if (c2 %in% names(merged)) merged[[c2]] <- NULL
  }

  merged[, c(meta_cols, sample_cols)]
}

# ================================
# Build merged table
# ================================
tables <- list()

for (b in bases) {
  f1 <- bowtie1_files[str_detect(bowtie1_files, paste0("/", b, "_bowtie1\\.csv$"))]
  f2 <- bowtie2_files[str_detect(bowtie2_files, paste0("/", b, "_bowtie2\\.csv$"))]
  if (length(f1)==1 && length(f2)==1)
    tables[[b]] <- merge_two(f1, f2)
}

single_counts <- dir_ls(COUNTS_DIR, recurse = TRUE, regexp = "_counts\\.csv$")
single_counts <- single_counts[!grepl("bowtie[12]", single_counts)]
single_counts <- single_counts[!grepl(EXCLUDE_RE, single_counts)]

for (f in single_counts) {
  message("  ➕ ", f)
  tables[[basename(f)]] <- read_csv(f, show_col_types = FALSE)
}

combined <- bind_rows(tables)

# Mt-tRF normalization
combined <- combined %>%
  mutate(
    biotype_2 = ifelse(
      biotype_1 == "Mt_tRNA" & !startsWith(biotype_2,"Mt-"),
      paste0("Mt-", biotype_2),
      biotype_2
    )
  )

sample_cols <- setdiff(names(combined), meta_cols)

combined <- combined %>%
  filter(rowSums(across(all_of(sample_cols), as.numeric)) > 0)

write_csv(combined, OUT_COUNTS)
message("💾 rna_counts.csv written")

# ================================
# Per-sample mapped reads
# ================================

orig_reads <- read_csv(LOG_FILE, show_col_types = FALSE) %>%
  rename(Original_Reads = Reads_Processed)

mapped_long <- combined %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Mapped_Reads"
  ) %>%
  group_by(Sample, biotype_2) %>%
  summarise(
    Mapped_Reads = sum(Mapped_Reads, na.rm = TRUE),
    .groups = "drop"
  )

reads_wide <- mapped_long %>%
  pivot_wider(
    names_from  = biotype_2,
    values_from = Mapped_Reads,
    names_glue  = "{biotype_2}_reads"
  )

reads_final <- orig_reads %>%
  left_join(reads_wide, by = "Sample")

# explicitly define read columns
read_cols <- grep("_reads$", names(reads_final), value = TRUE)

# replace NA → 0
reads_final[read_cols] <- lapply(reads_final[read_cols], function(x) {
  ifelse(is.na(x), 0, as.numeric(x))
})

# compute unmapped reads
reads_final$unmapped_reads <-
  reads_final$Original_Reads -
  rowSums(reads_final[read_cols])

write_csv(reads_final, OUT_READS)
message("💾 mapped_reads_by_biotype_2.csv written")

# ================================
# Per-sample percentages (WIDE)
# ================================

pct_final <- reads_final

# compute per-biotype percentages
for (col in read_cols) {
  pct_col <- sub("_reads$", "_pct", col)
  pct_final[[pct_col]] <-
    round(100 * pct_final[[col]] / pct_final$Original_Reads, 3)
}

# explicitly define pct columns
pct_cols <- grep("_pct$", names(pct_final), value = TRUE)

# compute unmapped percent (EXPLICIT)
pct_final$unmapped_pct <-
  round(100 - rowSums(pct_final[pct_cols]), 3)

pct_final <- pct_final %>%
  select(
    Sample,
    Original_Reads,
    all_of(pct_cols),
    unmapped_pct
  )

write_csv(pct_final, OUT_PCT)
message("💾 percentage_mapped_by_biotype_2.csv written")
