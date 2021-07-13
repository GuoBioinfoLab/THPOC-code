
# Library -----------------------------------------------------------------

library(magrittr)
library(doParallel)
library(DESeq2)

# src ---------------------------------------------------------------------

source(file = "src/doparallel.R")
source(file = "src/utils.R")

# Load data ---------------------------------------------------------------


# Functions ---------------------------------------------------------------

fn_load_data <- function(.path) {
  # .path = "data/raw/wuhan"
  .filenames <- list.files(path = .path, pattern = "htseq_count|htSeqCounts")

  fn_parallel_start(n_cores = 50)

  foreach(
    x = .filenames,
    .packages = c("magrittr"),
    .export = c("fn_load_htseq_count")
  ) %dopar% {
    fn_load_htseq_count(.txt = x, .path = .path)
  }  ->
    .d

  fn_parallel_stop()

  .d %>% purrr::reduce(.f = dplyr::inner_join, by = "Ensembl Gene ID")
}

fn_count_meta_2_se <- function(.count, .meta) {
  # .count <- tom_count
  # .meta <- tom_meta
  .samples <- intersect(colnames(.count), rownames(.meta))

  .clean <- fn_filter_annotation(.count)

  .matrix <- .clean$clean_count %>%
    dplyr::select(-`Ensembl Gene ID`) %>%
    as.matrix()
  .matrix <- .matrix[, .samples]
  rownames(.matrix) <- .clean$clean_count$`Ensembl Gene ID`

  .meta[colnames(.matrix), ] %>%
    as.data.frame() %>%
    dplyr::left_join(.clean$mapping_rate, by = "barcode") ->
    .coldata
  rownames(.coldata) <- .coldata$barcode

  SummarizedExperiment(assays = .matrix, colData = .coldata)
}
# load data ---------------------------------------------------------------

cc_count <- fn_load_data(.path = "data/data_cancer_cell_female")


# Load meta ---------------------------------------------------------------
meta <- readxl::read_xlsx(path = "data/data_cancer_cell_female/05-cancer-cell-female.xlsx", sheet = 2) %>%
  dplyr::mutate(barcode = srr) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "srr")


# count meta se -----------------------------------------------------------

cc.se <- fn_count_meta_2_se(.count = cc_count, .meta = meta)
readr::write_rds(x = cc.se, file = "data/rda/cc.se.rds")
# Save image --------------------------------------------------------------

save.image(file = "data/rda/cc.rda")
load(file = "data/rda/cc.rda")
