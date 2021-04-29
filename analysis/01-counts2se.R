
# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)
library(doParallel)

# src ---------------------------------------------------------------------

source(file = "src/doparallel.R")
source(file = "src/utils.R")

# Function ----------------------------------------------------------------

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
# Load data ---------------------------------------------------------------


# Load wuhan data ---------------------------------------------------------
wuhan_count <- fn_load_data(.path = "data/raw/wuhan")
readr::write_rds(x = wuhan_count, file = "data/rda/wuhan_count.rds.gz", compress = "gz")

tom_count <- fn_load_data(.path = "data/raw/tom")
readr::write_rds(x = tom_count, file = "data/rda/tom_count.rds.gz", compress = "gz")


# Load meta ---------------------------------------------------------------
wuhan_meta <- readr::read_rds(file = "data/rda/wuhan_meta.rds.gz")
wuhan_update_stage <- readxl::read_xlsx(path = "data/raw/wuhan/wuhan-stage.xlsx") %>%
  dplyr::slice(match(wuhan_meta$barcode, barcode))
wuhan_meta$figo_stage <- wuhan_update_stage$figo_stage

tom_meta <- readr::read_rds(file = "data/rda/tom_meta.rds.gz")




# Count meta se -----------------------------------------------------------
wuhan.se <- fn_count_meta_2_se(.count = wuhan_count, .meta = wuhan_meta)
readr::write_rds(x = wuhan.se, file = "data/rda/wuhan.se.rds.gz", compress = "gz")
tom.se <- fn_count_meta_2_se(.count = tom_count, .meta = tom_meta)
readr::write_rds(x = tom.se, file = "data/rda/tom.se.rds.gz", compress = "gz")


# Save image --------------------------------------------------------------

save.image(file = "data/rda/01-counts2se.rda")
load(file = "data/rda/01-counts2se.rda")
