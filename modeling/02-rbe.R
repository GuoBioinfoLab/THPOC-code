
# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)
library(sva)
library(doParallel)
library(limma)

# src ---------------------------------------------------------------------

source(file = "src/doparallel.R")
source(file = "src/utils.R")


# Load data ---------------------------------------------------------------

wuhan.se <- readr::read_rds(file = "data/rda/wuhan.se.rds.gz")
tom.se <- readr::read_rds(file = "data/rda/tom.se.rds.gz")
wuhan.tom.fs.fg.norm.se <- readr::read_rds(file = "data/rda/wuhan.tom.fs.fg.norm.se.rds.gz")

# Function ----------------------------------------------------------------
fn_unwanted_vairables <- function(.w, .t, .wt) {

  .w@colData %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(cohort = plyr::revalue(oc, c("OC521" = "TC", "OC44" = "DC", "OC79" = "VC1", "OC172" = "VC2"))) %>%
    dplyr::select(barcode, cohort, age, lib.size, library) ->
    .wcoldata

  .t@colData %>%
    as.data.frame() %>%
    dplyr::mutate(cohort = "Tom") %>%
    dplyr::mutate(age = ageAtBiopsy) %>%
    dplyr::mutate(library = "smart2") %>%
    dplyr::mutate(lib.size = X__mapped_reads) %>%
    dplyr::select(barcode, cohort, age, lib.size, library) ->
    .tcoldata

  .wt@colData %>%
    as.data.frame() %>%
    dplyr::left_join(dplyr::bind_rows(.wcoldata, .tcoldata), by = "barcode") ->
    .wtcoldata

  rownames(.wtcoldata) <- .wtcoldata$barcode

  .se <- SummarizedExperiment(assays = assay(.wt), colData = .wtcoldata)

  fn_parallel_start(n_cores = 50)
  .rbe.se <- fn_remove_unwanted_variables(.se = .se)
  fn_parallel_stop()

  .rbe.se

}

# Unwanted variables ------------------------------------------------------

wuhan.tom.fs.fg.norm.rbe.se <- fn_unwanted_vairables(.w = wuhan.se, .t = tom.se, .wt = wuhan.tom.fs.fg.norm.se)

readr::write_rds(wuhan.tom.fs.fg.norm.rbe.se, file = "data/rda/wuhan.tom.fs.fg.norm.rbe.se.rds.gz")


# Save image --------------------------------------------------------------

save.image(file = "data/rda/02-rbe.rda")

