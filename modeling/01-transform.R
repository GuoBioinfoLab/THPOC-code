
# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)


# src ---------------------------------------------------------------------
source(file = "src/utils.R")

# Function ----------------------------------------------------------------

fn_filter_samples_by_mapped_reads <- function(.se) {
  .lgl <- .se$`__mapped_reads` > 5e6
  .se[, .lgl]
}

fn_filter_genes <- function(.w, .t) {
  # .w <- wuhan.fs.se
  # .t <- tom.fs.se
  .se1 <- .w
  .se2 <- .t
  .se1@colData <- .se1@colData[, c("barcode", "class")]
  .se2@colData <- .se2@colData[, c("barcode", "class")]

  .se <- cbind(.se1, .se2)

  .keep_genes_c <- fn_filter_by_gene_count(.se = .se)

  .se[.keep_genes_c, ]$class %>%
    levels() %>%
    purrr::map(.f = fn_filter_by_hyper, .xse) %>%
    purrr::reduce(.f = dplyr::left_join, by = "name") %>%
    dplyr::filter_if(.predicate = is.numeric, .vars_predicate = dplyr::all_vars(. < thres)) %>%
    dplyr::pull(name) ->
    .keep_genes_ineq

  .keep_genes_ineq
  .se[.keep_genes_ineq, ]
}

fn_transform <- function(.se) {
  .vst <- varianceStabilizingTransformation(object = assay(.se), fitType = "local")

  SummarizedExperiment(assays = .vst, colData = .se@colData)
}

# Load data ---------------------------------------------------------------

wuhan.se <- readr::read_rds(file = "data/rda/wuhan.se.rds.gz")
tom.se <- readr::read_rds(file = "data/rda/tom.se.rds.gz")


# Filter samples by mapped reads ------------------------------------------


wuhan.fs.se <- fn_filter_samples_by_mapped_reads(wuhan.se)
tom.fs.se <- fn_filter_samples_by_mapped_reads(tom.se)

# Filter genes ------------------------------------------------------------

wuhan.tom.fs.fg.se <- fn_filter_genes(.w = wuhan.fs.se, .t = tom.fs.se)
readr::write_rds(x = wuhan.tom.fs.fg.se, file = "data/rda/wuhan.tom.fs.fg.se.rds.gz", compress = "gz")

# Transform ---------------------------------------------------------------
wuhan.tom.fs.fg.norm.se <- fn_transform(.se = wuhan.tom.fs.fg.se)
readr::write_rds(x = wuhan.tom.fs.fg.norm.se, file = "data/rda/wuhan.tom.fs.fg.norm.se.rds.gz", compress = "gz")


# Save image --------------------------------------------------------------

save.image(file = "data/rda/01-transform.rda")
