
# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)
# Arguments ---------------------------------------------------------------


args <- commandArgs(TRUE)
cutoff <- ifelse(is.na(args[1]), 0.3, as.numeric(args[1]))

# Load data ---------------------------------------------------------------

cc.se <- readr::read_rds(file = "data/rda/cc.se.rds")
wuhan.tom.fs.fg.norm.rbe.se <- readr::read_rds(file = "data/rda/wuhan.tom.fs.fg.norm.rbe.se.rds.gz")
confounding <- readr::read_rds(file = "data/rda/confounding-svobj.rds")
# Function ----------------------------------------------------------------

fn_filter_samples_by_mapped_reads <- function(.se) {
  .lgl <- .se$`__mapped_reads` > 5e6 & .se$mapping_rate >= cutoff
  .se$Disease
  .se$class <- factor(ifelse(.se$Disease == "Healthy", "B", "M"), levels = c("M", "B"))
  .se[, .lgl]
}

fn_filter_genes <- function(.se) {

  .keep_genes_c <- fn_filter_by_gene_count(.se = .se)

  .se[.keep_genes_c, ]$class %>%
    levels() %>%
    purrr::map(.f = fn_filter_by_hyper, .se[.keep_genes_c, ]) %>%
    purrr::reduce(.f = dplyr::left_join, by = "name") %>%
    dplyr::filter_if(.predicate = is.numeric, .vars_predicate = dplyr::all_vars(. < 3)) %>%
    dplyr::pull(name) ->
    .keep_genes_ineq

  .keep_genes_ineq
  .se[.keep_genes_ineq, ]
}
fn_transform <- function(.se) {
  .vst <- varianceStabilizingTransformation(object = assay(.se), fitType = "local")

  SummarizedExperiment(assays = .vst, colData = .se@colData)
}

fn_rbe <- function(.se, .lsc) {
  # .se <- cc.fs.fg.norm.se
  # .lsc <- confounding
  .newsvobj <- list(
    sv = .lsc$svobj$sv[, .lsc$confounding$confounding],
    pprob.gam = .lsc$svobj$pprob.gam,
    pprob.b = .lsc$svobj$pprob.b,
    n.sv = ncol(.lsc$svobj$sv[, .lsc$confounding$confounding])
  )
  .fsvaObj <- sva::fsva(dbdat = .lsc$dbdat, mod = .lsc$mod, sv = .newsvobj, newdat = assay(.se), method = 'fast')

  SummarizedExperiment(assays = .fsvaObj$new, colData = .se@colData[colnames(.fsvaObj$new), ])
}
# filter genes ------------------------------------------------------------
cc.se@colData %>%
  as.data.frame() %>%
  dplyr::group_by(Disease) %>%
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-n)
cc.fs.se <- fn_filter_samples_by_mapped_reads(cc.se)
cc.fs.se@colData %>%
  as.data.frame() %>%
  dplyr::group_by(Disease) %>%
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-n)


# Filter genes ------------------------------------------------------------


cc.fs.fg.se <- cc.fs.se[rownames(wuhan.tom.fs.fg.norm.rbe.se),]


# Transform ---------------------------------------------------------------

cc.fs.fg.norm.se <- fn_transform(.se = cc.fs.fg.se)
readr::write_rds(x = cc.fs.fg.norm.se, file = "data/rda/cc.fs.fg.norm.se.rds")

# RBE ---------------------------------------------------------------------

cc.fs.fg.norm.rbe.se <- fn_rbe(.se = cc.fs.fg.norm.se, .lsc = confounding)
readr::write_rds(x = cc.fs.fg.norm.rbe.se, file = "data/rda/cc.fs.fg.norm.rbe.se.rds")
# Save image --------------------------------------------------------------

save.image(file = "data/rda/cc-2.rda")

