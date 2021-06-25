# Meta info -----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-03-1 10:00:58
# @DESCRIPTION: 07-performance-epi.R
# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)
library(mlr)
library(doParallel)
library(ggplot2)

# src ---------------------------------------------------------------------

source(file = "src/doparallel.R")
source(file = "src/utils.R")
source(file = "src/performance.R")
source(file = "src/plots.R")

# Load data ---------------------------------------------------------------
wuhan.se <- readr::read_rds(file = "data/rda/wuhan.se.rds.gz")
tom.se <- readr::read_rds(file = "data/rda/tom.se.rds.gz")
wuhan.tom.fs.fg.norm.rbe.se <- readr::read_rds(file = "data/rda/wuhan.tom.fs.fg.norm.rbe.se.rds.gz")
panel <- readr::read_rds(file = "data/rda/panel.rds.gz")

wuhan.tom.panel.task <- readr::read_rds(file = "data/rda/wuhan.tom.panel.task.rds.gz")
wuhan.tom.panel.ca125.task <- readr::read_rds(file = "data/rda/wuhan.tom.panel.ca125.task.rds.gz")
panel.model <- readr::read_rds(file = "data/rda/panel.model.rds.gz")
ca125.model <- readr::read_rds(file = "data/rda/ca125.model.rds.gz")
panel_ca125.model <- readr::read_rds(file = "data/rda/panel_ca125.model.rds.gz")

# Function ----------------------------------------------------------------


fn_get_epi_task <- function(.list, .w, .t) {
  .w@colData %>%
    as.data.frame() %>%
    dplyr::filter(oc %in% c("OC44", "OC79", "OC172")) %>%
    dplyr::filter(!is.na(epi.non.epi)) %>%
    dplyr::mutate(epi.non.epi = ifelse(type == "normal", "Normal", epi.non.epi)) %>%
    dplyr::select(barcode, epi.non.epi) ->
    .wd

  .t@colData %>%
    as.data.frame() %>%
    dplyr::mutate(epi.non.epi = ifelse(class == "M", "epithelial", "Health")) %>%
    dplyr::mutate(epi.non.epi = ifelse(stageFourGroups == "benign", "Benign", epi.non.epi)) %>%
    dplyr::select(barcode, epi.non.epi) ->
    .td

  dplyr::bind_rows(.wd, .td) %>%
    dplyr::filter(epi.non.epi == "Benign") %>%
    dplyr::pull(barcode) ->
    .benign
  dplyr::bind_rows(.wd, .td) %>%
    dplyr::filter(epi.non.epi == "epithelial") %>%
    dplyr::pull(barcode) ->
    .epi
  dplyr::bind_rows(.wd, .td) %>%
    dplyr::filter(epi.non.epi == "non-epithelial") %>%
    dplyr::pull(barcode) ->
    .nepi

  .benign_epi <- c(.benign, .epi)
  .benign_nepi <- c(.benign, .nepi)

  .list.panel.samples <- .list$panel$samples %>% purrr::reduce(.f = c)
  .list$panel$samples <- list(
    "Epi" = .list.panel.samples[.benign_epi] %>% na.omit() %>% c(),
    "Non-epi" = .list.panel.samples[.benign_nepi] %>% na.omit() %>% c()
  )

  .list.ca125.samples <- .list$ca125$samples %>% purrr::reduce(c)
  .list$ca125$samples <- list(
    "Epi" = .list.ca125.samples[.benign_epi] %>% na.omit() %>% c(),
    "Non-epi" = .list.ca125.samples[.benign_nepi] %>% na.omit() %>% c()
  )
  .list$panel_ca125$samples <- list(
    "Epi" = .list.ca125.samples[.benign_epi] %>% na.omit() %>% c(),
    "Non-epi" = .list.ca125.samples[.benign_nepi] %>% na.omit() %>% c()
  )

  .list
}

fn_performance_epi <- function(.x, .y) {

  .perf <- fn_performance(.x, .y)
  .metrics <- fn_get_metrics(.perf)
  .plot <- fn_get_auc_plot(.perf, .metrics)

  list(
    performance = .perf,
    metrics = .metrics,
    plot = .plot
  )

}
# Prepare tasks -----------------------------------------------------------
epi.task <- fn_get_epi_task(
  .list = c("panel" = list(wuhan.tom.panel.task), wuhan.tom.panel.ca125.task),
  .w = wuhan.se,
  .t = tom.se
)
readr::write_rds(x = epi.task, file = "data/rda/epi.task.rds.gz", compress = "gz")

# Performance -------------------------------------------------------------

purrr::map2(
  .x = list(panel.model, panel_ca125.model, ca125.model),
  .y = epi.task,
  .f = fn_performance_epi
) ->
  epi.performance

names(epi.performance) <- c("panel", "panel_ca125", "ca125")
readr::write_rds(x = epi.performance, file = "data/rda/pei.performance.rds.gz", compress = "gz")

# save metrics and plot
purrr::walk2(
  .x = list("panel", "panel_ca125", "ca125"),
  .y = epi.performance,
  .f = function(.x, .y) {
    readr::write_tsv(x = .y$metrics, file = glue::glue("data/output/epi.{.x}.metrics.tsv"))
    writexl::write_xlsx(x = .y$metrics, path = glue::glue("data/output/epi.{.x}.metrics.xlsx"))
    fn_save_auc(.filename = glue::glue("epi.{.x}.aucplot.pdf"), .plot = .y$plot)
  }
)

# Merge plot --------------------------------------------------------------
merge_plots <- fn_get_merge_plots (
  .list = epi.performance,
  .datasets = as.list(names(epi.performance$panel$performance))
)

purrr::walk2(
  .x = as.list(names(epi.performance$panel$performance)),
  .y = merge_plots,
  .f = function(.x, .y) {
    fn_save_auc(
      .filename = glue::glue("EPI-{.x}-auc-merge.pdf"),
      .plot = .y
    )
  }
)

# Merge metrics -----------------------------------------------------------

merge_metrics <- fn_get_merge_metrics(.list = epi.performance)
readr::write_tsv(x = merge_metrics, file = glue::glue("data/output/EPI-metrics-merge.tsv"))
writexl::write_xlsx(x = merge_metrics, path = glue::glue("data/output/EPI-metrics-merge.xlsx"))

# Save image --------------------------------------------------------------

save.image(file = "data/rda/07-performance-epi.rda")
load("data/rda/07-performance-epi.rda")
