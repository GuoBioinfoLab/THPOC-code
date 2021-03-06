# Meta info -----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-03-1 10:00:58
# @DESCRIPTION: 06-performance-el.R
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

fn_get_el_task <- function(.list, .w, .t) {
  .w@colData %>%
    as.data.frame() %>%
    dplyr::filter(oc %in% c("OC44", "OC79", "OC172")) %>%
    dplyr::mutate(stage = ifelse(type == "normal", "H", stage)) %>%
    dplyr::select(barcode, stage) ->
    .wd

  .t@colData %>%
    as.data.frame() %>%
    dplyr::select(barcode, stageFourGroups, Stage) %>%
    dplyr::mutate(stage = plyr::revalue(x = stageFourGroups, replace = c(
      "benign" = "B",
      "healthy control" = "H",
      "I" = "E",
      "II" = "E",
      "III" = "L",
      "IV" = "L"
    ))) %>%
    dplyr::mutate(stage = ifelse(Stage == "IIB", "L", stage)) %>%
    dplyr::select(barcode, stage) ->
    .td

  dplyr::bind_rows(.wd, .td) -> .wtd

  .wtd %>%
    dplyr::filter(stage == "B") %>%
    dplyr::pull(barcode) ->
    .benign
  .wtd %>%
    dplyr::filter(stage == "E") %>%
    dplyr::pull(barcode) ->
    .early
  .wtd %>%
    dplyr::filter(stage == "L") %>%
    dplyr::pull(barcode) ->
    .late

  .benign_early <- c(.benign, .early)
  .benign_late <- c(.benign, .late)

  .list.panel.samples <- .list$panel$samples %>% purrr::reduce(.f = c)
  .list$panel$samples <- list(
    Early = .list.panel.samples[.benign_early] %>% na.omit() %>% c(),
    Late = .list.panel.samples[.benign_late] %>% na.omit() %>% c()
  )
  .list.ca125.samples <- .list$ca125$samples %>% purrr::reduce(c)
  .list$ca125$samples <- list(
    Early = .list.ca125.samples[.benign_early] %>% na.omit() %>% c(),
    Late = .list.ca125.samples[.benign_late] %>% na.omit() %>% c()
  )
  .list$panel_ca125$samples <- list(
    Early = .list.ca125.samples[.benign_early] %>% na.omit() %>% c(),
    Late = .list.ca125.samples[.benign_late] %>% na.omit() %>% c()
  )

  .list
}

fn_performance_el <- function(.x, .y) {

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


el.task <- fn_get_el_task(
  .list = c("panel" = list(wuhan.tom.panel.task), wuhan.tom.panel.ca125.task),
  .w = wuhan.se,
  .t = tom.se
)
readr::write_rds(x = el.task, file = "data/rda/el.task.rds.gz", compress = "gz")

# Performance -------------------------------------------------------------

purrr::map2(
  .x = list(panel.model, panel_ca125.model, ca125.model),
  .y = el.task,
  .f = fn_performance_el
) ->
  el.performance

names(el.performance) <- c("panel", "panel_ca125", "ca125")
readr::write_rds(x = el.performance, file = "data/rda/el.performance.rds.gz", compress = "gz")

# save metrics and plot
purrr::walk2(
  .x = list("panel", "panel_ca125", "ca125"),
  .y = el.performance,
  .f = function(.x, .y) {
    readr::write_tsv(x = .y$metrics, file = glue::glue("data/output/el.{.x}.metrics.tsv"))
    writexl::write_xlsx(x = .y$metrics, path = glue::glue("data/output/el.{.x}.metrics.xlsx"))
    fn_save_auc(.filename = glue::glue("el.{.x}.aucplot.pdf"), .plot = .y$plot)
  }
)

# Merge plot --------------------------------------------------------------
merge_plots <- fn_get_merge_plots (
  .list = el.performance,
  .datasets = as.list(names(el.performance$panel$performance))
  )

purrr::walk2(
  .x = as.list(names(el.performance$panel$performance)),
  .y = merge_plots,
  .f = function(.x, .y) {
    fn_save_auc(
      .filename = glue::glue("EL-{.x}-auc-merge.pdf"),
      .plot = .y
    )
  }
)

# Merge metrics -----------------------------------------------------------

merge_metrics <- fn_get_merge_metrics(.list = el.performance)
readr::write_tsv(x = merge_metrics, file = glue::glue("data/output/EL-metrics-merge.tsv"))
writexl::write_xlsx(x = merge_metrics, path = glue::glue("data/output/EL-metrics-merge.xlsx"))

# Save image --------------------------------------------------------------

save.image(file = "data/rda/06-performance-el.rda")
load(file = "data/rda/06-performance-el.rda")
