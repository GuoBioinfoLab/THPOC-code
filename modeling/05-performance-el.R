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

fn_get_el_task <- function(.list, .w) {
  .w@colData %>%
    as.data.frame() %>%
    dplyr::filter(oc %in% c("OC44", "OC79", "OC172")) %>%
    dplyr::filter(!is.na(stage)) %>%
    dplyr::select(barcode, CA125, stage) ->
    .wd

  .wd %>%
    dplyr::filter(stage == "B") %>%
    dplyr::pull(barcode) ->
    .benign
  .wd %>%
    dplyr::filter(stage == "E") %>%
    dplyr::pull(barcode) ->
    .early
  .wd %>%
    dplyr::filter(stage == "L") %>%
    dplyr::pull(barcode) ->
    .late

  .benign_early <- c(.benign, .early)
  .benign_late <- c(.benign, .late)

  .list.panel.samples <- .list$panel$samples %>% purrr::reduce(.f = c)
  .list$panel$samples <- list(
    early = .list.panel.samples[.benign_early],
    late = .list.panel.samples[.benign_late]
  )

  .list.ca125.samples <- .list$ca125$samples %>% purrr::reduce(c)
  .list$ca125$samples <- list(
    early = .list.ca125.samples[.benign_early],
    late = .list.ca125.samples[.benign_late]
  )
  .list$panel_ca125$samples <- list(
    early = .list.ca125.samples[.benign_early],
    late = .list.ca125.samples[.benign_late]
  )

  .list
}


# Prepare tasks -----------------------------------------------------------
el.task <- fn_get_el_task(
  .list = c("panel" = list(wuhan.tom.panel.task), wuhan.tom.panel.ca125.task),
  .w = wuhan.se
)
readr::write_rds(x = el.panel.task, file = "data/rda/el.task.rds.gz", compress = "gz")

# Performance -------------------------------------------------------------

el.performance <- fn_performance(.model = panel.model, .list = el.task$panel)
el.metrics <- fn_get_metrics(.perf = el.performance)
el.plot <- fn_get_auc_plot(.perf = el.performance, .metrics = el.metrics)

el.performance <- fn_performance(.model = ca125.model, .list = el.task$ca125)
el.metrics <- fn_get_metrics(.perf = el.performance)
el.plot <- fn_get_auc_plot(.perf = el.performance, .metrics = el.metrics)

# Save image --------------------------------------------------------------

save.image(file = "data/rda/05-performance-el.rda")
