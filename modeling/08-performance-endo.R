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


Endometriosis_samples <- readxl::read_xlsx(path = 'data/raw/wuhan/endoinfo-V5-metadata-1.xlsx', sheet = 3) %>%
  dplyr::pull(barcode)
# Function ----------------------------------------------------------------
fn_get_endo_task <- function(.list, .w, .endo) {
  # .list = c("panel" = list(wuhan.tom.panel.task), wuhan.tom.panel.ca125.task)
  # .w = wuhan.se
  # .endo = Endometriosis
  #

  .list.panel.samples <- .list$panel$samples %>%
    purrr::reduce(.f = c)
  .list$panel$samples <- list(Endometriosis = .list.panel.samples[.endo])

  .list.ca125.samples <- .list$ca125$samples %>%
    purrr::reduce(.f = c)
  .list$ca125$samples <- list(Endometriosis = .list.ca125.samples[.endo])

  .list.panel_ca125.samples <- .list$panel_ca125$samples %>%
    purrr::reduce(.f = c)
  .list$panel_ca125$samples <- list(Endometriosis = .list.panel_ca125.samples[.endo])

  .list
}

fn_performance_endo <- function(.x, .y) {

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

endo.task <- fn_get_endo_task(
  .list = c("panel" = list(wuhan.tom.panel.task), wuhan.tom.panel.ca125.task),
  .w = wuhan.se,
  .endo = Endometriosis_samples
)
readr::write_rds(x = endo.task, file = "data/rda/endo.task.rds.gz", compress = "gz")



# Performance -------------------------------------------------------------


purrr::map2(
  .x = list(panel.model, panel_ca125.model, ca125.model),
  .y = endo.task,
  .f = fn_performance_endo
) ->
  endo.performance

names(endo.performance) <- c("panel", "panel_ca125", "ca125")
readr::write_rds(x = endo.performance, file = "data/rda/endo.performance.rds.gz", compress = "gz")

# save metrics and plot
purrr::walk2(
  .x = list("panel", "panel_ca125", "ca125"),
  .y = endo.performance,
  .f = function(.x, .y) {
    readr::write_tsv(x = .y$metrics, file = glue::glue("data/output/endo.{.x}.metrics.tsv"))
    writexl::write_xlsx(x = .y$metrics, path = glue::glue("data/output/endo.{.x}.metrics.xlsx"))
    fn_save_auc(.filename = glue::glue("endo.{.x}.aucplot.pdf"), .plot = .y$plot)
  }
)

# Merge plot --------------------------------------------------------------
merge_plots <- fn_get_merge_plots (
  .list = endo.performance,
  .datasets = as.list(names(endo.performance$panel$performance))
)

purrr::walk2(
  .x = as.list(names(endo.performance$panel$performance)),
  .y = merge_plots,
  .f = function(.x, .y) {
    fn_save_auc(
      .filename = glue::glue("Endo-{.x}-auc-merge.pdf"),
      .plot = .y
    )
  }
)

# Save image --------------------------------------------------------------

save.image(file = "data/rda/08-performance-endo.rda")
