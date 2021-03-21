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
fn_performance_bm <- function(.x, .y) {

  .perf <- fn_performance(.x, .y)
  .metrics <- fn_get_metrics(.perf)
  .plot <- fn_get_auc_plot(.perf, .metrics)

  list(
    performance = .perf,
    metrics = .metrics,
    plot = .plot
  )

}


# Task --------------------------------------------------------------------
bm.task <- c("panel" = list(wuhan.tom.panel.task), wuhan.tom.panel.ca125.task)



# Performance --------------------------------------------------------------

purrr::map2(
  .x = list(panel.model, panel_ca125.model, ca125.model),
  .y = bm.task,
  .f = fn_performance_bm
) ->
  bm.performance
names(bm.performance) <- c("panel", "panel_ca125", "ca125")
readr::write_rds(x = bm.performance, file = "data/rda/bm.performance.rds.gz", compress = "gz")

# save metrics and plot
purrr::walk2(
  .x = list("panel", "panel_ca125", "ca125"),
  .y = bm.performance,
  .f = function(.x, .y) {
    readr::write_tsv(x = .y$metrics, file = glue::glue("data/output/bm.{.x}.metrics.tsv"))
    writexl::write_xlsx(x = .y$metrics, path = glue::glue("data/output/bm.{.x}.metrics.xlsx"))
    fn_save_auc(.filename = glue::glue("bm.{.x}.aucplot.pdf"), .plot = .y$plot)
  }
)

# Merge plot --------------------------------------------------------------
merge_plots <- fn_get_merge_plots(
  .list = bm.performance,
  .datasets = as.list(names(bm.performance$panel$performance))
)

purrr::walk2(
  .x = as.list(names(bm.performance$panel$performance)),
  .y = merge_plots,
  .f = function(.x, .y) {
    fn_save_auc(
      .filename = glue::glue("BM-{.x}-auc-merge.pdf"),
      .plot = .y
    )
  }
)

# Merge metrics -----------------------------------------------------------

merge_metrics <- fn_get_merge_metrics(.list = bm.performance)
readr::write_tsv(x = merge_metrics, file = glue::glue("data/output/BM-metrics-merge.tsv"))
writexl::write_xlsx(x = merge_metrics, path = glue::glue("data/output/BM-metrics-merge.xlsx"))


# Save image --------------------------------------------------------------

save.image(file = "data/rda/05-performance-bm.rda")
load(file = "data/rda/05-performance-bm.rda")
