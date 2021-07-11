# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Jul 11 08:46:06 2021
# @DESCRIPTION: 06-performance-new-bm.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(DESeq2)
library(mlr)

# src ---------------------------------------------------------------------

source(file = "src/utils.R")
source(file = "src/performance.R")
source(file = "src/plots.R")


# Load data ---------------------------------------------------------------

wuhan.se <- readr::read_rds(file = "data/rda/wuhan.se.rds.gz")
tom.se <- readr::read_rds(file = "data/rda/tom.se.rds.gz")
wuhan.tom.fs.fg.norm.rbe.se <- readr::read_rds(file = "data/rda/wuhan.tom.fs.fg.norm.rbe.se.rds.gz")

bm.task <- readr::read_rds(file = "data/rda/bm.task.rds.gz")

panel <- readr::read_rds(file = "data/rda/panel.rds.gz")

panel.model <- readr::read_rds(file = "data/rda/panel.model.rds.gz")
panel_ca125.model <- readr::read_rds(file = "data/rda/panel_ca125.model.rds.gz")
ca125.model <- readr::read_rds(file = "data/rda/ca125.model.rds.gz")
model.list <- readr::read_rds(file = "data/rda/model.list.rds.gz")

# Function ----------------------------------------------------------------

fn_get_bm_hc_task <- function(.task) {
  purrr::map(.x = .task, .f = function(.x) {
    .vc1_vc2 <- c(.x$samples$VC1, .x$samples$VC2)
    .vc1_vc2_tom <- c(.vc1_vc2, .x$samples$Tom)
    .x$samples$VC1_VC2 <- .vc1_vc2
    .x$samples$VC1_VC2_Tom <- .vc1_vc2_tom

    .x
  })
}

fn_get_bm_bam_task <- function(){}

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

bm.hc.task <- fn_get_bm_hc_task(.task = bm.task)
readr::write_rds(
  x = bm.hc.task,
  file = "data/rda/bm.hc.task.rds.gz",
  compress = "gz"
)

# Performance -------------------------------------------------------------

purrr::map2(
  .x = model.list,
  .y = bm.hc.task,
  .f = fn_performance_bm
) ->
  bm.hc.performance
names(bm.hc.performance) <- c("panel", "panel_ca125", "ca125")
readr::write_rds(
  x = bm.hc.performance,
  file = "data/rda/bm.hc.performance.rds.gz",
  compress = "gz"
)


# save metrics and plot
purrr::walk2(
  .x = list("panel", "panel_ca125", "ca125"),
  .y = bm.hc.performance,
  .f = function(.x, .y) {
    readr::write_tsv(
      x = .y$metrics,
      file = glue::glue("data/reviseoutput/02-BM/bm.hc.{.x}.metrics.tsv")
    )
    writexl::write_xlsx(
      x = .y$metrics,
      path = glue::glue("data/reviseoutput/02-BM/bm.hc.{.x}.metrics.xlsx")
    )
    ggsave(
      filename =  glue::glue("bm.hc.{.x}.aucplot.pdf"),
      plot = .y$plot,
      device = "pdf",
      path = "data/reviseoutput/02-BM",
      width = 5.2,
      height = 4.5
    )
  }
)
