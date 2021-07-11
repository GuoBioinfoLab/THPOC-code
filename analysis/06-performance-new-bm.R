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

fn_get_bm_bam_task <- function(.task, .w, .t){
  .w@colData %>%
    as.data.frame() %>%
    dplyr::filter(type != "normal") %>%
    dplyr::pull(barcode) ->
    .ws

  .t@colData %>%
    as.data.frame() %>%
    dplyr::filter(patientGroup != "HC") %>%
    dplyr::pull(barcode) ->
    .ts

  .wts <- c(.ws, .ts)

  purrr::map(.x = .task, .f = function(.x) {
    .samples <- purrr::map(.x$samples, .f = function(.x) {c(na.omit(.x[.wts]))})
    .x$samples <- .samples
    .x
  })
}

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

fn_save_auc <- function(.filename, .plot) {
  ggsave(
    filename = .filename,
    plot = .plot,
    device = "pdf",
    path = "data/reviseoutput/02-BM",
    width = 5.2,
    height = 4.5
  )
}

# Task --------------------------------------------------------------------

# hc task
bm.hc.task <- fn_get_bm_hc_task(.task = bm.task)
readr::write_rds(
  x = bm.hc.task,
  file = "data/rda/bm.hc.task.rds.gz",
  compress = "gz"
)

# bam task
bm.bam.task <- fn_get_bm_bam_task(.task = bm.hc.task, .w = wuhan.se, .t = tom.se)
readr::write_rds(
  x = bm.bam.task,
  file = "data/rda/bm.bam.task.rds.gz",
  compress = "gz"
)

# Performance -------------------------------------------------------------

# hc
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
    fn_save_auc(
      .filename = glue::glue("bm.hc.{.x}.aucplot.pdf"),
      .plot = .y$plot
    )
  }
)

# bam
purrr::map2(
  .x = model.list,
  .y = bm.bam.task,
  .f = fn_performance_bm
) ->
  bm.bam.performance
names(bm.bam.performance) <- c("panel", "panel_ca125", "ca125")
readr::write_rds(
  x = bm.bam.performance,
  file = "data/rda/bm.bam.performance.rds.gz",
  compress = "gz"
)

purrr::walk2(
  .x = list("panel", "panel_ca125", "ca125"),
  .y = bm.bam.performance,
  .f = function(.x, .y) {
    readr::write_tsv(
      x = .y$metrics,
      file = glue::glue("data/reviseoutput/02-BM/bm.bam.{.x}.metrics.tsv")
    )
    writexl::write_xlsx(
      x = .y$metrics,
      path = glue::glue("data/reviseoutput/02-BM/bm.bam.{.x}.metrics.xlsx")
    )
    fn_save_auc(
      .filename = glue::glue("bm.bam.{.x}.aucplot.pdf"),
      .plot = .y$plot
    )
  }
)

# Merge plot --------------------------------------------------------------

# hc
merge_plots_hc <- fn_get_merge_plots(
  .list = bm.hc.performance,
  .datasets = as.list(names(bm.hc.performance$panel$performance))
)

purrr::walk2(
  .x = as.list(names(bm.hc.performance$panel$performance)),
  .y = merge_plots_hc,
  .f = function(.x, .y) {
    fn_save_auc(
      .filename = glue::glue("BM.HC-{.x}-auc-merge.pdf"),
      .plot = .y
    )
  }
)

# bam

merge_plots_bam <- fn_get_merge_plots(
  .list = bm.bam.performance,
  .datasets = as.list(names(bm.bam.performance$panel$performance))
)

purrr::walk2(
  .x = as.list(names(bm.bam.performance$panel$performance)),
  .y = merge_plots_bam,
  .f = function(.x, .y) {
    fn_save_auc(
      .filename = glue::glue("BM.BAM-{.x}-auc-merge.pdf"),
      .plot = .y
    )
  }
)

# Merge metrics -----------------------------------------------------------

# hc
merge_metrics_hc <- fn_get_merge_metrics(.list = bm.hc.performance)
readr::write_tsv(x = merge_metrics_hc, file = glue::glue("data/reviseoutput/02-BM/BM.HC-metrics-merge.tsv"))
writexl::write_xlsx(x = merge_metrics_hc, path = glue::glue("data/reviseoutput/02-BM/BM.HC-metrics-merge.xlsx"))

# bam

merge_metrics_bam <- fn_get_merge_metrics(.list = bm.bam.performance)
readr::write_tsv(x = merge_metrics_bam, file = glue::glue("data/reviseoutput/02-BM/BM.BAM-metrics-merge.tsv"))
writexl::write_xlsx(x = merge_metrics_bam, path = glue::glue("data/reviseoutput/02-BM/BM.BAM-metrics-merge.xlsx"))

# Save image --------------------------------------------------------------

save.image(file = "data/rda/06-performance-new-bm.rda")
load(file = "data/rda/06-performance-new-bm.rda")
