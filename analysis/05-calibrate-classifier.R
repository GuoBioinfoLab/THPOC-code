# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Jul  5 08:59:57 2021
# @DESCRIPTION: 05-calibrate-classifier.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(DESeq2)
library(mlr)

# src ---------------------------------------------------------------------

source(file = "src/calibrate.R")

# Load data ---------------------------------------------------------------

bm.hc.task <- readr::read_rds(file = "data/rda/bm.hc.task.rds.gz")
bm.bam.task <- readr::read_rds(file = "data/rda/bm.bam.task.rds.gz")
model.list <- readr::read_rds(file = "data/rda/model.list.rds.gz")

# Function ----------------------------------------------------------------


# Performance --------------------------------------------------------------

# hc
purrr::map2(
  .x = model.list,
  .y = bm.hc.task,
  .f = fn_get_mlr_pred
) ->
  bm.hc.pred
names(bm.hc.pred) <- c("panel", "panel_ca125", "ca125")

# bam
purrr::map2(
  .x = model.list,
  .y = bm.bam.task,
  .f = fn_get_mlr_pred
) ->
  bm.bam.pred
names(bm.bam.pred) <- c("panel", "panel_ca125", "ca125")


# Calibration -------------------------------------------------------------

# hc
purrr::map(
  .x = bm.hc.pred$panel,
  .f = fn_get_calibrate
) ->
  bm.hc.calibrate

# bam
purrr::map(
  .x = bm.bam.pred$panel,
  .f = fn_get_calibrate
) ->
  bm.bam.calibrate

# Calibration plot --------------------------------------------------------

# hc
tibble::tibble(
  name = list(names(bm.hc.calibrate)),
  calib = list(c("before", "after"))
) %>%
  tidyr::unnest(cols = name) %>%
  tidyr::unnest(cols = calib) %>%
  dplyr::mutate(calib_plot = purrr::map2(
    .x = name,
    .y = calib,
    .f = function(.x, .y) {
      .title <- glue::glue("{.x} {.y} calibration")
      .d <- bm.hc.calibrate[[.x]][[.y]]
      fn_plot_calibration_curve(.x = .d, .title = .title)
    }
  )) ->
  bm.hc.calibrate.plot

# bam
tibble::tibble(
  name = list(names(bm.bam.calibrate)),
  calib = list(c("before", "after"))
) %>%
  tidyr::unnest(cols = name) %>%
  tidyr::unnest(cols = calib) %>%
  dplyr::mutate(calib_plot = purrr::map2(
    .x = name,
    .y = calib,
    .f = function(.x, .y) {
      .title <- glue::glue("{.x} {.y} calibration")
      .d <- bm.hc.calibrate[[.x]][[.y]]
      fn_plot_calibration_curve(.x = .d, .title = .title)
    }
  )) ->
  bm.bam.calibrate.plot

# Save ggplot -------------------------------------------------------------

# hc
purrr::pwalk(
  .l = bm.hc.calibrate.plot,
  .f = function(name, calib, calib_plot) {
    .filename <- glue::glue("BM.HC.{name}-{calib}-calibration-plot.pdf")

    ggsave(
      filename = .filename,
      plot = calib_plot,
      device = "pdf",
      path = "data/reviseoutput/01-model-calibration",
      width = 5.2,
      height = 5
    )
  }
)

# bam
purrr::pwalk(
  .l = bm.bam.calibrate.plot,
  .f = function(name, calib, calib_plot) {
    .filename <- glue::glue("BM.BAM.{name}-{calib}-calibration-plot.pdf")

    ggsave(
      filename = .filename,
      plot = calib_plot,
      device = "pdf",
      path = "data/reviseoutput/01-model-calibration",
      width = 5.2,
      height = 5
    )
  }
)

# Save image --------------------------------------------------------------

save.image(file = "data/rda/05-calibrate-classifier.rda")

