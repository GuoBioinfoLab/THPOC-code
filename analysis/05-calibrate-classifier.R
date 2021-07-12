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

model.list <- readr::read_rds(file = "data/rda/model.list.rds.gz")

bm.hc.task <- readr::read_rds(file = "data/rda/bm.hc.task.rds.gz")
bm.bam.task <- readr::read_rds(file = "data/rda/bm.bam.task.rds.gz")
bm.hc.bam.task <- list(bm.hc = bm.hc.task, bm.bam = bm.bam.task)

el.hc.bam.task <- readr::read_rds(file = "data/rda/el.hc.bam.task.rds.gz")
epi.hc.bam.task <- readr::read_rds(file = "data/rda/epi.hc.bam.task.rds.gz")
endo.hc.bam.task <- readr::read_rds(file = "data/rda/endo.hc.bam.task.rds.gz")
bl.hc.bam.task <- readr::read_rds(file = "data/rda/bl.hc.bam.task.rds.gz")
glh.hc.bam.task <- readr::write_rds(file = "data/rda/glh.hc.bam.task.rds.gz")
sr.hc.bam.task <- readr::write_rds(file = "data/rda/sr.hc.bam.task.rds.gz")

# Function ----------------------------------------------------------------

fn_calibrate_pipe <- function(.x, .prefix) {
  .prefix <- toupper(.prefix)
  # predict
  .pred <- purrr::map2(.x = model.list, .y = .x, .f = fn_get_mlr_pred)
  # calibration
  .calibrate <- purrr::map(.x = .pred$panel, .f = fn_get_calibrate )
  # plot
  .calibrate.plot <- tibble::tibble(
    name = list(names(.calibrate)),
    calib = list(c("before", "after"))
  ) %>%
    tidyr::unnest(cols = name) %>%
    tidyr::unnest(cols = calib) %>%
    dplyr::mutate(calib_plot = purrr::map2(
      .x = name,
      .y = calib,
      .f = function(.x, .y) {
        .title <- glue::glue("{.x} {.y} calibration")
        .d <- .calibrate[[.x]][[.y]]
        fn_plot_calibration_curve(.x = .d, .title = .title)
      }
    ))
  # save plot
  purrr::pwalk(
    .l = .calibrate.plot,
    .f = function(name, calib, calib_plot) {
      .filename <- glue::glue("{.prefix}.{name}-{calib}-calibration-plot.pdf")

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
}

# Calibrate ---------------------------------------------------------------

list(
  bm.hc.bam.task,
  el.hc.bam.task,
  epi.hc.bam.task,
  endo.hc.bam.task,
  bl.hc.bam.task,
  glh.hc.bam.task,
  sr.hc.bam.task
) %>%
  purrr::walk(.f = function(.l) {
    names(.l) %>%
      purrr::walk(.f = function(.name) {
        fn_calibrate_pipe(.x = .l[[.name]], .prefix = .name)
      })
  })


# Save image --------------------------------------------------------------

save.image(file = "data/rda/05-calibrate-classifier.rda")

