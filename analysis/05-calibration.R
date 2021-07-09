# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Jul  5 08:59:57 2021
# @DESCRIPTION: 05-calibration.R

# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)
library(mlr)
library(doParallel)
library(ggplot2)

# Load data ---------------------------------------------------------------

wuhan.tom.panel.task <- readr::read_rds(file = "data/rda/wuhan.tom.panel.task.rds.gz")
wuhan.tom.panel.ca125.task <- readr::read_rds(file = "data/rda/wuhan.tom.panel.ca125.task.rds.gz")

panel.model <- readr::read_rds(file = "data/rda/panel.model.rds.gz")
ca125.model <- readr::read_rds(file = "data/rda/ca125.model.rds.gz")
panel_ca125.model <- readr::read_rds(file = "data/rda/panel_ca125.model.rds.gz")


# Function ----------------------------------------------------------------

fn_predict <- function(.x, .y) {
  .model <- .x
  .task <- .y$task
  .samples <- .y$samples

  purrr::map(
    .x = names(.samples),
    .f = function(.x, .model, .task, .samples = .samples) {
      .pred <- predict(
        object = .model,
        task = .task,
        subset = .samples[[.x]]
      )
      .pred
    },
    .model = .model,
    .task = .task,
    .samples = .samples
  ) ->
    .perf

  names(.perf) <- names(.samples)
  .perf
}


# Task --------------------------------------------------------------------
bm.task <- c("panel" = list(wuhan.tom.panel.task), wuhan.tom.panel.ca125.task)


# Performance --------------------------------------------------------------

purrr::map2(
  .x = list(panel.model, panel_ca125.model, ca125.model),
  .y = bm.task,
  .f = fn_predict
) ->
  bm.predict

names(bm.predict) <- c("panel", "panel_ca125", "ca125")


# Calibration -------------------------------------------------------------


performance(
  pred = bm.predict$panel$TC,
  measures = list(auc, brier)
)
cal <- generateCalibrationData(bm.predict$panel$TC)
cal$proportion
plotCalibration(cal, smooth = TRUE)

