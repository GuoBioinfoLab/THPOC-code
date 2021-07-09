# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Jul  5 08:59:57 2021
# @DESCRIPTION: 05-calibrate-classifier.R

# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)
library(mlr)
library(doParallel)
library(ggplot2)

# src ---------------------------------------------------------------------

source(file = "src/calibrate.R")

# Load data ---------------------------------------------------------------

wuhan.tom.panel.task <- readr::read_rds(file = "data/rda/wuhan.tom.panel.task.rds.gz")
wuhan.tom.panel.ca125.task <- readr::read_rds(file = "data/rda/wuhan.tom.panel.ca125.task.rds.gz")

panel.model <- readr::read_rds(file = "data/rda/panel.model.rds.gz")
ca125.model <- readr::read_rds(file = "data/rda/ca125.model.rds.gz")
panel_ca125.model <- readr::read_rds(file = "data/rda/panel_ca125.model.rds.gz")


# Function ----------------------------------------------------------------



# Task --------------------------------------------------------------------
bm.task <- c("panel" = list(wuhan.tom.panel.task), wuhan.tom.panel.ca125.task)


# Performance --------------------------------------------------------------

purrr::map2(
  .x = list(panel.model, panel_ca125.model, ca125.model),
  .y = bm.task,
  .f = fn_get_mlr_pred
) ->
  bm.pred

names(bm.pred) <- c("panel", "panel_ca125", "ca125")


# Calibration -------------------------------------------------------------

purrr::map(
  .x = bm.pred$panel,
  .f = fn_get_calibrate
) ->
  bm.calibrate


bm.calibrate$TC$before_calib$pred %>%
  generateCalibrationData() %>%
  plotCalibration(smooth = TRUE)

bm.calibrate$VC1$before_calib$pred %>%
  generateCalibrationData() %>%
  plotCalibration(smooth = TRUE)

bm.calibrate$VC1$before_cali$auc_brier

bm.calibrate$VC2$before_calib$pred %>%
  generateCalibrationData() %>%
  plotCalibration(smooth = TRUE)
bm.calibrate$VC1$after_calib$auc_brier

