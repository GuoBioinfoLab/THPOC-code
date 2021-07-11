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

bm.task <- readr::read_rds(file = "data/rda/bm.task.rds.gz")
panel.model <- readr::read_rds(file = "data/rda/panel.model.rds.gz")
ca125.model <- readr::read_rds(file = "data/rda/ca125.model.rds.gz")
panel_ca125.model <- readr::read_rds(file = "data/rda/panel_ca125.model.rds.gz")


# Function ----------------------------------------------------------------

fn_get_task <- functon(.w, .t, .wt, .panel) {
  .se <- .wt[.panel, ]
  .task <- fn_se2task_panel(.se = .se, .id = "Panel-task")

  .datasets <- list("TC", "DC", "VC1", "VC2", "Tom", "VC1_VC2", "VC1_VC2_Tom")

  .w[.panel,]@colData %>%
    as.data.frame() %>%
}

# Task --------------------------------------------------------------------

