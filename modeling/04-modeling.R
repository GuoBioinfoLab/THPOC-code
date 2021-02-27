
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

# Function ----------------------------------------------------------------
fn_get_task_panel <- function(.wt, .panel) {
  .se <- .wt[.panel, ]

  .task <- fn_se2task_panel(.se = .se, .id = "Panel-task")

  .datasets <- list("TC", "DC", "VC1", "VC2", "Tom")
  .samples <- .datasets %>% purrr::map(.f = fn_task_ind, .se = .se)
  names(.samples) <- .datasets

  list(
    task = .task,
    samples = .samples
  )
}

fn_get_task_panel_ca125 <- function(.t, .w, .wt, .panel) {
  # .w <- wuhan.se
  # .wt <- wuhan.tom.fs.fg.norm.rbe.se
  # .panel <- panel
  .w@colData %>%
    as.data.frame() %>%
    dplyr::filter(barcode %in% .wt$barcode) %>%
    dplyr::group_by(oc) %>%
    dplyr::mutate(ca125 = scale(log2(CA125))[, 1]) %>%
    dplyr::ungroup() ->
    .wd

  .t@colData %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::filter(barcode %in% .wt$barcode) %>%
    dplyr::mutate(ca125 = as.numeric(CA125parameterTOC)) %>%
    dplyr::filter(!is.na(ca125), ca125 > 0) %>%
    dplyr::mutate(ca125 = scale(log2(ca125))[, 1]) ->
    .td

  .se <- .wt[.panel, c(.wd$barcode, .td$barcode)]
  .se@colData$CA125 <- c(.wd$ca125, .td$ca125)

  .task_list <- fn_se2task_panel_ca125(.se = .se)
  .datasets <- list("TC", "DC", "VC1", "VC2", "Tom")
  .samples <- .datasets %>%
    purrr::map(.f = fn_task_ind, .se)
  names(.samples) <- .datasets

  list(
    panel_ca125 = list(
      task = .task_list$task_panel_ca125,
      samples = .samples
    ) ,
    ca125 = list(
      task = .task_list$task_ca125,
      samples = .samples
    )
  )

}

fn_tune_hyper_parameters <- function(.list) {
  .task <- .list$task
  .samples <- .list$samples
  .task_id <- mlr::getTaskId(x = .task)

  .task_for_tunes <- mlr::subsetTask(task = .task, subset = .samples$TC)

  .tuned_model <- fn_tune_model(.tsk = .task_for_tunes)

  fn_plot_tune_path(.tune_result = .tuned_model$tune_result, .task_id = .task_id)

  mlr::train(learner = .tuned_model$learner, task = .task, subset = .samples$TC)

}


# Panel task --------------------------------------------------------------

wuhan.tom.panel.task <- fn_get_task_panel(.wt = wuhan.tom.fs.fg.norm.rbe.se, .panel = panel)
readr::write_rds(x = wuhan.tom.panel.task, file = "data/rda/wuhan.tom.panel.task.rds.gz")

# CA125 task --------------------------------------------------------------

wuhan.tom.panel.ca125.task <- fn_get_task_panel_ca125(.t = tom.se, .w = wuhan.se, .wt = wuhan.tom.fs.fg.norm.rbe.se, .panel = panel)
readr::write_rds(x = wuhan.tom.panel.ca125.task, file = "data/rda/wuhan.tom.panel.ca125.task.rds.gz")

# Panel modeling ----------------------------------------------------------


panel.model <- fn_tune_hyper_parameters(.list = wuhan.tom.panel.task)
readr::write_rds(x = panel.model, file = "data/rda/panel.model.rds.gz", compress = "gz")


# CA125 modeling ----------------------------------------------------------

ca125.model <- fn_tune_hyper_parameters(.list = wuhan.tom.panel.ca125.task$ca125)
readr::write_rds(x = ca125.model, file = "data/rda/ca125.model.rds.gz", compress = "gz")


# Panel CA125 modeling ----------------------------------------------------

panel_ca125.model <- fn_tune_hyper_parameters(.list = wuhan.tom.panel.ca125.task$panel_ca125)
readr::write_rds(x = panel_ca125.model, file = "data/rda/panel_ca125.model.rds.gz", compress = "gz")


# Save image --------------------------------------------------------------

save.image(file = "data/rda/04-modeling.rda")


