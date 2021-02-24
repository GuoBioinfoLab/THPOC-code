
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

# Load data ---------------------------------------------------------------

wuhan.tom.fs.fg.norm.rbe.se <- readr::read_rds(file = "data/rda/wuhan.tom.fs.fg.norm.rbe.se.rds.gz")
panel <- readr::read_rds(file = "data/rda/panel.rds.gz")

# Function ----------------------------------------------------------------
fn_task <- function(.wt, .panel) {
  .se <- .wt[.panel, ]

  .task <- fn_bm_se2task(.se = .se, .id = "Panel-task")

  .datasets <- list("TC", "DC", "VC1", "VC2", "Tom")
  .samples <- .datasets %>% purrr::map(.f = fn_task_ind, .se = .se)
  names(.samples) <- .datasets

  list(
    task = .task,
    samples = .samples
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

fn_performance <- function(.model, .list) {
  .task <- .list$task
  .samples <- .list$samples

  purrr::map(
    names(.samples),
    fn_predidct_performance_metrics,
    .model = .model,
    .task = .task,
    .samples = .samples
  ) ->
    .perf

  names(.perf) <- names(.samples)
  .perf
}

fn_get_metrics <- function(.perf) {
  .perf %>%
    purrr::map("metrics") %>%
    purrr::reduce(.f = dplyr::bind_rows)
}

fn_get_auc_plot <- function(.perf, .metrics) {

  .metrics %>%
    dplyr::mutate(auc = gsub(pattern = " ", replacement = "", x = `AUC (95% CI)`)) %>%
    dplyr::select(cohort, auc) %>%
    dplyr::mutate(label = glue::glue("{cohort} {auc}")) ->
    .labels

  .d <- .perf %>%
    purrr::map("perf") %>%
    purrr::reduce(.f = dplyr::bind_rows) %>%
    dplyr::mutate(cohort = factor(x = cohort, levels = .labels$cohort))

  fn_plot_auc(.d, .labels)

}
# Prepare task ------------------------------------------------------------
wuhan.tom.task <- fn_task(.wt = wuhan.tom.fs.fg.norm.rbe.se, .panel = panel)

panel.model <- fn_tune_hyper_parameters(.list = wuhan.tom.task)
readr::write_rds(x = panel.model, file = "data/rda/panel.model.rds.gz", compress = "gz")

panel.performance <- fn_performance(.model = panel.model, .list = wuhan.tom.task)

panel.metrics <- fn_get_metrics(.perf = panel.performance)
readr::write_tsv(x = panel.metrics, file = "data/output/panel.metrics.tsv")
writexl::write_xlsx(x = panel.metrics, path = "data/output/panel.metrics.xlsx")

panel.plot <- fn_get_auc_plot(.perf = panel.performance, .metrics = panel.metrics)
ggsave(
  filename = "panel.aucplot.pdf",
  plot = panel.plot,
  path = "data/output",
  device = "pdf",
  width = 6,
  height = 6
)

