
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
fn_task_panel <- function(.wt, .panel) {
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

fn_task_panel_ca125 <- function(.w, .wt, .panel) {
  # .w <- wuhan.se
  # .wt <- wuhan.tom.fs.fg.norm.rbe.se
  # .panel <- panel
  .w@colData %>%
    as.data.frame() %>%
    dplyr::mutate(a = ifelse(oc == "OC521", "TC", "Others")) %>%
    dplyr::group_by(a) %>%
    dplyr::mutate(ca125 = scale(log2(CA125))[, 1]) ->
    .wd
  .se <- .wt[.panel, .w$barcode]
  .se@colData$CA125 <- .wd$ca125

  .task_list <- fn_se2task_panel_ca125(.se = .se)
  .datasets <- list("TC", "DC", "VC1", "VC2")
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

fn_get_tom_plot <- function(.perf) {
  .legend_title <- glue::glue("AUC for Tom")
  .legend_text <- .perf %>%
    dplyr::mutate(auc_label = round(auc, digits = 3)) %>%
    dplyr::distinct(auc_label) %>%
    dplyr::pull(auc_label)
  .legend_color <- "#B22222"

  .perf %>%
    ggplot(aes(x = fpr, y = tpr, color = cohort)) +
    geom_path(size = 0.8) +
    scale_x_continuous(
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
      limits = c(0, 1), expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_color_manual(
      name = .legend_title,
      labels = .legend_text,
      values = .legend_color
    ) +
    guides(
      color = guide_legend(
        reverse = TRUE
      )
    ) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    fn_auc_theme()
}

# Panel -------------------------------------------------------------------


wuhan.tom.panel.task <- fn_task_panel(.wt = wuhan.tom.fs.fg.norm.rbe.se, .panel = panel)
readr::write_rds(x = wuhan.tom.panel.task, file = "data/rda/wuhan.tom.panel.task.rds.gz")

panel.model <- fn_tune_hyper_parameters(.list = wuhan.tom.panel.task)
readr::write_rds(x = panel.model, file = "data/rda/panel.model.rds.gz", compress = "gz")

panel.performance <- fn_performance(.model = panel.model, .list = wuhan.tom.panel.task)

panel.metrics <- fn_get_metrics(.perf = panel.performance)
readr::write_tsv(x = panel.metrics, file = "data/output/panel.metrics.tsv")
writexl::write_xlsx(x = panel.metrics, path = "data/output/panel.metrics.xlsx")

panel.plot <- fn_get_auc_plot(.perf = panel.performance, .metrics = panel.metrics)
fn_save_auc(.filename = "panel.aucplot.pdf", .plot = panel.plot)


# CA125 + panel-------------------------------------------------------------------

wuhan.tom.panel.ca125.task <- fn_task_panel_ca125(.w = wuhan.se, .wt = wuhan.tom.fs.fg.norm.rbe.se, .panel = panel)
readr::write_rds(x = wuhan.tom.panel.ca125.task, file = "data/rda/wuhan.tom.panel.ca125.task.rds.gz")
# CA125 -------------------------------------------------------------------


ca125.model <- fn_tune_hyper_parameters(.list = wuhan.tom.panel.ca125.task$ca125)

readr::write_rds(x = ca125.model, file = "data/rda/ca125.model.rds.gz", compress = "gz")

ca125.performance <- fn_performance(.model = ca125.model, .list = wuhan.tom.panel.ca125.task$ca125)

ca125.metrics <- fn_get_metrics(.perf = ca125.performance)
readr::write_tsv(x = ca125.metrics, file = "data/output/ca125.metrics.tsv")
writexl::write_xlsx(x = ca125.metrics, path = "data/output/ca125.metrics.xlsx")

ca125.plot <- fn_get_auc_plot(.perf = ca125.performance, .metrics = ca125.metrics)
fn_save_auc(.filename = "ca125.aucplot.pdf", .plot = ca125.plot)

# Panel CA125 -------------------------------------------------------------


panel_ca125.model <- fn_tune_hyper_parameters(.list = wuhan.tom.panel.ca125.task$panel_ca125)

readr::write_rds(x = panel_ca125.model, file = "data/rda/panel_ca125.model.rds.gz", compress = "gz")

panel_ca125.performance <- fn_performance(.model = panel_ca125.model, .list = wuhan.tom.panel.ca125.task$panel_ca125)

panel_ca125.metrics <- fn_get_metrics(.perf = panel_ca125.performance)
readr::write_tsv(x = panel_ca125.metrics, file = "data/output/panel_ca125.metrics.tsv")
writexl::write_xlsx(x = panel_ca125.metrics, path = "data/output/panel_ca125.metrics.xlsx")

panel_ca125.plot <- fn_get_auc_plot(.perf = panel_ca125.performance, .metrics = panel_ca125.metrics)
fn_save_auc(.filename = "panel_ca125.aucplot.pdf", .plot = panel_ca125.plot)


# Merge plot --------------------------------------------------------------

merge_plots <- fn_get_merge_plots(
  .list = list(
    "panel" = list(panel.performance),
    "ca125" = list(ca125.performance),
    "panel_ca125" = list(panel_ca125.performance)
  ),
  .datasets = list("TC", "DC", "VC1", "VC2")
  )

purrr::walk2(
  .x = list("TC", "DC", "VC1", "VC2"),
  .y = merge_plots,
  .f = function(.x, .y) {
    fn_save_auc(
      .filename = glue::glue("BM-{.x}-auc-merge.pdf"),
      .plot = .y
    )
  }
)

# Plot tom ----------------------------------------------------------------
tom_plot <- fn_get_tom_plot(panel.performance$Tom$perf)
fn_save_auc(
  .filename = glue::glue("BM-Tom-auc-merge.pdf"),
  .plot = tom_plot
)

# Save image --------------------------------------------------------------

save.image(file = "data/rda/04-modeling.rda")


