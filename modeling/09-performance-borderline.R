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

wuhan.tom.panel.task <- readr::read_rds(file = "data/rda/wuhan.tom.panel.task.rds.gz")
wuhan.tom.panel.ca125.task <- readr::read_rds(file = "data/rda/wuhan.tom.panel.ca125.task.rds.gz")
panel.model <- readr::read_rds(file = "data/rda/panel.model.rds.gz")
ca125.model <- readr::read_rds(file = "data/rda/ca125.model.rds.gz")
panel_ca125.model <- readr::read_rds(file = "data/rda/panel_ca125.model.rds.gz")

borderline <- readxl::read_xlsx(path = "data/raw/borderline.xlsx") %>%
  dplyr::filter(!is.na(grade)) %>%
  dplyr::filter(barcode != "name") %>%
  dplyr::filter(barcode != "barcode") %>%
  dplyr::pull(barcode)

# Function ----------------------------------------------------------------

fn_get_borderline_task <- function(.list, .w ,.t, .bl) {
  # .list = c("panel"= list(wuhan.tom.panel.task), wuhan.tom.panel.ca125.task)
  # .w = wuhan.se
  # .t = tom.se
  # .bl = borderline
  .w@colData %>%
    as.data.frame() %>%
    dplyr::filter(oc %in% c("OC44", "OC79", "OC172")) %>%
    dplyr::mutate(name = barcode) %>%
    dplyr::mutate(class = as.character(class)) %>%
    dplyr::select(barcode, name, class) ->
    .wd

  .t@colData %>%
    as.data.frame() %>%
    dplyr::mutate(class = as.character(class)) %>%
    dplyr::select(barcode, name, class) ->
    .td

  dplyr::bind_rows(.wd, .td) -> .wtd

  .wtd %>%
    dplyr::filter(name %in% .bl) %>%
    dplyr::pull(barcode) ->
    .blbarcode
  .wtd %>%
    dplyr::filter(class == "B") %>%
    dplyr::pull(barcode) ->
    .nblbarcode
  .samples <- c(.blbarcode, .nblbarcode)

  .list.panel.samples <- .list$panel$samples %>% purrr::reduce(.f = c)
  .list$panel$samples <- list(Borderline = .list.panel.samples[.samples] %>% na.omit() %>% c())
  .list.ca125.samples <- .list$ca125$samples %>% purrr::reduce(.f = c)
  .list$ca125$samples <- list(Borderline = .list.ca125.samples[.samples] %>% na.omit() %>% c())
  .list.panel_ca125.samples <- .list$panel_ca125$samples %>% purrr::reduce(.f = c)
  .list$panel_ca125$samples <- list(Borderline = .list.panel_ca125.samples[.samples] %>% na.omit() %>% c())
  .list
}

fn_performance_borderline <- function(.x, .y) {

  .perf <- fn_performance(.x, .y)
  .metrics <- fn_get_metrics(.perf)
  .plot <- fn_get_auc_plot(.perf, .metrics)

  list(
    performance = .perf,
    metrics = .metrics,
    plot = .plot
  )

}
# Prepare tasks -----------------------------------------------------------

borderline.task <- fn_get_borderline_task(
  .list = c("panel"= list(wuhan.tom.panel.task), wuhan.tom.panel.ca125.task),
  .w = wuhan.se,
  .t = tom.se,
  .bl = borderline
)

readr::write_rds(x = borderline.task, file = "data/rda/borderline.task.rds.gz", compress = "gz")


# Performance -------------------------------------------------------------

purrr::map2(
  .x = list(panel.model, panel_ca125.model, ca125.model),
  .y = borderline.task,
  .f = fn_performance_borderline
) ->
  borderline.performance

names(borderline.performance) <- c("panel", "panel_ca125", "ca125")
readr::write_rds(x = borderline.performance, file = "data/rda/borderline.performance.rds.gz", compress = "gz")

# save metrics and plot
purrr::walk2(
  .x = list("panel", "panel_ca125", "ca125"),
  .y = borderline.performance,
  .f = function(.x, .y) {
    readr::write_tsv(x = .y$metrics, file = glue::glue("data/output/borderline.{.x}.metrics.tsv"))
    writexl::write_xlsx(x = .y$metrics, path = glue::glue("data/output/borderline.{.x}.metrics.xlsx"))
    fn_save_auc(.filename = glue::glue("borderline.{.x}.aucplot.pdf"), .plot = .y$plot)
  }
)

# Merge plot --------------------------------------------------------------
merge_plots <- fn_get_merge_plots (
  .list = borderline.performance,
  .datasets = as.list(names(borderline.performance$panel$performance))
)

purrr::walk2(
  .x = as.list(names(borderline.performance$panel$performance)),
  .y = merge_plots,
  .f = function(.x, .y) {
    fn_save_auc(
      .filename = glue::glue("Borderline-{.x}-auc-merge.pdf"),
      .plot = .y
    )
  }
)

# Merge metrics -----------------------------------------------------------

merge_metrics <- fn_get_merge_metrics(.list = borderline.performance)
readr::write_tsv(x = merge_metrics, file = glue::glue("data/output/Borderline-metrics-merge.tsv"))
writexl::write_xlsx(x = merge_metrics, path = glue::glue("data/output/Borderline-metrics-merge.xlsx"))
# Save image --------------------------------------------------------------

save.image(file = "data/rda/09-performance-borderline.rda")
load(file = "data/rda/09-performance-borderline.rda")
