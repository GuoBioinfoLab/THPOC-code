# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Jul 11 10:56:08 2021
# @DESCRIPTION: 07-performance-new-subtype.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(DESeq2)
library(mlr)

# Load data ---------------------------------------------------------------
wuhan.se <- readr::read_rds(file = "data/rda/wuhan.se.rds.gz")
tom.se <- readr::read_rds(file = "data/rda/tom.se.rds.gz")

bm.hc.task <- readr::read_rds(file = "data/rda/bm.hc.task.rds.gz")
model.list <- readr::read_rds(file = "data/rda/model.list.rds.gz")

# src ---------------------------------------------------------------------

source(file = "src/utils.R")
source(file = "src/performance.R")
source(file = "src/plots.R")

# Function ----------------------------------------------------------------

fn_get_el_task <- function(.task, .w, .t) {
  .w@colData %>%
    as.data.frame() %>%
    dplyr::filter(oc %in% c("OC79", "OC172")) %>%
    dplyr::mutate(stage = ifelse(type == "normal", "H", stage)) %>%
    dplyr::select(barcode, stage) ->
    .wd

  .t@colData %>%
    as.data.frame() %>%
    dplyr::select(barcode, stageFourGroups, Stage) %>%
    dplyr::mutate(stage = plyr::revalue(
      x = stageFourGroups,
      replace = c(
        "benign" = "B",
        "healthy control" = "H",
        "I" = "E",
        "II" = "E",
        "III" = "L",
        "IV" = "L"
    ))) %>%
    dplyr::mutate(stage = ifelse(Stage == "IIB", "L", stage)) %>%
    dplyr::select(barcode, stage) ->
    .td

  .wtd <- dplyr::bind_rows(.wd, .td)

  .early.hc <- .wtd %>% dplyr::filter(stage != "L")
  .late.hc <- .wtd %>% dplyr::filter(stage != "E")
  purrr::map(
    .x = .task,
    .f = function(.x) {
      .v <- purrr::reduce(.x$samples, .f = c)
      .s <- .v[!duplicated(names(.v))]

      .e <- c(na.omit(.s[c(.early.hc$barcode)]))
      .l <- c(na.omit(.s[c(.late.hc$barcode)]))

      .x$samples <- list(Early = .e, Late = .l)
      .x
  }) ->
    .el.hc.task

  .early.bam <- .early.hc %>% dplyr::filter(stage != "H")
  .late.bam <- .late.hc %>% dplyr::filter(stage != "H")
  purrr::map(
    .x = .task,
    .f = function(.x) {
      .v <- purrr::reduce(.x$samples, .f = c)
      .s <- .v[!duplicated(names(.v))]

      .e <- c(na.omit(.s[c(.early.bam$barcode)]))
      .l <- c(na.omit(.s[c(.late.bam$barcode)]))

      .x$samples <- list(Early = .e, Late = .l)
      .x
    }) ->
    .el.bam.task

  list(
    el.hc = .el.hc.task,
    el.bam = .el.bam.task
  )
}

fn_performance_ensemble <- function(.x, .y) {
  .perf <- fn_performance(.x, .y)
  .metrics <- fn_get_metrics(.perf)
  .plot <- fn_get_auc_plot(.perf, .metrics)

  list(
    performance = .perf,
    metrics = .metrics,
    plot = .plot
  )
}


# Early Late --------------------------------------------------------------

# task
el.hc.bam.task <- fn_get_el_task(.task = bm.hc.task, .w = wuhan.se, .t = tom.se)
readr::write_rds(
  x = el.hc.bam.task,
  file = "data/rda/el.hc.bam.task.rds.gz",
  compress = "gz"
)

# performance
purrr::map(
  .x = el.hc.bam.task,
  .f = function(.x) {
    purrr::map2(
      .x = model.list,
      .y = .x,
      .f = fn_performance_ensemble
    )
  }
) ->
  el.hc.bam.perfromance

readr::write_rds(
  x = el.hc.bam.perfromance,
  file = "data/rda/el.hc.bam.perfromance.rds.gz",
  compress = "gz"
)

# save metrics and plot
purrr::walk(
  .x = names(el.hc.bam.perfromance),
  .f = function(.x) {
    .prefix <- .x
    .y <- el.hc.bam.perfromance[[.x]]
    purrr::walk2(
      .x = list("panel", "panel_ca125", "ca125"),
      .y = .y,
      .f = function(.x, .y, .prefix) {
        readr::write_tsv(
          x = .y$metrics,
          file = glue::glue("data/reviseoutput/03-EL/{.prefix}.{.x}.metrics.tsv")
        )
        writexl::write_xlsx(
          x = .y$metrics,
          path = glue::glue("data/reviseoutput/03-EL/{.prefix}.{.x}.metrics.xlsx")
        )
        ggsave(
          filename = glue::glue("{.prefix}.{.x}.aucplot.pdf"),
          plot = .y$plot,
          device = "pdf",
          path = "data/reviseoutput/03-EL",
          width = 5.2,
          height = 4.5
        )
      },
      .prefix = .prefix
    )
  }
)

# merge plots
merge_plots <- purrr::map(
  .x = el.hc.bam.perfromance,
  .f = function(.x) {
    fn_get_merge_plots(
      .list = .x,
      .datasets = as.list(names(.x$panel$performance))
    )
  }
)

purrr::walk(
  .x = names(merge_plots),
  .f = function(.x) {
    .prefix <- toupper(.x)
    .y <- merge_plots[[.x]]
    purrr::walk2(
      .x = list("Early", "Late"),
      .y = .y,
      .f = function(.x, .y, .prefix) {
        ggsave(
          filename = glue::glue("{.prefix}.{.x}-auc-merge.pdf"),
          plot = .y,
          device = "pdf",
          path = "data/reviseoutput/03-EL",
          width = 5.2,
          height = 4.5
        )
      },
      .prefix = .prefix
    )
  }
)

# merge metrics
purrr::walk(
  .x = names(el.hc.bam.perfromance),
  .f = function(.x) {
    .prefix <- toupper(.x)
    .l <- el.hc.bam.perfromance[[.x]]
    .mm <- fn_get_merge_metrics(.list = .l)

    readr::write_tsv(x = .mm, file = glue::glue("data/reviseoutput/03-EL/{.prefix}-metrics-merge.tsv"))
    writexl::write_xlsx(x = .mm, path = glue::glue("data/reviseoutput/03-EL/{.prefix}-metrics-merge.xlsx"))
  }
)
