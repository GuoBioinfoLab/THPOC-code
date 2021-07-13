# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Jul 13 08:25:39 2021
# @DESCRIPTION: 09-ca125-500.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(DESeq2)

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

fn_get_ca125_500_task <- function(.task, .w, .t) {
  .w@colData %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(ca125 = CA125) %>%
    dplyr::filter(!(ca125 > 500 & class == "B")) %>%
    dplyr::mutate(type = plyr::revalue(x = type, replace = c("benign" = "B", "malignant" = "M", "normal" = "H", "unkown" = "M"))) %>%
    dplyr::select(barcode, type, ca125) ->
    .wd

  .t@colData %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(ca125 = as.numeric(CA125parameterTOC)) %>%
    dplyr::mutate(type = plyr::revalue(x = patientGroup, replace = c("HC" = "H", "ovarianBenign" = "B", "ovarianBorderline" = "M", "ovarianCancer" = "M"))) %>%
    dplyr::filter(!(ca125 > 500 & class %in% c("H", "B"))) %>%
    dplyr::select(barcode, type, ca125) ->
    .td

  .wtd.hc <- dplyr::bind_rows(.wd, .td)
  .wtd.bam <- .wtd.hc %>% dplyr::filter(type != "H")

  purrr::map(
    .x = .task,
    .f = function(.x) {
      .samples <- purrr::map(.x = .x$samples, .f = function(.x){c(na.omit(.x[.wtd.hc$barcode]))})
      .x$samples <- .samples
      .x
    }
  ) -> .bm.hc

  purrr::map(
    .x = .task,
    .f = function(.x) {
      .samples <- purrr::map(.x = .x$samples, .f = function(.x){c(na.omit(.x[.wtd.bam$barcode]))})
      .x$samples <- .samples
      .x
    }
  ) -> .bm.bam

  list(
    bm.hc = .bm.hc,
    bm.bam = .bm.bam
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

fn_predict_subtype <- function(.task, .name, .out) {

  # performance
  purrr::map(
    .x = .task,
    .f = function(.x) {
      purrr::map2(
        .x = model.list,
        .y = .x,
        .f = fn_performance_ensemble
      )
    }
  ) ->
    .perfromance

  readr::write_rds(
    x = .perfromance,
    file = glue::glue("data/rda/{.name}.perfromance.rds.gz"),
    compress = "gz"
  )

  # save metrics and plot
  purrr::walk(
    .x = names(.perfromance),
    .f = function(.x) {
      .prefix <- .x
      .y <- .perfromance[[.x]]
      purrr::walk2(
        .x = list("panel", "panel_ca125", "ca125"),
        .y = .y,
        .f = function(.x, .y, .prefix) {
          readr::write_tsv(
            x = .y$metrics,
            file = glue::glue("data/reviseoutput/{.out}/{.prefix}.{.x}.metrics.tsv")
          )
          writexl::write_xlsx(
            x = .y$metrics,
            path = glue::glue("data/reviseoutput/{.out}/{.prefix}.{.x}.metrics.xlsx")
          )
          ggsave(
            filename = glue::glue("{.prefix}.{.x}.aucplot.pdf"),
            plot = .y$plot,
            device = "pdf",
            path = glue::glue("data/reviseoutput/{.out}"),
            width = 5.2,
            height = 4.5
          )
        },
        .prefix = .prefix
      )
    }
  )

  # merge plots
  .mps <- purrr::map(
    .x = .perfromance,
    .f = function(.x) {
      fn_get_merge_plots(
        .list = .x,
        .datasets = as.list(names(.x$panel$performance))
      )
    }
  )

  purrr::walk(
    .x = names(.mps),
    .f = function(.x) {
      .prefix <- toupper(.x)
      .y <- .mps[[.x]]
      purrr::walk2(
        .x =  as.list(names(.y)),
        .y = .y,
        .f = function(.x, .y, .prefix) {
          ggsave(
            filename = glue::glue("{.prefix}.{.x}-auc-merge.pdf"),
            plot = .y,
            device = "pdf",
            path = glue::glue("data/reviseoutput/{.out}"),
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
    .x = names(.perfromance),
    .f = function(.x) {
      .prefix <- toupper(.x)
      .l <- .perfromance[[.x]]
      .mm <- fn_get_merge_metrics(.list = .l)

      readr::write_tsv(x = .mm, file = glue::glue("data/reviseoutput/{.out}/{.prefix}-metrics-merge.tsv"))
      writexl::write_xlsx(x = .mm, path = glue::glue("data/reviseoutput/{.out}/{.prefix}-metrics-merge.xlsx"))
    }
  )
}


# Performance -------------------------------------------------------------

ca125.500.task <- fn_get_ca125_500_task(.task = bm.hc.task, .w = wuhan.se, .t = tom.se)

readr::write_rds(
  x = ca125.500.task,
  file = "data/rda/ca125.500.task.rds.gz",
  compress = "gz"
)

# predict subtype
fn_predict_subtype(
  .task = ca125.500.task,
  .name = "ca125.500",
  .out = "09-CA125-500"
)

# Save image --------------------------------------------------------------

save.image(file = "data/rda/09-ca125-500.rda")
