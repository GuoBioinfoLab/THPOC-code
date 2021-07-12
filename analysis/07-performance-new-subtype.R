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

Endometriosis_samples <- readxl::read_xlsx(path = 'data/raw/wuhan/endoinfo-V5-metadata-1.xlsx', sheet = 3)
Borderline_samples <- readxl::read_xlsx(path = "data/raw/borderline.xlsx") %>%
  dplyr::filter(!is.na(grade)) %>%
  dplyr::filter(barcode != "name") %>%
  dplyr::filter(barcode != "barcode") %>%
  dplyr::filter(oc != "OC44")

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

fn_get_epi_task <- function(.task, .w, .t) {
  .w@colData %>%
    as.data.frame() %>%
    dplyr::filter(oc %in% c("OC79", "OC172")) %>%
    dplyr::filter(!is.na(epi.non.epi)) %>%
    dplyr::mutate(epi.non.epi = ifelse(type == "normal", "H", epi.non.epi)) %>%
    dplyr::mutate(epi.non.epi = ifelse(epi.non.epi == "Benign", "B", epi.non.epi)) %>%
    dplyr::select(barcode, epi.non.epi) ->
    .wd

  .t@colData %>%
    as.data.frame() %>%
    dplyr::mutate(epi.non.epi = ifelse(class == "M", "epithelial", "H")) %>%
    dplyr::mutate(epi.non.epi = ifelse(stageFourGroups == "benign", "B", epi.non.epi)) %>%
    dplyr::select(barcode, epi.non.epi) ->
    .td

  .wtd <- dplyr::bind_rows(.wd, .td)

  .epi.hc <- .wtd %>% dplyr::filter(epi.non.epi != "non-epithelial")
  .nepi.hc <- .wtd %>% dplyr::filter(epi.non.epi != "epithelial")
  purrr::map(
    .x = .task,
    .f = function(.x) {
      .v <- purrr::reduce(.x$samples, .f = c)
      .s <- .v[!duplicated(names(.v))]

      .epi <- c(na.omit(.s[c(.epi.hc$barcode)]))
      .nepi <- c(na.omit(.s[c(.nepi.hc$barcode)]))

      .x$samples <- list("Epi" = .epi, "Non-epi" = .nepi)
      .x
    }) ->
    .epi.hc.task

  .epi.bam <- .epi.hc %>% dplyr::filter(epi.non.epi != "H")
  .nepi.bam <- .nepi.hc %>% dplyr::filter(epi.non.epi != "H")
  purrr::map(
    .x = .task,
    .f = function(.x) {
      .v <- purrr::reduce(.x$samples, .f = c)
      .s <- .v[!duplicated(names(.v))]

      .epi <- c(na.omit(.s[c(.epi.bam$barcode)]))
      .nepi <- c(na.omit(.s[c(.nepi.bam$barcode)]))

      .x$samples <- list("Epi" = .epi, "Non-epi" = .nepi)
      .x
    }) ->
    .epi.bam.task

  list(
    epi.hc = .epi.hc.task,
    epi.bam = .epi.bam.task
  )
}

fn_get_endo_task <- function(.task, .t, .endo) {

  .t@colData %>%
    as.data.frame() %>%
    dplyr::filter(patientGroup %in% c("ovarianCancer", "ovarianBorderline")) %>%
    dplyr::select(barcode, class = patientGroup) ->
    .te

  .endo %>%
    dplyr::filter(oc %in% c("OC79", "OC172")) %>%
    dplyr::select(barcode, class) ->
    .we

  .wtd <- dplyr::bind_rows(.we, .te)
  .wtd %>% dplyr::group_by(class) %>% dplyr::count()

  purrr::map(
    .x = .task,
    .f = function(.x) {
      .v <- purrr::reduce(.x$samples, .f = c)
      .s <- .v[!duplicated(names(.v))]

      .e <- c(na.omit(.s[c(.wtd$barcode)]))

      .x$samples <- list(Endometriosis = .e)
      .x
    }) ->
    .endo.task

  list(
    endo.hc = .endo.task,
    endo.bam = .endo.task
  )
}

fn_get_bl_task <- function(.task, .w, .t, .bl) {
  .w@colData %>%
    as.data.frame() %>%
    dplyr::filter(oc %in% c("OC79", "OC172")) %>%
    dplyr::mutate(class = ifelse(type == "normal", "H", stage)) %>%
    dplyr::select(barcode, class) %>%
    dplyr::filter(class %in% c("B", "H")) ->
    .wd

  .t@colData %>%
    as.data.frame() %>%
    dplyr::filter(patientGroup %in% c("HC", "ovarianBenign")) %>%
    dplyr::select(barcode, class = patientGroup) %>%
    dplyr::mutate(class = ifelse(class == "HC", "H", "B")) ->
    .td

  dplyr::bind_rows(
    .wd, .td,
    .bl %>% dplyr::select(barcode, class) %>% dplyr::mutate(class = "M")
  ) ->
    .wtd

  .bl.hc <- .wtd
  purrr::map(
    .x = .task,
    .f = function(.x) {
      .v <- purrr::reduce(.x$samples, .f = c)
      .s <- .v[!duplicated(names(.v))]

      .b <- c(na.omit(.s[c(.bl.hc$barcode)]))

      .x$samples <- list(Borderline = .b)
      .x
    }) ->
    .bl.hc.task

  .bl.bam <- .bl.hc %>% dplyr::filter(class != "H")
  purrr::map(
    .x = .task,
    .f = function(.x) {
      .v <- purrr::reduce(.x$samples, .f = c)
      .s <- .v[!duplicated(names(.v))]

      .b <- c(na.omit(.s[c(.bl.bam$barcode)]))

      .x$samples <- list(Borderline = .b)
      .x
    }) ->
    .bl.bam.task

  list(
    bl.hc = .bl.hc.task,
    bl.bam = .bl.bam.task
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

# Early Late --------------------------------------------------------------

# task
el.hc.bam.task <- fn_get_el_task(.task = bm.hc.task, .w = wuhan.se, .t = tom.se)
readr::write_rds(
  x = el.hc.bam.task,
  file = "data/rda/el.hc.bam.task.rds.gz",
  compress = "gz"
)

# predict subtype
fn_predict_subtype(
  .task = el.hc.bam.task,
  .name = "el.hc.bam",
  .out = "03-EL"
)

# Epi ---------------------------------------------------------------------
# task
epi.hc.bam.task <- fn_get_epi_task(.task = bm.hc.task, .w = wuhan.se, .t = tom.se)
readr::write_rds(
  x = epi.hc.bam.task,
  file = "data/rda/epi.hc.bam.task.rds.gz",
  compress = "gz"
)

# predict subtype
fn_predict_subtype(
  .task = epi.hc.bam.task,
  .name = "epi.hc.bam",
  .out = "04-EPI"
)


# Endo --------------------------------------------------------------------
endo.hc.bam.task <- fn_get_endo_task(.task = bm.hc.task, .t = tom.se, .endo = Endometriosis_samples)
readr::write_rds(
  x = endo.hc.bam.task,
  file = "data/rda/endo.hc.bam.task.rds.gz",
  compress = "gz"
)

# predict subtype
fn_predict_subtype(
  .task = endo.hc.bam.task,
  .name = "endo.hc.bam",
  .out = "05-ENDO"
)

# Borderline --------------------------------------------------------------
bl.hc.bam.task <- fn_get_bl_task(.task = bm.hc.task, .w = wuhan.se, .t = tom.se, .bl = Borderline_samples)
readr::write_rds(
  x = bl.hc.bam.task,
  file = "data/rda/bl.hc.bam.task.rds.gz",
  compress = "gz"
)

# predict subtype
fn_predict_subtype(
  .task = bl.hc.bam.task,
  .name = "bl.hc.bam",
  .out = "06-BL"
)


# Grade -------------------------------------------------------------------


# Serous ------------------------------------------------------------------


# Save image --------------------------------------------------------------

save.image(file = "data/rda/07-performance-new-subtype.rda")
load(file = "data/rda/07-performance-new-subtype.rda")
