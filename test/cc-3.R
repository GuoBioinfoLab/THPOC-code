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

cc.fs.fg.norm.rbe.se <- readr::read_rds(file = "data/rda/cc.fs.fg.norm.rbe.se.rds")

# Load model --------------------------------------------------------------


panel.model <- readr::read_rds(file = "data/rda/panel.model.rds.gz")
panel <- readr::read_rds(file = "data/rda/panel.rds.gz")

# Function ----------------------------------------------------------------
fn_se2df <- function(.se, .panel) {
  .sub_se <- .se[.panel,]
  .d <- cbind(
    t(assay(.sub_se)),
    data.frame(class = factor(x = as.character(colData(.sub_se)[, 'class']), levels = c('M', 'B')))
  )
  .d
}
fn_auc_ci <- function(.p) {
  .truth <- mlr::getPredictionTruth(pred = .p)
  .prob <- mlr::getPredictionProbabilities(pred = .p)
  .roc <- pROC::roc(response = .truth, predictor = .prob, plot = FALSE, ci = TRUE)
  .ci <- as.numeric(.roc$ci)
  names(.ci) <- c('lower', 'auc', 'upper')
  .data <- data.frame(
    fpr = 1 - .roc$specificities,
    tpr = .roc$sensitivities
  )
  list(
    data = .data,
    auc = .ci
  )
}
fn_perf <- function(.pred) {
  .pred_perf <- fn_auc_ci(.p = .pred)
  .pred_perf
}

fn_roc_95ci <- function(.p) {
  # acc - Accuracy
  # tpr - True positive rate (Sensitivity, Recall)
  # tnr - True negative rate (Specificity)
  # ppv - Positive predictive value (Precision)
  # npv - Negative predictive value
  # kappa - Cohen's kappa
  # f1 - F1 measure
  # mlr::performance(pred = .pred$panel, measures = list(mlr::acc, mlr::tpr, mlr::tnr, mlr::ppv, mlr::npv, mlr::kappa, mlr::f1))
  .auc <- fn_auc_ci(.p = .p)
  .rocm <- mlr::calculateROCMeasures(.p)
  .kappa_f1 <- mlr::performance(pred = .p, measures = list(mlr::kappa, mlr::f1))
  .rocm_epi <- epiR::epi.tests(dat = t(.rocm$confusion.matrix))
  .pred_metrics <- tibble::tibble(
    auc = list(unlist(c(.auc$auc[2], .auc$auc[1], .auc$auc[3]))),
    acc = list(unlist(.rocm_epi$rval$diag.acc)),
    tpr = list(unlist(.rocm_epi$rval$se)),
    tnr = list(unlist(.rocm_epi$rval$sp)),
    ppv = list(unlist(.rocm_epi$rval$ppv)),
    npv = list(unlist(.rocm_epi$rval$npv)),
    kappa = .kappa_f1['kappa'],
    f1 = .kappa_f1['f1']
  )
  .table_metrics <- .pred_metrics %>%
    dplyr::mutate_if(
      .predicate = rlang::is_list,
      .funs = ~purrr::map(.x = ., .f = function(.x) {glue::glue('{sprintf("%.3f", .x[1])} ({sprintf("%.3f", .x[2])} - {sprintf("%.3f", .x[3])})')})) %>%
    dplyr::mutate_if(
      .predicate = rlang::is_double,
      .funs = ~purrr::map(.x = ., .f = function(.x) {sprintf("%.3f", .x)})
    ) %>%
    dplyr::mutate_if(
      .predicate = rlang::is_list,
      .funs = unlist
    )
  colnames(.table_metrics) <- c('AUC (95% CI)', 'Accuracy (95% CI)', 'SN (95% CI)', 'SP (95% CI)', 'PPV (95% CI)', 'NPV (95% CI)', 'Kappa', 'F1')
  .table_metrics
}
fn_performance <- function(.se, .model, .panel) {
  # .se <- cc.fs.fg.norm.rbe.se
  # .model <- panel.model
  # .panel <- panel

  .se.df <- fn_se2df(.se = .se, .panel = .panel)
  .pred <- predict(object = .model, newdata = .se.df)
  .pred_data <- .pred$data %>%
    dplyr::mutate(predictright = ifelse(truth == response, TRUE, FALSE)) %>%
    dplyr::group_by(truth, predictright) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(`truth-predictright` = glue::glue("{truth} {predictright}")) %>%
    dplyr::select(-c(truth, predictright)) %>%
    tidyr::spread(key = `truth-predictright`, value = n) %>%
    dplyr::mutate(total_sample = sum(dplyr::c_across()))
  .pred_metrics <- fn_roc_95ci(.pred)
  .metrics <- .pred_data %>%
    dplyr::bind_cols(.pred_metrics)

  list(
    metrics = .pred_metrics,
    pred_data = .pred$data
  )
  .metrics
}

fn_perf_disease <- function(.x) {
  cc.fs.fg.norm.rbe.se@colData %>%
    as.data.frame() %>%
    dplyr::filter(Disease %in% c(.x, "Healthy")) %>%
    dplyr::pull(barcode) ->
    .barcode
  .se <- cc.fs.fg.norm.rbe.se[, .barcode]
  fn_performance(.se = .se, .model = panel.model, .panel = panel)
}
# Predict -----------------------------------------------------------------
cc.fs.fg.norm.rbe.se@colData %>%
  as.data.frame() %>%
  dplyr::filter(Disease != "Healthy") %>%
  dplyr::group_by(Disease) %>%
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-n) %>%
  dplyr::filter(n >=10) ->
  disease_count

disease_count %>%
  dplyr::mutate(perf = purrr::map(.x = Disease, .f = fn_perf_disease)) %>%
  tidyr::unnest(perf) ->
  perf_disease

perf_disease %>%
  writexl::write_xlsx(path = "data/output/cc-diease-perf.xlsx")


# Save --------------------------------------------------------------------

save.image(file = "data/rda/cc-3.rda")
load(file = "data/rda/cc-3.rda")
