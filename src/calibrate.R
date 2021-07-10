# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Jul  9 08:18:54 2021
# @DESCRIPTION: calibrate.R

#' Function to predict mlr task by mlr model
#'
#' @param .x sample subset id
#' @param .model mlr model
#' @param .task mlr task
#' @param .samples mlr task samples
#' @return mlr pred object
fn_mlr_predict <- function(.x, .model, .task, .samples) {
  predict(
    object = .model,
    task = .task,
    subset = .samples[[.x]]
  )
}

#' Input function for purrr::map
#' to return pred
#'
#' @param .x mlr model
#' @param .y list of mlr task and sample id
#' @return list pred
fn_get_mlr_pred <- function(.x, .y) {
  .model <- .x
  .task <- .y$task
  .samples <- .y$samples

  purrr::map(
    .x = names(.samples),
    .f = fn_mlr_predict,
    .model = .model,
    .task = .task,
    .samples = .samples
  ) ->
    .perf
  names(.perf) <- names(.samples)
  .perf
}


#' Generate before calibrate and after calibrate data
#'
#' @param .x mlr prediction object
#' @return list of before and after calibrate
fn_get_calibrate <- function(.x) {
  .bc_pred <- .x
  .bc_auc_brier <- performance(
    pred = .bc_pred,
    measures = list(auc, brier)
  )

  .ac_pred <- .x
  .ac_pred$data <- fn_platt_scaling(.ac_pred$data)
  .ac_auc_brier <- performance(
    pred = .ac_pred,
    measures = list(auc, brier)
  )

  list(
    before_calib = list(pred = .bc_pred, auc_brier = .bc_auc_brier),
    after_calib = list(pred = .ac_pred, auc_brier = .ac_auc_brier)
  )
}


#' Platt Scaling - Linear logistic calibration
#'
#' @param .x prediction data before calibration
#' @return prediction data after calibration
fn_platt_scaling <- function(.d) {
  .d %>%
    dplyr::mutate(y = ifelse(truth == "M", 1, 0)) ->
    .dd

  .fit <- glm(y~prob.M, data = .dd, family = binomial())
  .p <- predict(.fit, .dd[3], type = "response")

  .dd %>%
    dplyr::mutate(prob.M = .p) %>%
    dplyr::mutate(prob.B = 1 - prob.M) %>%
    dplyr::mutate(response = ifelse(prob.M >=0.5, "M", "B")) %>%
    dplyr::select(-y)
}


#' Plot calibration curve
#'
#' @param calibdata the calibration data prediction and auc brier score
#' @return calibration curve ggplot figure
fn_plot_calibration_curve <- function(.x) {
  .x <- bm.calibrate$VC1$after_calib
  .pred <- .x$pred
  .pred_data <- .pred$data
  .auc_brier <- .x$auc_brier

  .auc_brier %>%
    tibble::enframe() %>%
    dplyr::mutate(name = c("AUC", "Brier")) %>%
    dplyr::mutate(label = glue::glue("{name}:{round(value, digits = 3)}")) %>%
    dplyr::mutate(x = c(0.1, 0.1), y = c(0.9, 0.85)) ->
    .labels


  .pred_data %>%
    dplyr::mutate(pred = prob.M) %>%
    dplyr::mutate(prob = ifelse(truth == "M", 1, 0)) %>%
    dplyr::mutate(bin = dplyr::ntile(pred, 10)) %>%
    dplyr::group_by(bin) %>%
    dplyr::mutate(
      n = dplyr::n(),
      bin_pred = mean(pred),
      bin_prob = mean(prob),
      se = sqrt((bin_prob * (1 - bin_prob)) / n),
      ul = bin_prob + 1.96 * se,
      ll = bin_prob - 1.96 * se
    ) %>%
    dplyr::ungroup() ->
    .for_plot

  .for_plot %>%
    ggplot(aes(x = bin_pred, y = bin_prob)) +
    geom_point(shape = 2) +
    geom_abline(linetype = 17) +
    geom_smooth(aes(x = pred, y = prob), color = "red", se = FALSE, method = "loess") +
    geom_text(data = .labels, aes(x = x, y = y, label = label)) +
    labs(
      x = "Predicted probability",
      y = "Actual probability"
    )

  .for_plot %>%
    ggplot(aes(x = pred)) +
    geom_histogram(fill = "black", bins = 10) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))


}




