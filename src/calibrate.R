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
