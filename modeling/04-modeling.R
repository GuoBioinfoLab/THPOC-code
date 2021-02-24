
# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)
library(mlr)
library(doParallel)

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

  set.seed(123)
  .cv10d <- mlr::makeResampleDesc(method = "CV", iters = 10, stratify = TRUE)
  .cv10i <- mlr::makeResampleInstance(desc = .cv10d, task = .task_for_tunes)
  .learner <- mlr::makeLearner(
    cl = "classif.ksvm",
    id = "svm-learner",
    predict.type = "prob"
  )
  .hyperparameters <- ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("C", lower = -10, upper = 10, trafo = function(x) 10 ^ x),
    ParamHelpers::makeNumericParam("sigma", lower = -10, upper = 10, trafo = function(x) 10 ^ x)
  )
  .tune_algorithm <-  mlr::makeTuneControlRandom(
    same.resampling.instance = TRUE, maxit = 500L
  )
  mlr::configureMlr(
    show.info = FALSE,
    on.learner.error = 'warn',
    on.measure.not.applicable = 'warn'
  )

  parallelMap::parallelStart(mode = 'multicore', cpus = 100)
  .tune_result <-  mlr::tuneParams(
    learner = .learner, task = .task_for_tunes,
    resampling = .cv10i,
    measures = list(mlr::auc, mlr::mmce, mlr::acc),
    par.set = .hyperparameters,
    control = .tune_algorithm
  )
  parallelMap::parallelStop()

  .hped <- mlr::generateHyperParsEffectData(.tune_result, trafo = TRUE)
  .plot_hpe <- mlr::plotHyperParsEffect(.hped, x = "iteration", y = "auc.test.mean", plot.type = "line")
  
  .learner_tuned <- mlr::setHyperPars(learner = .learner, par.vals = .tune_result$x)
  .model_tuned <- mlr::train(learner = .learner_tuned, task = .task, subset = .samples$train)

}

# Prepare task ------------------------------------------------------------
wuhan.tom.task <- fn_task(.wt = wuhan.tom.fs.fg.norm.rbe.se, .panel = panel)

