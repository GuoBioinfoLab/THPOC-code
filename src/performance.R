

fn_bm_se2task <- function(.se, .id = "task") {
  .d <- cbind(
    t(assay(.se)),
    data.frame(
      class = .se$class
      )
    )

  mlr::makeClassifTask(
    id = .id,
    data = .d,
    target = "class",
    positive = "M"
  )
}


fn_task_ind <- function(.x, .se) {
  .ind <- which(.se$cohort == .x)
  setNames(.ind, nm = .se@colData[.ind, "barcode"])
}

fn_tune_model <- function(.tsk) {
  set.seed(123)
  .cv10d <- mlr::makeResampleDesc(method = "CV", iters = 5, stratify = TRUE)
  .cv10i <- mlr::makeResampleInstance(desc = .cv10d, task = .tsk)
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
    same.resampling.instance = TRUE, maxit = 10L
  )
  mlr::configureMlr(
    show.info = FALSE,
    on.learner.error = 'warn',
    on.measure.not.applicable = 'warn'
  )

  parallelMap::parallelStart(mode = 'multicore', cpus = 100)
  .tune_result <-  mlr::tuneParams(
    learner = .learner, task = .tsk,
    resampling = .cv10i,
    measures = list(mlr::auc, mlr::mmce, mlr::acc),
    par.set = .hyperparameters,
    control = .tune_algorithm
  )
  parallelMap::parallelStop()


  list(
    tune_result = .tune_result,
    learner = mlr::setHyperPars(learner = .learner, par.vals = .tune_result$x)
  )
}

fn_plot_tune_path <- function(.tune_result, .task_id) {
  .task_name <- gsub(pattern = '-task', replacement = '', x = .task_id)
  .hped <-  mlr::generateHyperParsEffectData(tune.result = .tune_result, trafo = TRUE)

  .plot_hpe <- mlr::plotHyperParsEffect(
    hyperpars.effect.data = .hped,
    x = "iteration", y = "auc.test.mean",
    plot.type = "line"
  ) +
    theme_bw() +
    theme() +
    guides(shape = guide_legend(title = 'Learner Status'), colour = guide_legend(title = 'Learner Status')) +
    labs(
      x = 'Iteration',
      y = 'Accuracy test mean',
      title = glue::glue('Random search iteration for training {.task_name}')
    )
  ggsave(
    filename = glue::glue('01-Tune-parameter-{.task_name}.pdf'),
    plot = .plot_hpe,
    device = 'pdf',
    path = "data/output",
    width = 8,
    height = 4
  )

  .plot_hpe
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
  tibble::tibble(
    auc = list(unlist(c(.auc$auc[2], .auc$auc[1], .auc$auc[3]))),
    acc = list(unlist(.rocm_epi$rval$diag.acc)),
    tpr = list(unlist(.rocm_epi$rval$se)),
    tnr = list(unlist(.rocm_epi$rval$sp)),
    ppv = list(unlist(.rocm_epi$rval$ppv)),
    npv = list(unlist(.rocm_epi$rval$npv)),
    kappa = .kappa_f1['kappa'],
    f1 = .kappa_f1['f1']
  )

}

fn_predidct_performance_metrics <- function(.x, .model, .task, .samples) {
  .pred <- predict(
    object = .model,
    task = .task,
    subset = .samples[[.x]]
    )

  .pred_perf <- mlr::generateThreshVsPerfData(
    obj = .pred,
    measures = list(mlr::fpr, mlr::tpr, mlr::auc),
    gridsize = length(.samples[[.x]])
  ) %>%
    .$data %>%
    tibble::add_column(cohort = .x, .before = 1)


  .pred_metrics <- fn_roc_95ci(.pred) %>%
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
    ) %>%
    tibble::add_column(cohort = .x, .before = 1)

  names(.pred_metrics)  <- c("cohort", "AUC (95% CI)", 'Accuracy (95% CI)', 'SN (95% CI)', 'SP (95% CI)', 'PPV (95% CI)', 'NPV (95% CI)', 'Kappa', 'F1')
  list(
    perf = .pred_perf,
    metrics = .pred_metrics
  )
}


fn_plot_auc <- function(.d, .labels) {
  .d %>%
    ggplot(aes(x = fpr, y = tpr, color = cohort)) +
    geom_path(size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = 11) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1), expand = c(0, 0)) +
    scale_color_manual(
      name = 'AUC',
      labels = .labels$label,
      values = RColorBrewer::brewer.pal(n=5, name = 'Set1')
    ) +
    theme(
      panel.background = element_rect(fill = NA),
      panel.grid = element_blank(),
      panel.border = element_blank(),

      axis.line.x.bottom = element_line(color = 'black'),
      axis.line.y.left = element_line(color = 'black'),
      axis.ticks.length = unit(x = 0.2, units = 'cm'),
      axis.text = element_text(color = 'black', size = 14),
      axis.title = element_text(color = 'black', size = 16),

      legend.position = c(0.68, 0.2),
      legend.background = element_rect(fill = NA),
      legend.key = element_rect(fill = NA),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.key.width = unit(1.5, units = 'cm'),
      legend.spacing = unit(c(0,0,0,0), units = 'cm'),
      legend.title.align = 1,

      plot.margin = unit(c(1,1,0.5,0.5), units = 'cm'),
      plot.title = element_text(hjust = 0.5, size = 18)
    ) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity"
    )
}
