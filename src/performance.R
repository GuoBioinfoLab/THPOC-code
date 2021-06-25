# Meta info -----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-03-1 10:00:58
# @DESCRIPTION: performance.R


fn_se2task_panel <- function(.se, .id = "task") {
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


fn_se2task_panel_ca125 <- function(.se) {
  .d <- cbind(
    t(assay(.se)),
    data.frame(
      CA125 = .se$CA125,
      class = .se$class
    )
  )
  .task_panel_ca125 <- mlr::makeClassifTask(
    id = "Panel-CA125-task",
    data = .d,
    target = "class",
    positive = "M"
  )
  .task_ca125 <- mlr::makeClassifTask(
    id = "CA125-task",
    data = .d[, c("class", "CA125")],
    target = "class",
    positive = "M"
  )
  list(
    task_panel_ca125 = .task_panel_ca125,
    task_ca125 = .task_ca125
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
    ParamHelpers::makeNumericParam("C", lower = -10, upper = 10, trafo = function(x) 10^x),
    ParamHelpers::makeNumericParam("sigma", lower = -10, upper = 10, trafo = function(x) 10^x)
  )
  .tune_algorithm <- mlr::makeTuneControlRandom(
    same.resampling.instance = TRUE, maxit = 100L
  )
  mlr::configureMlr(
    show.info = FALSE,
    on.learner.error = "warn",
    on.measure.not.applicable = "warn"
  )
  parallelMap::parallelStart(mode = "multicore", cpus = 100)
  .tune_result <- mlr::tuneParams(
    learner = .learner, task = .tsk,
    resampling = .cv10i,
    measures = list(mlr::acc, mlr::mmce, mlr::auc),
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
  .task_name <- gsub(pattern = "-task", replacement = "", x = .task_id)
  .hped <- mlr::generateHyperParsEffectData(tune.result = .tune_result, trafo = TRUE)

  .plot_hpe <- mlr::plotHyperParsEffect(
    hyperpars.effect.data = .hped,
    x = "iteration", y = "acc.test.mean",
    plot.type = "line"
  ) +
    theme_bw() +
    theme() +
    guides(shape = guide_legend(title = "Learner Status"), colour = guide_legend(title = "Learner Status")) +
    labs(
      x = "Iteration",
      y = "Accuracy test mean",
      title = glue::glue("Random search iteration for training {.task_name}")
    )
  ggsave(
    filename = glue::glue("01-Tune-parameter-{.task_name}.pdf"),
    plot = .plot_hpe,
    device = "pdf",
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
  names(.ci) <- c("lower", "auc", "upper")
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
    kappa = .kappa_f1["kappa"],
    f1 = .kappa_f1["f1"]
  )
}


fn_predidct_performance_metrics <- function(.x, .model, .task, .samples) {
  .pred <- predict(
    object = .model,
    task = .task,
    subset = .samples[[.x]]
  )
  .pred$data %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tibble::add_column(cohort = .x, .before = 1) ->
    .pred_data

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
      .funs = ~ purrr::map(.x = ., .f = function(.x) {
        glue::glue('{sprintf("%.3f", .x[1])} ({sprintf("%.3f", .x[2])} - {sprintf("%.3f", .x[3])})')
      })
    ) %>%
    dplyr::mutate_if(
      .predicate = rlang::is_double,
      .funs = ~ purrr::map(.x = ., .f = function(.x) {
        sprintf("%.3f", .x)
      })
    ) %>%
    dplyr::mutate_if(
      .predicate = rlang::is_list,
      .funs = unlist
    ) %>%
    tibble::add_column(cohort = .x, .before = 1)

  names(.pred_metrics) <- c("cohort", "AUC (95% CI)", "Accuracy (95% CI)", "SN (95% CI)", "SP (95% CI)", "PPV (95% CI)", "NPV (95% CI)", "Kappa", "F1")
  list(
    perf = .pred_perf,
    metrics = .pred_metrics,
    pred_data = .pred_data
  )
}


fn_performance <- function(.model, .list) {
  .task <- .list$task
  .samples <- .list$samples

  purrr::map(
    names(.samples),
    fn_predidct_performance_metrics,
    .model = .model,
    .task = .task,
    .samples = .samples
  ) ->
  .perf

  names(.perf) <- names(.samples)
  .perf
}

fn_get_metrics <- function(.perf) {
  .perf %>%
    purrr::map("metrics") %>%
    purrr::reduce(.f = dplyr::bind_rows)
}


fn_get_auc_plot <- function(.perf, .metrics) {
  .metrics %>%
    dplyr::mutate(auc = gsub(pattern = " ", replacement = "", x = `AUC (95% CI)`)) %>%
    dplyr::select(cohort, auc) %>%
    dplyr::mutate(label = glue::glue("{cohort} {auc}")) ->
  .labels

  .d <- .perf %>%
    purrr::map("perf") %>%
    purrr::reduce(.f = dplyr::bind_rows) %>%
    dplyr::mutate(cohort = factor(x = cohort, levels = .labels$cohort))

  fn_plot_auc(.d, .labels)
}

fn_get_merge_plots <- function(.list, .datasets) {
  .panel <- .list$panel$performance
  .panel_ca125 <- .list$panel_ca125$performance
  .ca125 <- .list$ca125$performance

  purrr::map2_df(
    .x = list("CA125", "THPOC", "THPOC + CA125"),
    .y = list(.ca125, .panel, .panel_ca125),
    .f = function(.x, .y) {
      .y %>%
        purrr::map("perf") %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(type = .x) %>%
        # dplyr::filter(cohort != "Tom") %>%
        dplyr::select(-threshold)
    }
  ) ->
  .merge_data

  purrr::map(
    .x = .datasets,
    .f = fn_plot_merge_auc,
    .d = .merge_data
  ) ->
  .plots

  names(.plots) <- .datasets
  .plots
}

fn_delong_test <- function(.x) {
  .x %>% dplyr::filter(name == "panel") -> .x.panel
  .x %>% dplyr::filter(name == "ca125") -> .x.ca125
  .x %>% dplyr::filter(name == "panel_ca125") -> .x.panel_ca125

  dplyr::bind_rows(
    pROC::roc.test(
      pROC::roc(.x.panel$truth, .x.panel$prob.M),
      pROC::roc(.x.ca125$truth, .x.ca125$prob.M),
      method = 'delong'
    ) %>%
      broom::tidy() %>%
      dplyr::mutate(`P value` = format(p.value, digits = 3,scientific = TRUE)) %>%
      dplyr::select(`P value`) %>%
      tibble::add_column(name = "THPOC", .before = 1),
    tibble::tibble(name = "CA125", `P value` = "-"),
    pROC::roc.test(
      pROC::roc(.x.panel_ca125$truth, .x.panel_ca125$prob.M),
      pROC::roc(.x.ca125$truth, .x.ca125$prob.M),
      method = 'delong'
    ) %>%
      broom::tidy() %>%
      dplyr::mutate(`P value` = format(p.value, digits = 3,scientific = TRUE)) %>%
      dplyr::select(`P value`) %>%
      tibble::add_column(name = "THPOC + CA125", .before = 1)
  ) %>%
    dplyr::mutate(`P value` = gsub(pattern = "e", replacement = "x10^", x = `P value`))
}

fn_get_merge_metrics <- function(.list) {
  .list %>%
    purrr::map(.f = function(.x) {
      .x$performance %>%
        purrr::map("pred_data")
    }) %>%
    tibble::enframe() %>%
    tidyr::unnest(value) %>%
    tidyr::unnest(value) %>%
    dplyr::group_by(cohort) %>%
    tidyr::nest() %>%
    dplyr::ungroup() ->
    .for_delong_test

  .for_delong_test %>%
    dplyr::mutate(p.value = purrr::map(.x = data, .f = fn_delong_test)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(p.value) %>%
    dplyr::mutate(name = factor(name, levels = c("THPOC", "CA125", "THPOC + CA125")))->
    .pvalue

  .list %>%
    purrr::map("metrics") %>%
    tibble::enframe() %>%
    tidyr::unnest(cols = value) %>%
    dplyr::mutate(name = plyr::revalue(x = name, replace = c("panel" = "THPOC", "ca125" = "CA125", "panel_ca125" = "THPOC + CA125"))) %>%
    dplyr::mutate(name = factor(name, levels = c("THPOC", "CA125", "THPOC + CA125"))) %>%
    dplyr::arrange(cohort, name) %>%
    dplyr::select(2, 1, 3:10) %>%
    dplyr::inner_join(.pvalue, by = c("cohort", "name"))
}
