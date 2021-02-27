
# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(ggridges)
library(mlr)
library(DESeq2)

# src ---------------------------------------------------------------------

source(file = "src/performance.R")

# path_tmp ----------------------------------------------------------------

path_tmp <- "/home/liucj/tmp/THPOC-cutoff"


# Function ----------------------------------------------------------------

fn_compare_metrics <- function(.x) {
  print(.x)
  .path_output <- file.path(path_tmp, .x , "output")

  tryCatch(
    expr = {
      purrr::map_df(
        .x = list("ca125", "panel", "panel_ca125"),
        .f = function(.x) {
          .filename <- glue::glue("bm.{.x}.metrics.tsv")
          readr::read_tsv(file = file.path(.path_output, .filename)) %>%
            tibble::add_column(predictor = .x, .before = 1)
        }
      )
    },
    error = function(e) NULL
  )
}

fn_perf <- function(.x, .y) {
  .task <- .y$task
  .sample <- .y$samples$Tom
  # .rms <- c("Vumc-HD-92-TR920", "TR4299-OVA-CATH", "Cath-Ova-CZE-001-TR2270", "Cath-Ova-CZE-032-TR2717", "TR3662-OVA-CATH", "TR3664-OVA-CATH", "TR4294-OVA-CATH", "TR4300-OVA-CATH", "TR4010-OVA-GDANSK", "TR4214-OVA-Supernat", "Vumc-HD-15-1-TR1154")
  # .rms <- c("TR4299-OVA-CATH","TR3664-OVA-CATH","TR4294-OVA-CATH","TR4300-OVA-CATH","TR4010-OVA-GDANSK")
  # .rms <- c("TR3664-OVA-CATH","TR4294-OVA-CATH")
  # .rms <- c("TR3664-OVA-CATH","TR4300-OVA-CATH")
  .rms <- c("Cath-Ova-CZE-024-TR2771","TR4320-OVA-GDANSK","TR4377-OVA-GDANSK","TR4379-OVA-GDANSK")
  .rms <- c("Cath-Ova-CZE-024-TR2771","TR4320-OVA-GDANSK","TR4377-OVA-GDANSK")
  .rms <- c("Cath-Ova-CZE-024-TR2771")
  # .rms <- c()

  .rm <- setdiff(names(.sample), .rms)

  .sample
  .pred <- predict(
    object = .x,
    task = .task,
    subset = .sample[.rm]
  )
  .pred
}

fn_compare_predict <- function(.x) {
  print(.x)
  .path_rda <- file.path(path_tmp, .x, "rda")

  .tasks <- c("panel" = list(readr::read_rds(file = file.path(.path_rda, "wuhan.tom.panel.task.rds.gz"))), readr::read_rds(file = file.path(.path_rda, "wuhan.tom.panel.ca125.task.rds.gz")))

  .models <- purrr::map(
    .x = names(.tasks),
    .f = function(.x) {
      .modelname <- glue::glue("{.x}.model.rds.gz")
      readr::read_rds(file = file.path(.path_rda, .modelname))
    }
  )
  names(.models) <- names(.tasks)

  purrr::map2(
    .x = .models,
    .y = .tasks,
    .f = fn_perf
  ) ->
    .bm.performance
  .bm.performance %>%
    tibble::enframe(name = "predictor", value = "pred")
}


fn_compare_samples <- function(.x) {

  .d <- .x$data %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(predictright = truth == response)

  .d %>%
    dplyr::group_by(truth, predictright) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(name = glue::glue("{truth}_{predictright}")) %>%
    dplyr::select(n, name) %>%
    tidyr::spread(key = name, value = n) ->
    .dtn

  .d %>%
    dplyr::filter(!predictright) %>%
    dplyr::pull(barcode) ->
    .es

  tibble::tibble(
    auc = performance(.x, measures = mlr::auc),
    acc = performance(.x, measures = mlr::acc),
    n = nrow(.x$data),
    predictright = sum(.d$predictright),
    predicterror = sum(!.d$predictright),
    error_samples = list(.es)
  ) %>%
    dplyr::bind_cols(.dtn)
}
# cutoff --------------------------------------------------------------------

list.files(path = path_tmp) %>%
  tibble::enframe(value = "cutoff") %>%
  dplyr::select(-name) %>%
  dplyr::mutate(
    metrics = purrr::map(
      .x = cutoff,
      .f = fn_compare_metrics
    )
  ) %>%
  dplyr::filter(purrr::map_lgl(.x = metrics, .f = Negate(is.null))) ->
  all_metrics

all_metrics %>%
  tidyr::unnest(metrics) %>%
  dplyr::filter(cohort == "Tom") %>%
  dplyr::mutate(
    auc = gsub(
      pattern = " \\(.*$",
      replacement = "",
      x = `AUC (95% CI)`)
  ) %>%
  dplyr::mutate(
    acc = gsub(
      pattern = " \\(.*$",
      replacement = "",
      x = `Accuracy (95% CI)`)
  ) %>%
  dplyr::select(cutoff, predictor, auc, acc) %>%
  dplyr::mutate(
    cutoff = as.numeric(cutoff),
    auc = as.numeric(auc),
    acc = as.numeric(acc)
  ) %>%
  tidyr::replace_na(replace = list(auc = 0, acc = 0)) ->
  all_metrics_tom
all_metrics_tom %>%
  ggplot(aes(x = cutoff, y = auc)) +
  geom_point() +
  facet_wrap(~ predictor) +
  geom_vline(xintercept = 0.4) +
  geom_hline(yintercept = 0.844)
all_metrics_tom %>%
  ggplot(aes(x = cutoff, y = acc)) +
  geom_point() +
  facet_wrap(~ predictor) +
  geom_vline(xintercept = 0.4)

# Predict right -----------------------------------------------------------


all_metrics %>%
  dplyr::mutate(
    predictright = purrr::map(
      .x = cutoff,
      .f = fn_compare_predict
    )
  ) ->
  all_predicts


all_predicts %>%
  dplyr::select(-metrics) %>%
  tidyr::unnest(predictright) %>%
  dplyr::mutate(
    result = purrr::map(.x = pred, .f = fn_compare_samples)
  ) %>%
  dplyr::select(-pred) %>%
  tidyr::unnest(cols = result) ->
  all_predicts_stat

all_predicts_stat %>%
  dplyr::mutate(cutoff = as.numeric(cutoff)) %>%
  dplyr::mutate(predictor = as.factor(predictor)) ->
  all_predicts_stat_d

all_predicts_stat_d %>%
  ggplot(aes(x = cutoff, y = auc, color = predictor))+
  geom_point() +
  geom_line() +
  scale_color_manual(values =  RColorBrewer::brewer.pal(n=3, name = 'Set1')[c(3,1,2)]) +
  scale_x_continuous(breaks = all_predicts_stat_d$cutoff) +
  geom_hline(yintercept = 0.8, linetype = 11) +
  geom_hline(yintercept = 0.85, color = "red", linetype = 11) +
  theme_bw() +
  labs(
    title = "Tom data AUC",
    x = "Mapping rate",
    y = "AUC"
  ) -> p_auc

all_predicts_stat_d %>%
  ggplot(aes(x = cutoff, y = acc, color = predictor))+
  geom_point() +
  geom_line() +
  scale_color_manual(values =  RColorBrewer::brewer.pal(n=3, name = 'Set1')[c(3,1,2)]) +
  scale_x_continuous(breaks = all_predicts_stat_d$cutoff) +
  geom_hline(yintercept = 0.8, linetype = 11) +
  geom_hline(yintercept = 0.85, color = "red", linetype = 11) +
  theme_bw() +
  labs(
    title = "Tom data ACC",
    x = "Mapping rate",
    y = "ACC"
  ) -> p_acc

merge_plots <- cowplot::plot_grid(
  plotlist = list(p_auc, p_acc),
  ncol = 1,
  align = 'v',
  rel_heights = c(1,1)
);merge_plots

ggsave(
  filename = "mergeplots1.pdf",
  path = path_tmp,
  plot = merge_plots,
  device = 'pdf',
  width = 8,
  height = 6
  )

all_predicts_stat %>%
  dplyr::mutate(cutoff = as.numeric(cutoff)) %>%
  dplyr::filter(predictor == "panel") %>%
  dplyr::filter(cutoff <= 0.33) %>%
  dplyr::pull(error_samples) %>%
  purrr::reduce(.f = intersect) -> .v1

all_predicts_stat %>%
  dplyr::mutate(cutoff = as.numeric(cutoff)) %>%
  dplyr::filter(predictor == "ca125") %>%
  dplyr::filter(cutoff <= 0.33) %>%
  dplyr::pull(error_samples) %>%
  purrr::reduce(.f = intersect) -> .v2

all_predicts_stat %>%
  dplyr::mutate(cutoff = as.numeric(cutoff)) %>%
  dplyr::filter(predictor == "panel_ca125") %>%
  dplyr::filter(cutoff <= 0.33) %>%
  dplyr::pull(error_samples) %>%
  purrr::reduce(.f = intersect) -> .v3

setdiff(intersect(.v1, .v3), .v2)

tom.se <- readr::read_rds(file = "/workspace/liucj/project/09-THPOC/rda/tom.se.rds.gz")
wuhan.se <- readr::read_rds(file = "/workspace/liucj/project/09-THPOC/rda/wuhan.se.rds.gz")


tom.se@colData %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  # dplyr::filter(barcode %in% .rms) %>%
  dplyr::mutate(ca125 = as.numeric(CA125parameterTOC)) %>%
  dplyr::select(barcode, class, mapping_rate, X__mapped_reads, ca125) %>%
  dplyr::filter(!is.na(ca125), ca125 > 0) %>%
  dplyr::mutate(ca125 = ifelse(ca125 > 10000, 10000, ca125)) %>%
  dplyr::mutate(cohort = "tom") %>%
  dplyr::mutate(ca125 = scale(log2(ca125))[, 1]) ->
  .tom_d

wuhan.se@colData %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  dplyr::mutate(class = as.factor(as.character(class))) %>%
  dplyr::mutate(cohort = oc) %>%
  dplyr::select(barcode, class, mapping_rate, X__mapped_reads, ca125 = CA125, cohort) %>%
  dplyr::group_by(cohort) %>%
  dplyr::mutate(ca125 = scale(log2(ca125))[, 1]) %>%
  dplyr::ungroup() ->
  .wuhan_d

dplyr::bind_rows(.tom_d, .wuhan_d) %>%
  ggplot(aes(x = cohort, y = ca125, color = class)) +
  geom_boxplot()
