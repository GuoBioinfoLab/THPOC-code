# Meta info -----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-03-1 10:00:58
# @DESCRIPTION: 03-validation-data-use.R

library(magrittr)
library(mlr)
library(DESeq2)

# Load data ---------------------------------------------------------------

el.task <- readr::read_rds(file = "data/rda/el.task.rds.gz")
epi.task <- readr::read_rds(file = "data/rda/epi.task.rds.gz")
endo.task <- readr::read_rds(file = "data/rda/endo.task.rds.gz")
wuhan.se <- readr::read_rds(file = "data/rda/wuhan.se.rds.gz")
tom.se <- readr::read_rds(file = "data/rda/tom.se.rds.gz")
wuhan.tom.se <- readr::read_rds(file = "data/rda/wuhan.tom.fs.fg.norm.rbe.se.rds.gz")

metadata <- wuhan.tom.se@colData %>%
  as.data.frame() %>%
  tibble::as_tibble()


metadata %>%
  dplyr::mutate(cohort = factor(cohort, levels = c("TC", "DC", "VC1", "VC2", "Tom"))) %>%
  dplyr::group_by(cohort, class) %>%
  dplyr::count() %>%
  dplyr::ungroup()

el.task %>%
  purrr::map(.f = "samples") %>%
  purrr::map(.f = function(.x) {
    .x %>%
      purrr::map(.f = function(.xx) {
      .xx %>%
        tibble::enframe(name = "barcode", value = "index") %>%
        dplyr::left_join(metadata, by = "barcode") %>%
        dplyr::select(barcode, class, cohort) %>%
          dplyr::group_by(cohort, class) %>%
          dplyr::count() %>%
          dplyr::ungroup() %>%
          dplyr::mutate(group = glue::glue("{cohort}-{class}")) %>%
          dplyr::select(group, n) %>%
          tidyr::spread(key = group, value = n)
    }) %>%
      tibble::enframe(name = "stage")
  }) %>%
  tibble::enframe() %>%
  tidyr::unnest(cols = value) %>%
  tidyr::unnest(cols = value)




epi.task %>%
  purrr::map(.f = "samples") %>%
  purrr::map(.f = function(.x) {
    .x %>%
      purrr::map(.f = function(.xx) {
        .xx %>%
          tibble::enframe(name = "barcode", value = "index") %>%
          dplyr::left_join(metadata, by = "barcode") %>%
          dplyr::select(barcode, class, cohort) %>%
          dplyr::group_by(cohort, class) %>%
          dplyr::count() %>%
          dplyr::ungroup() %>%
          dplyr::mutate(group = glue::glue("{cohort}-{class}")) %>%
          dplyr::select(group, n) %>%
          tidyr::spread(key = group, value = n)
      }) %>%
      tibble::enframe(name = "stage")
  }) %>%
  tibble::enframe() %>%
  tidyr::unnest(cols = value) %>%
  tidyr::unnest(cols = value)

endo.task %>%
  purrr::map(.f = "samples") %>%
  purrr::map(.f = function(.x) {
    .x %>%
      purrr::map(.f = function(.xx) {
        .xx %>%
          tibble::enframe(name = "barcode", value = "index") %>%
          dplyr::left_join(metadata, by = "barcode") %>%
          dplyr::select(barcode, class, cohort) %>%
          dplyr::group_by(cohort, class) %>%
          dplyr::count() %>%
          dplyr::ungroup() %>%
          dplyr::mutate(group = glue::glue("{cohort}-{class}")) %>%
          dplyr::select(group, n) %>%
          tidyr::spread(key = group, value = n)
      }) %>%
      tibble::enframe(name = "stage")
  }) %>%
  tibble::enframe() %>%
  tidyr::unnest(cols = value) %>%
  tidyr::unnest(cols = value)

# Save image --------------------------------------------------------------

save.image(file = "data/rda/03-validation-data-use.rda")
