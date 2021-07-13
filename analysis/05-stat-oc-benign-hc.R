# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Jul 12 01:59:18 2021
# @DESCRIPTION: 05-stat-oc-benign-hc.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(DESeq2)
library(mlr)

# Load data ---------------------------------------------------------------

wuhan.se <- readr::read_rds(file = "data/rda/wuhan.se.rds.gz")
tom.se <- readr::read_rds(file = "data/rda/tom.se.rds.gz")
wuhan.tom.fs.fg.norm.rbe.se <- readr::read_rds(file = "data/rda/wuhan.tom.fs.fg.norm.rbe.se.rds.gz")

# Function ----------------------------------------------------------------


# Stat --------------------------------------------------------------------

wuhan.se@colData %>%
  as.data.frame() %>%
  dplyr::mutate(type = plyr::revalue(x = type, replace = c("benign" = "Benign", "malignant" = "Malignant", "normal" = "Health", "unkown" = "Malignant"))) %>%
  dplyr::select(barcode, oc = oc, type) %>%
  dplyr::mutate(oc = plyr::revalue(x = oc, replace = c("OC521" = "TC", "OC44" = "DC", "OC79" = "VC1", "OC172" = "VC2"))) ->
  wuhan.samples


tom.se@colData %>%
  as.data.frame() %>%
  dplyr::mutate(type = plyr::revalue(x = patientGroup, replace = c("HC" = "Health", "ovarianBenign" = "Benign", "ovarianBorderline" = "Malignant", "ovarianCancer" = "Malignant"))) %>%
  dplyr::mutate(oc = "VC3") %>%
  dplyr::select(barcode, oc, type) ->
  tom.samples

wuhan.tom.samples <- dplyr::bind_rows(wuhan.samples, tom.samples)


wuhan.tom.fs.fg.norm.rbe.se@colData %>%
  as.data.frame() %>%
  dplyr::inner_join(wuhan.tom.samples, by = "barcode") ->
  wuhan.tom.used.samples

wuhan.tom.used.samples %>%
  dplyr::group_by(oc, type) %>%
  tidyr::nest() %>%
  dplyr::mutate(age = purrr::map(.x = data, .f = function(.x) {
    .age <- .x$age
    .age %>% quantile() -> .q
    .all <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    .q <- .age[.age<=45] %>% quantile()
    .less45 <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    .q <- .age[.age>45] %>% quantile()
    .greater45 <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    tibble::tibble(
      name = c('Age', '<=45', '>45'),
      iqr = c(.all, .less45, .greater45)
    )
  })) %>%
  dplyr::ungroup() %>%
  dplyr::select(-3) %>%
  tidyr::unnest(age)  ->
  wuham.tom.used.samples.age

wuham.tom.used.samples.age %>%
  dplyr::filter(name == "Age") %>%
  dplyr::select(-name) %>%
  tidyr::spread(key = oc, value = iqr) ->
  wuham.tom.used.samples.age.spread
writexl::write_xlsx(x = wuham.tom.used.samples.age.spread, path = glue::glue("data/reviseoutput/age.xlsx"))


wuhan.tom.used.samples %>%
  dplyr::group_by(type) %>%
  tidyr::nest() %>%
  dplyr::mutate(age = purrr::map(.x = data, .f = function(.x) {
    .age <- .x$age
    .age %>% quantile() -> .q
    .all <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    .q <- .age[.age<=45] %>% quantile()
    .less45 <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    .q <- .age[.age>45] %>% quantile()
    .greater45 <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    tibble::tibble(
      name = c('Age', '<=45', '>45'),
      iqr = c(.all, .less45, .greater45)
    )
  })) %>%
  dplyr::ungroup() %>%
  dplyr::select(-2) %>%
  tidyr::unnest(age)  ->
  wuham.tom.used.samples.age.type

wuham.tom.used.samples.age.type %>%
  dplyr::filter(name == "Age") %>%
  dplyr::select(-name) %>%
  tidyr::spread(key = type, value = iqr) ->
  wuham.tom.used.samples.age.spread.type
writexl::write_xlsx(x = wuham.tom.used.samples.age.spread.type, path = glue::glue("data/reviseoutput/age.type.xlsx"))

# Save image --------------------------------------------------------------

save.image(file = "data/rda/05-stat-oc-benign-hc.rda")
