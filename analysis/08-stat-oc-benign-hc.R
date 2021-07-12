# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Jul 12 01:59:18 2021
# @DESCRIPTION: 08-stat-oc-benign-hc.R

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
    quantile(.x$age)
  }))
