# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Jul 12 06:44:06 2021
# @DESCRIPTION: merge-grade-histotype.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(DESeq2)


# Load data ---------------------------------------------------------------

# tom.se
tom.se <- readr::read_rds(file = "data/rda/tom.se.rds.gz")
tom.se@colData %>%
  as.data.frame() %>%
  dplyr::select(barcode, name) ->
  tom.df

wuhan.grade <- purrr::map(.x = c(1,2,3,4), .f = function(.x) {
  readxl::read_xlsx(path = "data/raw/grade-serous/grade-new.xlsx", sheet = .x)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::select(barcode, oc, grade)

wuhan.histotype <- readxl::read_xlsx(path = "data/raw/grade-serous/serous.xlsx") %>%
  dplyr::select(barcode, histotype)

wuhan.grade.histotype <- wuhan.grade %>% dplyr::inner_join(wuhan.histotype, by = "barcode")

wuhan.grade.histotype %>% dplyr::group_by(grade, histotype) %>% dplyr::count()

tom.grade.histotype <- readxl::read_xlsx(path = "data/raw/grade-serous/Tom-grade.xlsx") %>%
  dplyr::select(name, grade = Grade, histotype = morphologicalDiagnosis) %>%
  dplyr::left_join(tom.df, by = "name") %>%
  dplyr::mutate(oc = "Tom") %>%
  dplyr::select(barcode, oc, grade, histotype)

dplyr::bind_rows(
  wuhan.grade.histotype,
  tom.grade.histotype
) %>%
  readr::write_rds(file = "data/rda/wuhan.tom.grade.histotype.rds.gz", compress = "gz")

# Save image --------------------------------------------------------------

save.image(file = "data/rda/merge-grade-histotype.rda")
