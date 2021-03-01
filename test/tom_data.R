tom.se <- readr::read_rds(file = "data/rda/tom.se.rds.gz")
wuhan.tom.fs.fg.norm.rbe.se <- readr::read_rds(file = "data/rda/wuhan.tom.fs.fg.norm.rbe.se.rds.gz")
panel <- readr::read_rds(file = "data/rda/panel.rds.gz")

tom.se@colData %>%
  as.data.frame()


wuhan.tom.fs.fg.norm.rbe.se@colData %>%
  as.data.frame() %>%
  dplyr::filter(cohort == "Tom") ->
  .d

tom.se@colData %>%
  as.data.frame() %>%
  dplyr::filter(barcode %in% .d$barcode) %>%
  dplyr::select(!dplyr::starts_with("X__")) ->
  tom_clinical
tom_clinical %>%
  writexl::write_xlsx("data/raw/tom_meta.xlsx")

table(tom_clinical$hospital.of.blood.collection, tom_clinical$type)

tom_clinical %>%
  dplyr::mutate(ca125 = as.numeric(CA125parameterTOC)) %>%
  dplyr::filter(ca125 > 0) ->
  tom_clinical_ca125


table(tom_clinical_ca125$hospital.of.blood.collection, tom_clinical_ca125$class)
