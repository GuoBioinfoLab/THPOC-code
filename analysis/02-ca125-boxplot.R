
# Library -----------------------------------------------------------------

library(magrittr)
library(mlr)

# src ---------------------------------------------------------------------

source(file = "src/utils.R")

# Load data ---------------------------------------------------------------

el.task <- readr::read_rds(file = "data/rda/el.task.rds.gz")
epi.task <- readr::read_rds(file = "data/rda/epi.task.rds.gz")
endo.task <- readr::read_rds(file = "data/rda/endo.task.rds.gz")
borderline.task <- readr::read_rds(file = "data/rda/borderline.task.rds.gz")

wuhan.se <- readr::read_rds(file = "data/rda/wuhan.se.rds.gz")

panel_ca125.model <- readr::read_rds(file = "data/rda/panel_ca125.model.rds.gz")
panel <- readr::read_rds(file = "data/rda/panel.rds.gz")


# CA125 data --------------------------------------------------------------
ca125_data <- mlr::getTaskData(task = el.task$panel_ca125$task) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "barcode") %>%
  dplyr::select(barcode, CA125, class)


# Training data -----------------------------------------------------------


wuhan.se@colData %>%
  as.data.frame() %>%
  dplyr::filter(oc == "OC521") %>%
  dplyr::select(barcode, stage, epi.non.epi) ->
  training_data

training_data %>%
  dplyr::left_join(ca125_data, by = "barcode") %>%
  dplyr::mutate(group = plyr::revalue(x = stage, replace = c(
    "B" = "B",
    "E" = "Early-stage",
    "L" = "Late-stage"
  ))) %>%
  dplyr::select(barcode, CA125, group) ->
  train.el.ca125


training_data %>%
  dplyr::left_join(ca125_data, by = "barcode") %>%
  dplyr::mutate(group = plyr::revalue(x = epi.non.epi, replace = c(
    "Benign" = "B",
    "epithelial" = "Epithelial",
    "non-epithelial" = "Non-epithelial"
  ))) %>%
  dplyr::select(barcode, CA125, group) ->
  train.epi.ca125

dplyr::bind_rows(train.el.ca125, train.epi.ca125) %>%
  dplyr::distinct() ->
  train.el.epi.ca125

# Testing data ------------------------------------------------------------


el.task$panel_ca125$samples %>%
  purrr::map(.f = tibble::enframe, name = "barcode", value = "index") %>%
  tibble::enframe(name = "type") %>%
  tidyr::unnest(cols = value) %>%
  dplyr::left_join(ca125_data, by = "barcode") %>%
  dplyr::mutate(group = ifelse(class == "B", "B", type)) %>%
  dplyr::mutate(group = plyr::revalue(x = group, replace = c(
    "B" = "B",
    "Early" = "Early-stage",
    "Late" = "Late-stage"
  ))) %>%
  dplyr::select(barcode, CA125, group)->
  test.el.ca125

epi.task$panel_ca125$samples %>%
  purrr::map(.f = tibble::enframe, name = "barcode", value = "index") %>%
  tibble::enframe(name = "type") %>%
  tidyr::unnest(cols = value) %>%
  dplyr::left_join(ca125_data, by = "barcode") %>%
  dplyr::mutate(group = ifelse(class == "B", "B", type)) %>%
  dplyr::mutate(group = plyr::revalue(x = group, replace = c(
    "B" = "B",
    "Epi" = "Epithelial",
    "Non-epi" = "Non-epithelial"
  ))) %>%
  dplyr::select(barcode, CA125, group) ->
  test.epi.ca125

borderline.task$panel_ca125$samples %>%
  purrr::map(.f = tibble::enframe, name = "barcode", value = "index") %>%
  tibble::enframe(name = "type") %>%
  tidyr::unnest(cols = value) %>%
  dplyr::left_join(ca125_data, by = "barcode") %>%
  dplyr::mutate(group = ifelse(class == "B", "B", type)) %>%
  dplyr::select(barcode, CA125, group) ->
  test.borderline.ca125


endo.task$panel_ca125$samples %>%
  purrr::map(.f = tibble::enframe, name = "barcode", value = "index") %>%
  tibble::enframe(name = "type") %>%
  tidyr::unnest(cols = value) %>%
  dplyr::left_join(ca125_data, by = "barcode") %>%
  dplyr::mutate(group = ifelse(class == "B", "B", type)) %>%
  dplyr::select(barcode, CA125, group) ->
  test.endo.ca125

dplyr::bind_rows(
  test.el.ca125,
  test.epi.ca125,
  test.borderline.ca125,
  test.endo.ca125
) %>%
  dplyr::distinct() ->
  test.el.epi.borderline.endo.ca125

# merge -------------------------------------------------------------------

dplyr::bind_rows(train.el.epi.ca125, test.el.epi.borderline.endo.ca125) %>%
  dplyr::filter(!is.na(CA125)) %>%
  dplyr::mutate(CA125 = ifelse(CA125 < -3, -3, CA125)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(group = ifelse(group == "B", "Non-OC", group)) %>%
  dplyr::mutate(group = factor(x = group, levels = c("Late-stage", "Early-stage", "Epithelial", "Non-epithelial","Endometriosis", "Borderline", "Non-OC"))) %>%
  dplyr::filter(group != "Endometriosis") ->
  clean.el.epi.borderline.endo.ca125

clean.el.epi.borderline.endo.ca125 %>%
  dplyr::filter(group == "Non-OC") ->
  non_oc

clean.el.epi.borderline.endo.ca125 %>%
  dplyr::filter(group != "Non-OC") %>%
  dplyr::group_by(group) %>%
  tidyr::nest() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pval = purrr::map2_dbl(.x = data, .y = group, .f = function(.x, .y) {
    .x %>%
      dplyr::mutate(group = .y) %>%
      dplyr::bind_rows(non_oc) %>%
      t.test(CA125~group, data = .) %>%
      broom::tidy() %>%
      dplyr::pull(p.value)
  })) %>%
  dplyr::select(-data) ->
  pval

clean.el.epi.borderline.endo.ca125 %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  dplyr::mutate(label = glue::glue("{group}\n(n={n})")) %>%
  dplyr::left_join(pval, by = "group") %>%
  dplyr::mutate(pval = ifelse(is.na(pval), 1, pval)) %>%
  dplyr::mutate(pval_label = purrr::map(.x = pval, .f = function(.x) {
    human_read_latex_pval(human_read(.x))
  })) ->
  clean.el.epi.borderline.endo.ca125.label

# Plot --------------------------------------------------------------------

# color_epi_lebn <- c('Late stage' = "#e31a1c", 'Early stage' = "#1f78b4", 'Epithelial' = "#399205", 'Non-epithelial' = '#e8a733', "Endometriosis" = "#00FFFF", "Borderline" = "#8A2BE2", 'Non-cancer' = '#000000')
color_epi_lebn <- c('Late stage' = "#e31a1c", 'Early stage' = "#1f78b4", 'Epithelial' = "#399205", 'Non-epithelial' = '#e8a733', "Endometriosis" = "#00FFFF", 'Non-cancer' = '#000000')
clean.el.epi.borderline.endo.ca125 %>%
  ggplot(aes(x = group, y = CA125, color = group))+
  stat_boxplot(geom = 'errorbar', width = 0.2) +
  geom_boxplot(outlier.colour = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6, size = 1.3) +
  scale_x_discrete(name = 'Class', labels = clean.el.epi.borderline.endo.ca125.label$label) +
  scale_color_manual(values = unname(color_epi_lebn)) +
  labs(x = 'Class', y = 'Normalized CA125 level') +
  theme(
    panel.background = element_rect(fill = NA, color = 'black'),
    axis.text = element_text(color = 'black'),
    legend.position = 'none'
  ) +
  annotate(
    geom = "text",
    x = 1:5,
    y = 3,
    label = simplify2array(clean.el.epi.borderline.endo.ca125.label$pval_label[1:5])
  ) ->
  clean.el.epi.borderline.endo.ca125_plot;clean.el.epi.borderline.endo.ca125_plot


ggsave(
  filename = 'CA125-plot.pdf',
  plot = clean.el.epi.borderline.endo.ca125_plot,
  device = 'pdf',
  path = "data/output",
  width = 8, height = 4
)


# Save image --------------------------------------------------------------

save.image(file = "data/rda/02-ca125-boxplot.rda")
load(file = "data/rda/02-ca125-boxplot.rda")
