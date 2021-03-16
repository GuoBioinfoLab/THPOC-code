
# Library -----------------------------------------------------------------

library(magrittr)
library(mlr)

# Load data ---------------------------------------------------------------

el.task <- readr::read_rds(file = "data/rda/el.task.rds.gz")
epi.task <- readr::read_rds(file = "data/rda/epi.task.rds.gz")
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
  dplyr::select(barcode, stage, epi.non.epi) %>%
  dplyr::left_join()->
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

dplyr::bind_rows(
  test.el.ca125,
  test.epi.ca125
) %>%
  dplyr::distinct() ->
  test.el.epi.ca125

# merge -------------------------------------------------------------------

dplyr::bind_rows(train.el.epi.ca125, test.el.epi.ca125) %>%
  dplyr::filter(!is.na(CA125)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(group = ifelse(group == "B", "Non-OC", group)) %>%
  dplyr::mutate(group = factor(x = group, levels = c("Late-stage", "Early-stage", "Epithelial", "Non-epithelial", "Non-OC"))) ->
  clean.el.epi.ca125

clean.el.epi.ca125 %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  dplyr::mutate(label = glue::glue("{group}\n(n={n})")) ->
  clean.el.epi.ca125.label

# Plot --------------------------------------------------------------------
color_epi_lebn <- c('Late stage' = "#e31a1c", 'Early stage' = "#1f78b4", 'Epithelial' = "#399205", 'Non-epithelial' = '#e8a733', 'Non-cancer' = '#000000')
clean.el.epi.ca125 %>%
  ggplot(aes(x = group, y = CA125, color = group))+
  stat_boxplot(geom = 'errorbar', width = 0.2) +
  geom_boxplot(outlier.colour = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6, size = 1.3) +
  scale_x_discrete(name = 'Class', labels = clean.el.epi.ca125.label$label) +
  scale_color_manual(values = unname(color_epi_lebn)) +
  labs(x = 'Class', y = 'Normalized CA125 level') +
  theme(
    panel.background = element_rect(fill = NA, color = 'black'),
    axis.text = element_text(color = 'black'),
    legend.position = 'none'
  ) ->
  clean.el.epi.ca125_plot
ggsave(
  filename = 'CA125-plot.pdf',
  plot = clean.el.epi.ca125_plot,
  device = 'pdf',
  path = "data/output",
  width = 5, height = 3.5
)


# Save image --------------------------------------------------------------

save.image(file = "data/rda/02-ca125-boxplot.rda")
