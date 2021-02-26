fn_auc_theme <- function(.legend.position = c(0.8, 0.2)) {

  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),

    axis.line.x.bottom = element_line(color = "black"),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.length = unit(x = 0.2, units = 'cm'),
    axis.text = element_text(color = 'black', size = 14),
    axis.title = element_text(color = 'black', size = 18),

    legend.position = .legend.position,
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.key.width = unit(1.8, units = 'cm'),
    legend.spacing = unit(c(0,0,0,0), units = 'cm'),
    legend.title.align = 0,

    plot.margin = unit(c(1,1,0.5,0.5), units = 'cm'),
    plot.title = element_text(hjust = 0.5, size = 18)
  )
}


fn_plot_auc <- function(.d, .labels) {
  .lgl <-  any(grepl(pattern = "-", .labels$label))
  .legend.position <- if(.lgl) c(0.5, 0.5) else c(0.8, 0.2)

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
    labs(
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    fn_auc_theme(.legend.position)

}


fn_plot_merge_auc <- function(.x, .d) {
  .d %>%
    dplyr::filter(cohort == .x) %>%
    dplyr::mutate(type = as.factor(type))->
    .dd

  .legend_title <-  glue::glue("AUC for {.x}")
  .legend_text <- .dd %>%
    dplyr::select(auc, type) %>%
    dplyr::distinct() %>%
    dplyr::mutate(auc_label = round(auc, digits = 3)) %>%
    dplyr::pull(auc_label)
  .legend_color <- c("#006400", "#B22222", "#00008B")

  .dd %>%
    ggplot(aes(x = fpr, y = tpr, color = type)) +
    geom_abline(intercept = 0, slope = 1, linetype = 11) +
    geom_path(size = 0.8) +
    scale_x_continuous(
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
      limits = c(0, 1), expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_color_manual(
      name = .legend_title,
      labels = .legend_text,
      values = .legend_color
    ) +
    guides(
      color = guide_legend(
        reverse = TRUE
      )
    ) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    fn_auc_theme()
}


fn_save_auc <- function(.filename, .plot) {
  ggsave(
    filename = .filename,
    plot = .plot,
    device = "pdf",
    path = "data/output",
    width = 5.2,
    height = 4.5
  )
}



fn_get_tom_plot <- function(.perf) {
  .legend_title <- glue::glue("AUC for Tom")
  .legend_text <- .perf %>%
    dplyr::mutate(auc_label = round(auc, digits = 3)) %>%
    dplyr::distinct(auc_label) %>%
    dplyr::pull(auc_label)
  .legend_color <- "#B22222"

  .perf %>%
    ggplot(aes(x = fpr, y = tpr, color = cohort)) +
    geom_path(size = 0.8) +
    scale_x_continuous(
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
      limits = c(0, 1), expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_color_manual(
      name = .legend_title,
      labels = .legend_text,
      values = .legend_color
    ) +
    guides(
      color = guide_legend(
        reverse = TRUE
      )
    ) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    fn_auc_theme()
}
