
fn_load_htseq_count <- function(.txt, .path) {
  .name <- gsub(pattern = "\\..*", replacement = "", x = .txt)

  readr::read_tsv(file = file.path(.path, .txt), col_names = c("Ensembl Gene ID", .name))
}


fn_filter_annotation <- function(.count) {
  .count %>%
    dplyr::filter(grepl(pattern = "__", x = `Ensembl Gene ID`)) ->
    .a

  .count %>%
    dplyr::filter(!grepl(pattern = "__", x = `Ensembl Gene ID`)) ->
    .b

  .b %>%
    dplyr::select(-`Ensembl Gene ID`) %>%
    dplyr::summarise_all(.funs = sum) %>%
    tidyr::gather(key = "barcode", value = "__mapped_reads") ->
    .bb

  .a %>%
    tidyr::gather(key = "barcode", value = "count", -`Ensembl Gene ID`) %>%
    tidyr::spread(key = `Ensembl Gene ID`, value = count) %>%
    dplyr::left_join(.bb, by = "barcode") %>%
    dplyr::mutate(`__unmapped_reads` = `__alignment_not_unique` + `__ambiguous` + `__no_feature` + `__not_aligned` + `__too_low_aQual`) %>%
    dplyr::mutate(`__reads_total` = `__mapped_reads` + `__unmapped_reads`) %>%
    dplyr::mutate(mapping_rate = `__mapped_reads` / `__reads_total`) ->
    .aa

  list(
    mapping_rate = .aa,
    clean_count = .b
  )
}


fn_filter_by_gene_count <- function(.se, minCounts = 10, fSample = 0.9, thres = 3) {

  apply(X = assay(.se), MARGIN = 1, FUN = function(x) {
    sum(x >= minCounts) / length(x) > fSample
  }) -> .a

  rownames(.se)[.a]
}


fn_filter_by_hyper <- function(.x, .xse) {
  assay(.xse[, .xse$class == .x]) %>% apply(MARGIN = 1, FUN = ineq::ineq, type = "var") %>%
    tibble::enframe()
}
