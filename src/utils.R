# Meta info -----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-03-1 10:00:58
# @DESCRIPTION: utils.R

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
  assay(.xse[, .xse$class == .x]) %>%
    apply(MARGIN = 1, FUN = ineq::ineq, type = "var") %>%
    tibble::enframe()
}


fn_remove_unwanted_variables <- function(.se, .vars = c("cohort", "age", "lib.size", "library"), .th_pval = 0.05, .th_var_strong = 0.1) {
  .vars <- if (!"class" %in% .vars) c(.vars, "class") else .vars
  .mod <- model.matrix(~class, data = .se@colData)
  .mod0 <- model.matrix(~1, data = .se@colData)
  .svobj <- sva(dat = assay(.se), mod = .mod, mod0 = .mod0, n.sv = 100)

  .matrix_w_corr_vars <- foreach(i = 1:ncol(.svobj$sv), .combine = rbind, .packages = c("magrittr")) %dopar% {
    .vars %>%
      purrr::map_dbl(
        .f = function(.x) {
          if (.x %in% c("class", "cohort", "library")) {
            aov(formula = .svobj$sv[, i] ~ .se@colData[, .x]) %>%
              broom::tidy() %>%
              dplyr::slice(1) %>%
              dplyr::pull(p.value)
          } else {
            cor.test(
              x = .svobj$sv[, i],
              y = .se@colData[, .x],
              method = "pearson"
            ) %>%
              broom::tidy() %>%
              dplyr::slice(1) %>%
              dplyr::pull(p.value)
          }
        }
      ) -> .vars_pval

    names(.vars_pval) <- .vars

    if (.vars_pval["class"] > .th_pval) {
      .strong_var <- names(sort(x = .vars_pval))[1]
    } else {
      .strong_var <- "class"
    }
    if (.vars_pval[.strong_var] < .th_var_strong) {
      .verdict <- names(.vars_pval[.strong_var])
    } else {
      .verdict <- NA
    }

    c(.vars_pval, verdict = .verdict)
  }

  rownames(.matrix_w_corr_vars) <- 1:ncol(.svobj$sv)

  .vars %>%
    purrr::map(.f = function(.x) {
      unname(which(.matrix_w_corr_vars[, "verdict"] == .x))
    }) -> confounding
  names(confounding) <- .vars

  confounding$na <- unname(which(is.na(.matrix_w_corr_vars[, "verdict"])))

  confounding$confounding <- setdiff(1:ncol(.svobj$sv), c(confounding$class, confounding$na))

  readr::write_rds(x = list(confounding = confounding,svobj = .svobj,mod = .mod, dbdat = assay(.se)), file = "data/rda/confounding-svobj.rds")

  .data_rm_be <- removeBatchEffect(assay(.se), design = .mod, covariates = .svobj$sv[, confounding$confounding])

  SummarizedExperiment(assays = .data_rm_be, colData = .se@colData[colnames(.data_rm_be), ])
}

fn_sd_median <- function(.se) {
  .sd <- assay(.se) %>%
    apply(MARGIN = 1, FUN = sd)
  .md <- assay(.se) %>%
    apply(MARGIN = 1, FUN = median)

  names(.md[.md > 0])
}

fn_se2df <- function(.se) {
  .df <- as.data.frame(t(assay(.se)))
  .target <- .se$class
  .df$class <- .target
  .df
}


human_read <- function(.x){
  .sign = ifelse(.x < 0 , TRUE, FALSE)
  .x <- abs(.x)

  if (.x >= 0.1) {
    .x %>% signif(digits = 2) %>% toString() -> .xx
  } else if (.x < 0.1 && .x >= 0.001 ) {
    .x %>% signif(digits = 2) %>% toString() -> .xx
  } else if (.x < 0.001 && .x > 0) {
    .x %>% format(digits = 3, scientific = TRUE) -> .xx
  } else {
    .xx <- '0'
  }

  ifelse(.sign, paste0('-',.xx), .xx)
}

human_read_latex_pval <- function(.x, .s = NA) {

  if (is.na(.s)) {
    if (grepl(pattern = "e", x = .x)) {
      sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
      latex2exp::TeX(glue::glue("$<<.xx[1]>> \\times 10^{<<.xx[2]>>}$", .open = "<<", .close = ">>"))
    } else {
      latex2exp::TeX(glue::glue("$<<.x>>$", .open = "<<", .close = ">>"))
    }
  } else {
    if (grepl(pattern = "e", x = .x)) {
      sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
      latex2exp::TeX(glue::glue("<<.s>>, $<<.xx[1]>> \\times 10^{<<.xx[2]>>}$", .open = "<<", .close = ">>"))
    } else {
      latex2exp::TeX(glue::glue("<<.s>>, $<<.x>>$", .open = "<<", .close = ">>"))
    }
  }

}
