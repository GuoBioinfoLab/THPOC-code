
# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)
library(caret)
library(glmnet)
library(mRMRe)

# src ---------------------------------------------------------------------

source(file = "src/doparallel.R")
source(file = "src/utils.R")


# Load data ---------------------------------------------------------------

wuhan.tom.fs.fg.norm.rbe.se <- readr::read_rds(file = "data/rda/wuhan.tom.fs.fg.norm.rbe.se.rds.gz")

# Function ----------------------------------------------------------------

fn_select_features <- function(.wt) {
  .se <- .wt[, .wt$cohort == "TC"]
  .feats1 <- fn_sd_median(.se = .se)

  .df <- fn_se2df(.se = .se[.feats1, ])

  fn_parallel_start(n_cores = 50)
  set.seed(1234)
  .glm <- caret::train(
    class ~ .,
    data = .df,
    method = "glmnet",
    family = "binomial",
    trainControl = caret::trainControl("cv", number = 10),
    tuneGrid =  expand.grid(
      alpha = 1,
      lambda = 10^seq(-3, 3,length = 100)
    )
  )

  fn_parallel_stop()

  coef(.glm$finalModel, .glm$bestTune$lambda) %>%
    as.matrix() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'ensid') %>%
    dplyr::filter(`1` != 0) %>%
    dplyr::filter(grepl(pattern = "ENSG", x = ensid)) %>%
    dplyr::pull(ensid)
}

# Select features ---------------------------------------------------------


panel <- fn_select_features(.wt = wuhan.tom.fs.fg.norm.rbe.se)
readr::write_rds(x = panel, file = "data/rda/panel.rds.gz", compress = "gz")



# Save .image -------------------------------------------------------------
save.image(file = "data/rda/03-feature-selection.rda")

