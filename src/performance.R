

fn_bm_se2task <- function(.se, .id = "task") {
  .d <- cbind(
    t(assay(.se)),
    data.frame(
      class = .se$class
      )
    )

  mlr::makeClassifTask(
    id = .id,
    data = .d,
    target = "class",
    positive = "M"
  )
}


fn_task_ind <- function(.x, .se) {
  .ind <- which(.se$cohort == .x)
  setNames(.ind, nm = .se@colData[.ind, "barcode"])
}
