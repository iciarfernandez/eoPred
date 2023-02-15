#' @title predictPreeclampsia
#'
#' @description Uses 45 CpGs to predict early preeclampsia (PE delivered before or at 34 weeks of gestation)
#' on placental DNA methylation microarray data.
#'
#' @details Assigns the class labels "early-PE" or "normotensive" to each sample
#' and returns a class probability.
#'
#' It is recommended that users apply beta-mixture quantile normalization (BMIQ) to their data
#' prior to prediction. This was the normalization method used on the training data.
#'
#' @param betas matrix or array of methylation values on the beta scale (0, 1),
#' where the variables are arranged in rows, and samples in columns.
#'
#' @return produces a list with components detailed in the `mixOmics::predict`
#'
#' @examples
#'
#' To predict early preeclampsia on 450k/850k samples
#'
#' Load data
#' data(testBetas)
#' predictPreeclampsia(testBetas, dist = "max.dist")
#'
#' @importFrom mixOmics dplyr
#' @importFrom tibble magrittr
#'
#' @export predictPreeclampsia
#'

predictPreeclampsia <- function(betas, ...){

  # read in data
  mod <- eoPred:::mod

  trainCpGs <- colnames(mod$X)
  peCpGs <- mixOmics::selectVar(mod)$name

  # check that there are no NAs in the predictors (or if there are, how many)
  pp <- intersect(colnames(betas), peCpGs)

  if(length(pp) < length(peCpGs)){
    warning(paste(
      "Only", length(pp), "out of 45 predictors present."
    ))
  } else {
    message(paste(length(pp), "of 45 predictors present."))
    message("BMIQ normalization is recommended for best results.")
  }

  # set up data for prediction

  # if input data is missing any of the cpgs present in the training data, this function
  # adds the ones that are missing as NAs
  # necessary for `mixOmics::predict` to work

  outersect = function(x, y) {
    sort(c(x[!x%in%y],
           y[!y%in%x]))
  }

  if(inherits(betas, 'matrix')){
  } else if (inherits(betas, 'array')) {
  } else {

    # throw an error
    print(paste0("Input data must be a matrix or an array"))
  }

  subset <- betas[,colnames(betas) %in% trainCpGs]
  cpgs <- outersect(colnames(subset), trainCpGs)
  cpgs_df <- data.frame(1:length(cpgs))
  row.names(cpgs_df) <- cpgs
  xx <- rownames(subset)
  cpgs_df[xx] <- NA
  cpgs_df[,1] <- NULL
  subset <- cbind(subset, t(cpgs_df))

  # order
  subset <- subset[drop=FALSE,, trainCpGs]

  if(all(colnames(subset) == trainCpGs) == FALSE){
    stop()
  } else

    # predict
    out <- mixOmics:::predict.mixo_spls(mod, subset)

  # get class probabilities
  CP <- out$predict[,,1]
  CP <- t(apply(as.matrix(CP), 1, function(data) exp(data)/sum(exp(data))))
  CP <- as.data.frame(CP) %>% tibble::rownames_to_column("Sample_ID")
  CP$Pred_Class <- CP$comp1
  CP <- CP %>%
    dplyr::mutate(Pred_Class = dplyr::case_when(EOPE > `Non-PE Preterm` ~ "Early-PE",
                                                EOPE < `Non-PE Preterm` ~ "Normotensive")) %>%
    dplyr::rename(`Early-PE` = EOPE,
                  Normotensive = `Non-PE Preterm`)

  return(CP)
}
