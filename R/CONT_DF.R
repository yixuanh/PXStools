#' Toy data with continuous outcome
#'
#' Artificially generated toy data set of 5000 individuals with a numeric outcome
#'  and a set of categorical and continuous exposure variables
#'
#' @docType data
#'
#' @usage data(CONT_DF)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#' #' \describe{
#'  \item{IID}{ID's of individuals}
#'  \item{SEX}{sex of individuals (character)}
#'  \item{AGE}{age (numeric)}
#'  \item{COV_Q_OTHER}{quantitative covariates (numeric)}
#'  \item{COV_C_OTHER}{categorical covariate (categorical)}
#'  \item{PC_6-14}{genetic principal components (randomly generated)}
#'  \item{VAR1-18}{numeric exposure variables (numeric)}
#'  \item{VAR19-33}{categorical exposure variables (categorical)}
#'  \item{PHENO}{continuous phenotype (numeric)}
#' }
#' @references This data set was artificially created for the PXStools package.
#' @keywords datasets
#' @examples
#'
#' data(BINARY_DF)
#' head(BINARY_DF)
#'
"BINARY_DF"
