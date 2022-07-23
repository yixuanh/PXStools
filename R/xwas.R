#' @title xwas
#' @description  this is a function to run an exposure wide association study (XWAS) with any phenotype and a set of exposure variables

#' @author Yixuan He, \email{yixuan_he@@hms.harvard.edu}

#' @param df the data frame input
#' @param X column name of exposure variables to run XWAS
#' @param cov column name of covariates
#' @param mod type of model to run; 'lm' for linear regression, 'logistic' for logistic regression; 'cox' for Cox regression
#' @param IDA list of IDs to include in XWAS
#' @param removes any exposure response, categorical or numerical, to remove from XWAS. This should be in the form of a list
#' @param adjust method for adjusting for multiple comparison, see ?p.adjust to see other options
#' @param intermediate saves an intermediate file containing the coefficients of covariates. Default is False
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @import naniar
#' @import glmnet
#' @import broom
#' @import ggplot2
#' @import logger


xwas = function(df,
                X,
                cov,
                mod,
                IDA,
                intermediate=F,
                removes = NULL,
                adjust = 'BY') {

  `%notin%` <- Negate(`%in%`)

  df = df[which(df$ID %in% IDA), ]
  ci <- which(colnames(df) %in% X)
  df$place = 0
  mat = c()
  if(intermediate==T){
    coeffs=c()
  }

  pb <- txtProgressBar(0, length(ci), style = 3)
  stepi = 0


  for (i in ci) {
    stored <-
      data.frame(df[, c(which(colnames(df) %in% c('PHENO', cov)), which(colnames(df) ==
                                                                          'place'), i)])

    if (mod == 'cox') {
      stored <-
        data.frame(df[, c(which(colnames(df) %in% c('PHENO', 'TIME', cov)), which(colnames(df) ==
                                                                                    'place'), i)])
    }

    #omit missing and unwanted X var
    stored <- na.omit(stored)
    r = which(stored[, ncol(stored)] %in% removes)
    if (length(r) != 0) {
      stored <- stored[-r, ]
    }

    #skip if X only has one value
    if(nrow(unique(stored[ncol(stored)]))==1){
      log_warn(paste(colnames(stored)[ncol(stored)],'has only one value, skipped'))
      next
    }

    #run regression models

    if (mod == 'lm') {
      fit <- lm(PHENO ~ 0 + ., data = stored)
    }
    if (mod == 'logistic') {
      fit <- glm(PHENO ~ 0 + ., data = stored, family = 'binomial')
    }
    if (mod == 'cox') {
      fit <-
        survival::coxph(survival::Surv(TIME, PHENO) ~ 0 + ., data = stored)
    }
    if (mod %notin% c('cox', 'lm', 'logistic')) {
      log_warn('please specificy a regression model: lm, logsitic, or cox ')
    }

    #get index for which coeffs to extract
    temp = data.frame(coef(fit))
    tempp = temp[which(row.names(temp) == 'place') + 1, ]

    #case for if there is NA in the first categorical factor, will iterate until non-NA
    if (is.na(tempp) == TRUE) {
      for (k in 2:(nrow(temp) - which(row.names(temp) == 'place'))) {
        tempp = temp[which(row.names(temp) == 'place') + k, ]
        if (is.na(tempp) == FALSE) {
          break
        }
      }

    }

    summary_fit <- data.frame(rbind(summary(fit)$coefficients),nrow(stored))

    if(intermediate==T){
      if(length(cov)==0){
        log_warn('there are no covariates.')
        break
      }
      if(mod=='cox'){
        coefftemp <- summary_fit[1:(which(row.names(summary_fit) == 'place')-1), ]
        rownames(coefftemp) <-paste(colnames(df)[i],rownames(coefftemp),sep='_')
      }
      else{
        coefftemp <- summary_fit[1:which(summary_fit[, 1] == tempp), ]
        rownames(coefftemp) <-paste(colnames(df)[i],rownames(coefftemp),sep='_')
      }
      coeffs = rbind(coeffs,coefftemp)
    }
    summary_fit <-
      summary_fit[(which(summary_fit[, 1] == tempp):nrow(summary_fit)), ]
    mat = rbind(mat, summary_fit)

    stepi = stepi + 1
    setTxtProgressBar(pb, stepi)

  }

  mat = data.frame(mat)

  mat$fdr = p.adjust(mat[, ncol(mat)-1],method=adjust)

  if(intermediate==T){
    write.csv(coeffs,'XWAS_covariate_coefficients.csv')
  }
  return(mat)
}
