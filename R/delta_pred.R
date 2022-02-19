#' @title delta_pred
#' @description  this is a function that claculates the change in predictive ability between two models. For linear mdoels,
#' a change in R2 wil be reported; for logistic regression models, a change in AUC will be reported; for Cox regresison modelx,
#' a change in C index will be reported. The column name of the Y variable must be "PHENO". For Cox regression models, the time to event column name must be "TIME"/

#' @author Yixuan He, \email{yixuan_he@@hms.harvard.edu}

#' @param df the data frame inpt
#' @param xvarsA column name of variables to include in first model
#' @param xvarsB column name of variables to include in second  model
#' @param mod type of model to run; 'lm' for linear regression, 'logistic' for logistic regression; 'cox' for Cox regression
#' @param boot number of boostrap samples, default is 100
#' @export
#'
#' @importFrom magrittr %>%
#' @import naniar
#' @import glmnet
#' @import broom
#' @import ggplot2
#' @import logger
#' @import boot
#' @import pROC
#'
delta_pred = function(df,
                      varsA,
                      varsB,
                      mod,
                      boot=100){
  if (mod=='lm'){
      fA=as.formula(paste('PHENO',paste(varsA,collapse='+'),sep='~'))
      fB=as.formula(paste('PHENO',paste(varsB,collapse='+'),sep='~'))

      difffunc <- function(data, i){
          d1 <- summary(lm(fA,data[i,]))$r.squared
          d2 <- summary(lm(fB,data[i,]))$r.squared
          return(d2-d1)
      }
      difffunc1 <- function(data, i){
        d1 <- summary(lm(fA,data[i,]))$r.squared
        return(d1)
      }
      difffunc2 <- function(data, i){
        d2 <- summary(lm(fB,data[i,]))$r.squared
        return(d2)
      }
      f1=boot(df,difffunc1,R=boot)
      f2=boot(df,difffunc2,R=boot)

      diffs=invisible(boot(df,difffunc,R=boot))
      log_info(paste('change in R2: ',round(diffs$t0,3),' (',round(quantile(diffs$t,c(0.025,0.975))[1],3),', ',round(quantile(diffs$t,c(0.025,0.975))[2],3),')',sep=''))
      log_info(paste('first model R2: ', paste(round(f1$t0,3),' (',round(quantile(f1$t,c(0.025,0.975))[1],3),', ',round(quantile(f1$t,c(0.025,0.975))[2],3),')',sep='')))
      log_info(paste('second model R2: ', paste(round(f2$t0,3),' (',round(quantile(f2$t,c(0.025,0.975))[1],3),', ',round(quantile(f2$t,c(0.025,0.975))[2],3),')',sep='')))
  }
  if (mod=='logistic'){
      fA=as.formula(paste('PHENO',paste(varsA,collapse='+'),sep='~'))
      fB=as.formula(paste('PHENO',paste(varsB,collapse='+'),sep='~'))
      difffunc <- function(data, i){
        d=data[i,]
        fit1 <- glm(fA, d,family='binomial')
        fit2 <- glm(fB, d,family='binomial')
        prob1=predict(fit1,type=c("response"))
        prob2=predict(fit2,type=c("response"))

        d$prob1=prob1
        d$prob2=prob2
        g <- roc(PHENO~ prob1, data = d)
        g2 <- roc(PHENO~ prob2, data = d)

        return(g2$auc-g$auc)
      }
      difffunc1 <- function(data, i){
        d=data[i,]
        fit1 <- glm(fA, d,family='binomial')
        prob1=predict(fit1,type=c("response"))

        d$prob1=prob1
        g <- roc(PHENO~ prob1, data = d)
        return(g$auc)
      }
      difffunc2 <- function(data, i){
        d=data[i,]
        fit2 <- glm(fB, d,family='binomial')
        prob2=predict(fit2,type=c("response"))

        d$prob2=prob2
        g2 <- roc(PHENO~ prob2, data = d)
        return(g2$auc)
      }
      diffs=invisible(boot(df,difffunc,R=boot))
      f1=invisible(boot(df,difffunc1,R=boot))
      f2=invisible(boot(df,difffunc2,R=boot))
      log_info(paste('change in AUC: ',round(diffs$t0,3),' (',round(quantile(diffs$t,c(0.025,0.975))[1],3),', ',round(quantile(diffs$t,c(0.025,0.975))[2],3),')',sep=''))
      log_info(paste('first model AUC: ', paste(round(f1$t0,3),' (',round(quantile(f1$t,c(0.025,0.975))[1],3),', ',round(quantile(f1$t,c(0.025,0.975))[2],3),')',sep='')))
      log_info(paste('second model AUC: ', paste(round(f2$t0,3),' (',round(quantile(f2$t,c(0.025,0.975))[1],3),', ',round(quantile(f2$t,c(0.025,0.975))[2],3),')',sep='')))
 }
  if (mod=='cox'){
      fA=as.formula(paste('survival::Surv(TIME, PHENO)',paste(0,paste(varsA,collapse='+'),sep='+'),sep='~'))
      fB=as.formula(paste('survival::Surv(TIME, PHENO)',paste(0,paste(varsB,collapse='+'),sep='+'),sep='~'))
      difffunc <- function(data, i){
          d1 <- summary(survival::coxph(fA,data[i,]))$concordance[1]
          d2 <- summary(survival::coxph(fB,data[i,]))$concordance[1]
          return(d2-d1)
      }
      difffunc1 <- function(data, i){
        d1 <- summary(survival::coxph(fA,data[i,]))$concordance[1]
        return(d1)
      }
      difffunc2 <- function(data, i){
        d2 <- summary(survival::coxph(fB,data[i,]))$concordance[1]
        return(d2)
      }
      diffs=boot(df,difffunc,R=boot)
      f1=boot(df,difffunc1,R=boot)
      f2=boot(df,difffunc2,R=boot)
      log_info(paste('change in C index: ',round(diffs$t0,3),' (',round(quantile(diffs$t,c(0.025,0.975))[1],3),', ',round(quantile(diffs$t,c(0.025,0.975))[2],3),')',sep=''))
      log_info(paste('first model C index: ', paste(round(f1$t0,3),' (',round(quantile(f1$t,c(0.025,0.975))[1],3),', ',round(quantile(f1$t,c(0.025,0.975))[2],3),')',sep='')))
      log_info(paste('second model C index: ', paste(round(f2$t0,3),' (',round(quantile(f2$t,c(0.025,0.975))[1],3),', ',round(quantile(f2$t,c(0.025,0.975))[2],3),')',sep='')))

  }
}
