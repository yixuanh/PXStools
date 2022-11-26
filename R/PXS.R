#' @title PXS
#' @description  this function builds a polyexposure risk score
#' @source
#' @author Yixuan He, \email{yixuan_he@@hms.harvard.edu}

#' @param df the data frame input
#' @param X column name of significant exposure variables from XWAS
#' @param cov column name of covariates
#' @param mod type of model to run; 'lm' for linear regression, 'logistic' for logistic regression; 'cox' for Cox regression
#' @param IDA list of IDs to from XWAS procedure
#' @param IDB list of IDs for testing set
#' @param IDC list of IDs in the final prediction set
#' @param seed setting a seed
#' @param removes any exposure response, categorical or numerical, to remove from the analysis. This should be in the form of a list
#' @param fdr whether or not to adjust for multiple hypothesis correction
#' @param intermediate whether or not to save intermediate files
#' @param folds number of folds for glmnet cross validation, default is 10
#' @param alph the alpha value used in glmnet, alpha = 1 is assumed by default (lasso),
#' setting alpha = 0 for ridge, and anything in between 0 and 1 for elastic net.
#' please refer to glmnet documentation for more details

#' @export
#'
#' @importFrom magrittr %>%
#' @import naniar
#' @import glmnet
#' @import broom
#' @import ggplot2
#' @import logger
#'
PXS = function(df,
               X,
               cov,
               mod,
               IDA,
               IDB,
               IDC,
               seed,
               alph=1,
               folds=10,
               removes = NULL,
               intermediate = F) {


  `%notin%` <- Negate(`%in%`)
  set.seed(seed)

  #check to make sure each column has more than 1 unique value
  logger::log_info(paste('intiating PXS procedure with', length(X), 'variables'))


  dfA = df[which(df$ID %in% IDA), ]
  keep <- data.frame(dfA[, c(which(colnames(dfA) %in% c('PHENO', cov, X)))])

  if(mod=='cox'){
    keep <- data.frame(dfA[, c(which(colnames(dfA) %in% c('PHENO','TIME', cov, X)))])

  }
  one = which(sapply(keep, function(x)
    length(unique(x)) > 1) == FALSE)
  if (length(one) != 0) {
    logger::log_warn(paste(
      colnames(keep)[one],
      'has only one unique value, please double check or remove.'
    ))
    stop()
  }

  if (length(removes) != 0) {
    logger::log_info('excluding individuals...')
    b=apply(keep, 1, function(r) any(r %in% removes))
    if(length(which(b==TRUE))!=0){
      keep=keep[-which(b==TRUE),]
    }
    else{
      logger::log_info('no responses to remove')
    }
  }
  keep = na.omit(keep)
  logger::log_info(paste(nrow(keep),'individuals remain'))


  responsetab=function(dff){
    resptab=c()
    resptab2=c()
    nums = data.frame(sapply(colnames(dff), function(x)
      class(dff[[x]])))
    if(nrow(nums)!=2|ncol(dff)==2){
      nums=t(nums)
      logger::log_warn('transformed responsetab')
    }

    cati = colnames(nums)[which(nums[1,] != 'numeric')]
    if(length(cati)!=0){
    resptab = plyr::ldply(dff[which(colnames(dff)%in%cati)], function(x)
      t(rbind(names(table(
        x
      )), table(x))))
    resptab$varsrespon = paste(resptab[, 1], resptab[, 2], sep = '')

    colnames(resptab) = c('Var', 'Response', 'N','VR')
    resptab$Response=as.character(resptab$Response)
    }

    nums= colnames(nums)[which(nums[1,] == 'numeric')]
    if(length(nums)!=0){
      resptab2=cbind(nums,nrow(dff),nrow(dff),nums)
      colnames(resptab2)=c('Var', 'Response', 'N','VR')
    }

    resptab=rbind(resptab,resptab2)
    return(resptab)
  }

  Xtemp = keep[, which(colnames(keep) %in% X)]
  rt=responsetab(Xtemp)
  nn=which(rt$N==0)
  if(length(nn)!=0){
  rt=rt[-which(rt$N==0),]
  }

  ##REGULARIZATION
  ############

  x_vars = model.matrix(keep$PHENO ~ ., keep[, -(which(colnames(keep) ==
                                                           'PHENO'))]) ##change here#####

  y_var = keep$PHENO

  lambdav <- NULL

  if(alph==0){
    logger::log_info('ridge regression initiating...')
  }
  if(alph==1){
    logger::log_info('LASSO initiating...')
  }
  if(alph>0&alph<1){
    logger::log_info('elastic net initiating...')
  }
  if(alph<0|alph>1){
    logger::log_info('please use an alpha value between 0 and 1')
    break
  }
  if (mod == 'lm') {
    cv_output <- glmnet::cv.glmnet(x_vars, y_var,nfolds=folds, alpha=alph,)
    best_lamb <- cv_output$lambda.min
    logger::log_info ('cross validation complete')
  }

  if (mod == 'logistic') {
    cv_output <- glmnet::cv.glmnet(x_vars, y_var, nfolds=folds, alpha=alph, family = "binomial")
    best_lamb <- cv_output$lambda.min
    logger::log_info ('cross validation complete')
  }

  if (mod == 'cox') {
    y_var = cbind(keep$PHENO, keep$TIME)
    colnames(y_var) = c( 'status','time')
    y_var=as.matrix(y_var)
    x_vars = model.matrix(keep$PHENO ~ ., keep[, -which(colnames(keep)%in%c('PHENO', 'TIME'))])

    cv_output <- glmnet::cv.glmnet(x_vars, y_var,nfolds=folds, alpha=alph,family='cox')
    best_lamb<-cv_output$lambda.min
    logger::log_info (paste('cross validation complete'))

  }
  if (mod %notin% c('cox', 'lm', 'logistic')) {
    logger::log_warn('please specificy a regression model: lm, logsitic, or cox ')
  }

  logger::log_info(paste('the  min lamda  is:', best_lamb))

  if(mod=='lm'){
    lasso_best <- glmnet::glmnet(x_vars, y_var,  lambda = best_lamb, alpha=alph,)
  }
  if(mod=='logistic'){
    lasso_best <- glmnet::glmnet(x_vars, y_var,  lambda = best_lamb, alpha=alph,family = "binomial")
  }
  if(mod=='cox'){
    lasso_best <- glmnet::glmnet(x_vars, y_var,  lambda = best_lamb, alpha=alph,family='cox')

  }
  if (intermediate == TRUE) {
    saveRDS(lasso_best, 'regularization_best.rds')
  }

  tmp_coeffs <- coef(lasso_best)
  #save variables with non-zero coeffs
  M <-
    data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
  M <- unique(rt$Var[which(rt$VR %in% M$name)])

  if (length(M) == 0) {
    logger::log_warn('no variables remain after regularization')
    break
  }

  if (length(M) != 0) {
    logger::log_info(paste(length(M), 'variables remain after regularization'))
  }

  ################
  #stepwise procedure
  dfB = df[which(df$ID %in% IDB), c(which(colnames(df) %in% c('PHENO', cov, M)))]

  if(mod=='cox'){
    dfB = df[which(df$ID %in% IDB), c(which(colnames(df) %in% c('PHENO', 'TIME',cov, M)))]
  }

  if (length(removes) != 0) {
    logger::log_info('excluding individuals...')
    b=apply(dfB, 1, function(r) any(r %in% removes))
    if(length(which(b==TRUE))!=0){
      dfB=dfB[-which(b==TRUE),]
    }
    else{
      logger::log_info('no responses to remove')
    }
    logger::log_info(paste(nrow(dfB),'individuals remain'))

  }

  if (mod == 'lm') {
    B_temp <-
      data.frame(dfB[, c(which(colnames(dfB) %in% c('PHENO', cov, M)))])
    B_temp<-na.omit(B_temp)
    fit <- broom::tidy(lm(PHENO ~ 0 + ., data = B_temp))
  }

  if (mod == 'logistic') {
    B_temp <-
      data.frame(dfB[, c(which(colnames(dfB) %in% c('PHENO', cov, M)))])
    B_temp<-na.omit(B_temp)
    fit <- broom::tidy(glm(PHENO ~ 0 + ., data = B_temp, family = 'binomial'))
  }

  if (mod == 'cox') {
    B_temp <-
      data.frame(dfB[, c(which(colnames(dfB) %in% c('PHENO', 'TIME', cov, M)))])
    B_temp<-na.omit(B_temp)
    fit <-
      broom::tidy(survival::coxph(survival::Surv(TIME, PHENO) ~ 0 + ., data = B_temp))
  }

  sig = fit$term[which(fit$p.value < 0.05)]
  sig = unique(rt$Var[which(rt$VR %in% sig)])

  logger::log_info(paste(length(sig), 'remain after BackS iteration 1'))

  if (length(sig) == 0) {
    break
  }

  initial = sig

  for (i in 1:(length(M) * 5)) {
    if (mod == 'lm') {
      B_temp <-
        data.frame(dfB[, c(which(colnames(dfB) %in% c('PHENO', cov, sig)))])
      B_temp<-na.omit(B_temp)
      fit <- broom::tidy(lm(PHENO ~ 0 + ., data = B_temp))
    }

    if (mod == 'logistic') {
      B_temp <-
        data.frame(dfB[, c(which(colnames(dfB) %in% c('PHENO', cov, sig)))])
      B_temp<-na.omit(B_temp)
      fit <- broom::tidy(glm(PHENO ~ 0 + ., data = B_temp, family = 'binomial'))
    }

    if (mod == 'cox') {
      B_temp <-
        data.frame(dfB[, c(which(colnames(dfB) %in% c('PHENO', 'TIME', cov, sig)))])
      B_temp<-na.omit(B_temp)
      fit <-
        broom::tidy(survival::coxph(survival::Surv(TIME, PHENO) ~ 0 + ., data = B_temp))
    }

    sig = fit$term[which(fit$p.value < 0.05)]
    sig = unique(rt$Var[which(rt$VR %in% sig)])

    if (length(setdiff(sig, initial)) == 0) {
      logger::log_info(cat(length(sig),"remain after final BackS iteration, they are: ", sig,"\n",sep=" "))
      break
    }

    logger::log_info(paste(length(sig), 'remain after BackS iteration', i))

  }
  #################
  ##FINAL PREDICTION MODEL
  if(mod=='lm'|mod=='logistic'){
    dfBC=df[which(df$ID%in%c(IDB,IDC)),c(which(colnames(df) %in% c('ID','PHENO', cov, sig)))]
  }
  if(mod=='cox'){
    dfBC=df[which(df$ID%in%c(IDB,IDC)),c(which(colnames(df) %in% c('ID','PHENO', 'TIME',cov, sig)))]
  }
  dfBC=na.omit(dfBC)


  if (length(removes) != 0) {
    b=apply(dfBC, 1, function(r) any(r %in% removes))
    dfBC=dfBC[-which(b==TRUE),]

  }



  #final model fit
  if (mod == 'lm') {
    B_temp <-
      data.frame(dfBC[which(dfBC$ID%in%IDB), c(which(colnames(dfBC) %in% c('PHENO', cov, sig)))])
    B_temp<-na.omit(B_temp)
    fit <- lm(PHENO ~ 0 + ., data = B_temp)

  }

  if (mod == 'logistic') {
    B_temp <-
      data.frame(dfBC[which(dfBC$ID%in%IDB), c(which(colnames(dfBC) %in% c('PHENO', cov, sig)))])
    B_temp<-na.omit(B_temp)
    fit <- glm(PHENO ~ 0 + ., data = B_temp, family = 'binomial')
  }

  if (mod == 'cox') {
    B_temp <-
      data.frame(dfBC[which(dfBC$ID%in%IDB), c(which(colnames(dfBC) %in% c('PHENO', 'TIME', cov, sig)))])
    B_temp<-na.omit(B_temp)
    fit <-
      survival::coxph(survival::Surv(TIME, PHENO) ~ 0 + ., data = B_temp)
  }
  coeffs=broom::tidy(fit)

  if(intermediate==TRUE){
    write.csv(coeffs,'coefficients.csv')
  }

  C_temp <-dfBC[which(dfBC$ID%in%IDC),]

  templength=nrow(C_temp)
  for(s in sig){
    if(data.class(B_temp[,which(colnames(B_temp)%in%s)])=='factor'){
      id <- which(!(C_temp[,which(colnames(C_temp)%in%s)] %in% levels(droplevels(B_temp[,which(colnames(B_temp)%in%s)]))))
    }
    if(data.class(B_temp[,which(colnames(B_temp)%in%s)])=='character'){
      id <- which(!(C_temp[,which(colnames(C_temp)%in%s)] %in% unique(B_temp[,which(colnames(B_temp)%in%s)])))
    }

    if(length(id)!=0){
      C_temp=C_temp[-id,]
    }
  }
  logger::log_warn(paste((templength-nrow(C_temp)),'individuals removed due to factor having a new level'))


  if(mod=='lm'){
    C_temp$PXS=predict(fit,C_temp)
  }
  if (mod=='logistic'){
    C_temp$PXS=predict(fit,C_temp,type='response')
  }
  if (mod=='cox'){
    C_temp$PXS=predict(fit,C_temp,type='risk')
  }

  return(C_temp)
}


