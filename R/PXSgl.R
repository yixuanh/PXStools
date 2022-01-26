#' @title PXSgl
#' @description  this function builds a polyexposure risk score considering interactions between exposures using the group lasso method
#' @source
#' @author Yixuan He, \email{yixuan_he@@hms.harvard.edu}

#' @param df the data frame inpt
#' @param X column name of signficant exposure variables from XWAS
#' @param cov column name of covariates
#' @param mod type of model to run; 'lm' for linear regression, 'logistic' for logistic regression; 'cox' for Cox regression
#' @param IDA list of IDs to from XWAS procedure
#' @param IDB list of IDs for testing set
#' @param IDC list of IDs in the final prediction set
#' @param seed setting a seed
#' @param removes any exposure response to remove from XWAS, in the form of a list
#' @param fdr whether or not to adjust for multiple hypothesis correctin
#' @param intermediate whether or not to save intermediate files
#' @param folds number of folds for the cross validation step, default is 10.

#' @export
#'
#' @importFrom magrittr %>%
#' @import naniar
#' @import glmnet
#' @import broom
#' @import ggplot2
#' @import glinternet

#'
PXSgl = function(df,
               X,
               cov,
               mod,
               IDA,
               IDB,
               IDC,
               seed,
               removes = NULL,
               intermediate = F,
               folds=10) {


  `%notin%` <- Negate(`%in%`)
  set.seed(seed)

  #check to make sure each column has more than 1 unique value
  print(paste('intiating group lasso PXS procedure with', length(X), 'variables'))

  df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],as.factor) #change all non-numeric columns to factor form

  dfA = df[which(df$ID %in% IDA), ]
  keep <- data.frame(dfA[, c(which(colnames(dfA) %in% c('PHENO', cov, X)))])

  one = which(sapply(keep, function(x)
    length(unique(x)) > 1) == FALSE)
  if (length(one) != 0) {
    print(paste(
      colnames(keep)[one],
      'has only one unique value, please double check or remove.'
    ))
    stop()
  }

  if (length(removes) != 0) {
    print('excluding individuals...')
    b=apply(keep, 1, function(r) any(r %in% removes))
    if(length(which(b==TRUE))!=0){
      keep=keep[-which(b==TRUE),]
    }
    else{
      print('no responses to remove')
    }
  }
  keep = na.omit(keep)
  print(paste(nrow(keep),'individuals remain'))

  ###GROUP LASSO WITH GLINTERNET
  #code adapted from tutorial https://strakaps.github.io/post/glinternet/
  Y=keep$PHENO
  keep=keep %>% dplyr::select(-PHENO)

  i_num <- sapply(keep, is.numeric)
  X <- keep

  # get the numLevels vector containing the number of categories
  numLevels <- X %>% sapply(nlevels)
  numLevels[numLevels==0] <- 1

  # make the categorical variables take integer values starting from 0
  X[, !i_num] <- apply(X[, !i_num], 2, function(col) as.integer(as.factor(col)) - 1)

  print(paste('cross validation with',folds,'folds'))

  if(mod %in% c('cox','logistic')){
    cv_fit <- glinternet::glinternet.cv(X,Y, nFolds =  folds,numLevels,family='binomial')
  }
  if(mod %in% c('lm')){
    cv_fit <- glinternet::glinternet.cv(X,Y, nFolds =  folds,numLevels)
  }

  lamb <- which.min(cv_fit$cvErr) #picking lambda with min error
  print(paste('the  min lamda  is:', cv_fit$lambda[lamb]))

  coefs <- coef(cv_fit$glinternetFit)[[lamb]]

  ##names of main/interaction effects
  ###names
  idx_num <- (1:length(i_num))[i_num]
  idx_cat <- (1:length(i_num))[!i_num]
  catnames=names(numLevels)[idx_cat[coefs$mainEffects$cat]] #categorical
  contnames=names(numLevels)[idx_num[coefs$mainEffects$cont]] #continuous

  if(length(catnames)==0){ #if no cat variables, set it to empty list
    catnames=c()
  }
  if(length(contnames)==0){ #if no numeric variables, set it to empty list
    contnames=c()
  }

  if(length(catnames) ==0 & length(contnames)==0){
    print ('group lasso did not select any variables, procedure terminated')
    stop()
  }

  mains=coefs$mainEffects #get main interactions only (gives indices of cat and cont variable)
  interacts=coefs$interactions #interactions (in pairs [i,j])

  ##################
  ## re-calibrate in group B
  print('recalibrating model in group B...')

  keepBC = df[which(df$ID %in% c(IDB,IDC)), c(which(colnames(df) %in% c('ID','PHENO', catnames,contnames,cov)))]

  if(mod=='cox'){
    keepBC = df[which(df$ID %in% c(IDB,IDC)), c(which(colnames(df) %in% c('ID','PHENO','TIME', catnames,contnames,cov)))]
  }
  if (length(removes) != 0) {
    print('excluding individuals...')
    b=apply(keepBC, 1, function(r) any(r %in% removes))
    if(length(which(b==TRUE))!=0){
      keepBC=keepBC[-which(b==TRUE),]
    }
    else{
      print('no responses to remove')
    }
  }
  keepBC = na.omit(keepBC)
  print(paste(nrow(keepBC),'individuals remain'))

  one = which(sapply(keepBC, function(x)
    length(unique(x)) > 1) == FALSE)
  if (length(one) != 0) {
    print(paste(
      colnames(keepBC)[one], #####may need to only have just one here not this whole thing
      'has only one unique value, removed in calibration step.'
    ))
  }


  catcat=data.frame(interacts$catcat)
  contcont=data.frame(interacts$contcont)
  catcont=data.frame(interacts$catcont)

  if(length(catcat)!=0){
  catcat$one_name=names(numLevels)[idx_cat[catcat[,1]]]
  catcat$two_name=names(numLevels)[idx_cat[catcat[,2]]]
  }

  if(length(contcont)!=0){
  contcont$one_name=names(numLevels)[idx_num[contcont[,1]]]
  contcont$two_name=names(numLevels)[idx_num[contcont[,2]]]
  }

  if(length(catcont)!=0){
  catcont$one_name=names(numLevels)[idx_cat[catcont[,1]]]
  catcont$two_name=names(numLevels)[idx_num[catcont[,2]]]
  }

  interacts=rbind(catcat,contcont,catcont)
  interacts$interact=paste(interacts$one_name,interacts$two_name,sep='*')

  if(length(one!=0)){
    interacts=interacts[!apply(interacts, 1, function(r) any(r %in%  names(one))),] #remove any interactions with variable that has only one response
  }

  mains=c(catnames,contnames)
  if(length(one!=0)){
   mains=mains[which(mains%notin%names(one))] #remove any mains with variable that has only one response
  }
  form_inter=paste(interacts$interact,collapse = '+')
  form_main=paste(mains,collapse = '+')

  cov=paste(cov,collapse = '+')
  if(mod%in%c('lm','logistic')){
    f=as.formula(paste('PHENO',paste(cov,form_main,form_inter,sep='+'),sep='~'))
  }

  if(mod%in%c('cox')){
    f=as.formula(paste('survival::Surv(TIME, PHENO)',paste(0,form_main,form_inter,sep='+'),sep='~'))
  }

  keepB <-keepBC[which(keepBC$ID%in%IDB),]
  keepB$ID=NULL

  if(mod=='lm'){
    model=lm(f,data=keepB)
  }
  if(mod=='logistic'){
    model=glm(f,data=keepB,family='binomial')
  }
  if(mod=='cox'){
    model=survival::coxph(f,data=keepB)
  }

  dfC <-keepBC[which(keepBC$ID%in%IDC),]

  if(mod=='lm'){
    dfC$pred=predict(model,dfC)
  }
  if (mod=='logistic'){
    dfC$pred=predict(model,dfC,type='response')
  }
  if (mod=='cox'){
    dfC$pred=predict(model,dfC,type='risk')
  }

  return(dfC)

}


