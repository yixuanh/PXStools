# PXStools
PXStools is a R software package to provide tools for conducting exposure association studies. The accompanying paper can be found at: 

## Installation
The package can be directly downloaded from R: 
```R
install.packages("devtools")
```

The development version from GitHub can be downloaded with: 
```R
devtools::install_github("yixuanh/PXStools")
```

## Functions
The package contains four functions: 

``xwas`` :conduct an exposure-wide association study (XWAS) across given exposure for a single phenotype. Please refer to https://doi.org/10.1371/journal.pone.0010746 for more details.

``plot_coeff_xwas`` : plots the beta coefficients from the XWAS results. 

``manhattan_xwas`` : plots the p values from the XWAS results on a -log scale. Be careful fof any associations that have a p value of close to zero as that will approach infinity in the -log scale. 

``PXS`` : conducts the LASSO-based selection procedure on a set of given exposures to build a poly-exposure risk score for a single phenotype. It is recommended that the inputed exposures for ``PXS`` are the signficant associations from the XWAS to minimize sample loss. 

## Options 
In both the ``xwas`` and ``PXS`` functions, the user can input any set of exposures of interest. It is also possible to run different types of regression analysis including ``lm`` for linear models, ``logistic`` for binary phenotypes, and ``cox`` for cox regression. The user can choose a set of covariates to adjust for at each stage of the analysis as well as which exposure factors to remove from the analysis. 

## Requirements
The input data frame must have the following columns: 
``ID``: ID of indivduals in dataframe
``PHENO``: phenotype of interest (binary (0/1) or continuous)
if running survival analysis, it must also have 
``TIME``: time to event or censoring
In addition to the the final prediction group, two other groups are needed to train the model. 

## Example

This will be an example using the ``CONT_DF.RData`` dataset provided in the package. The CONT_DF dataset contains the individual ID, sex, gender, continuous and categorical variables, and a continuous phenotype. We will use SEX, AGE, COV_Q_OTHER, and COV_C_OTHER, as our covariate. The initial set of exposures that we are interested in are VAR_1 through VAR_33. We will be using a linear model. 


Store variable names: 
```R
set.seed(7)

COV=colnames(CONT_DF)[2:5] #covariate names
XVAR=colnames(CONT_DF)[16:48] #exposure names
REM='B' #remove the response 'B' from our analysis 

#randomly sort data into three equal sized group, group C will contain individuals with a final predicted PXS
ss <- sample(1:3,size=nrow(CONT_DF),replace=TRUE,prob=c(1/3,1/3,1/3))
id_A<-CONT_DF$ID[ss==1]
id_B<-CONT_DF$ID[ss==2]
id_C<-CONT_DF$ID[ss==3]

```
Run XWAS: 
```R
XWAS_results=xwas(df=CONT_DF,X=XVAR,cov = COV,mod = 'lm',IDA = id_A,removes = REM)
head(XWAS_results)
#       Estimate Std..Error   t.value      Pr...t.. nrow.stored.           fdr
#VAR_1 501.99790   18.86722 26.606885 4.129432e-130         1661 7.391683e-128
#VAR_2 270.70631   24.82412 10.904970  8.711046e-27         1661  1.376345e-24
#VAR_3  71.63699   26.65408  2.687656  7.267780e-03         1661  1.000000e+00
#VAR_4  93.61497   25.91737  3.612055  3.128642e-04         1661  4.724249e-02
#VAR_5  55.64746   26.38567  2.109003  3.509444e-02         1661  1.000000e+00
#VAR_6 446.23995   20.25330 22.032946  1.570251e-94         1661  2.795047e-92

#obtain significant X's
sigx=row.names(XWAS_results)[which(XWAS_results$fdr<0.05)]
sigx[12:length(sigx)]=substr(sigx[12:length(sigx)],1,nchar(sigx[12:length(sigx)])-1) #remove levels and only keep name of variable
sigx=unique(sigx)

sigx
# [1] "VAR_1"  "VAR_2"  "VAR_4"  "VAR_6"  "VAR_8"  "VAR_9"  "VAR_10" "VAR_14" "VAR_16" "VAR_17"
# [11] "VAR_18" "VAR_22" "VAR_25"

```
Visualize results from XWAS: 
```R
manhattan_xwas(xdff = XWAS_results,pval = 'fdr',thresh = 0.05) #plots p values on -log10 scale
plot_coeff_xwas(xdff = XWAS_results,pval = 'fdr',coeff='Estimate',thresh = 0.05) #plots coefficients of signficant results, set all=TRUE to plot all results
```
![Untitled 3 001](https://user-images.githubusercontent.com/54297194/140362768-f24f7990-0e37-47b4-91ea-675f62d64bec.png)

 
Run PXS (with only signficant exposures): 
```R
PXSS=PXS(df=CONT_DF,X=sigx,cov=COV,removes = REM,mod = 'lm',IDA = id_A,IDB = id_B,IDC = id_C,seed=5)

# "intiating PXS procedure with 13 variables"
# "excluding individuals..."
# "1019 individuals remain"
# "transformed responsetab"
# "LASSO step initiating..."
# "cross validated LASSO complete"
# "the  min lamda  is: 0.0311043365158031"
# "12 variables remain after LASSO"
# "excluding individuals..."
# "74 individuals remain"
# "3 remain after FS iteration 1"
3 remain after final FS iteration, they are:  VAR_22 VAR_25 VAR_9 
# "0 individuals removed due to factor having a new level"

nrow(PXSS) #number of individuals with PXS
# 1026  

head(PXSS)
#      ID    SEX AGE COV_Q_OTHER COV_C_OTHER        VAR_9 VAR_22 VAR_25 PHENO     pred
# 7  3976   MALE  66    1.470289  CATEGORY_2  0.018172056      G      C   102 94.68780
# 10  947 FEMALE  42   11.257224  CATEGORY_4  0.018919003      G      Q    87 97.59953
# 15 3605 FEMALE  46   14.650950  CATEGORY_1  0.012576732      D      K    89 83.15233
# 17 3979   MALE  62    8.133292  CATEGORY_2 -0.002686714      E      D    91 97.28691
# 21 3243   MALE  45    9.615657  CATEGORY_2 -0.018024730      K      I    87 94.38959
# 27 4345   MALE  47    1.258407  CATEGORY_4 -0.009366692      D      A    69 74.13231
```
     



