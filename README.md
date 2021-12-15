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
REM='L' #remove the response 'L' from our analysis 

#randomly sort data into three equal sized group, group C will contain individuals with a final predicted PXS
ss <- sample(1:3,size=nrow(CONT_DF),replace=TRUE,prob=c(1/5,1/5,3/5))
id_A<-CONT_DF$ID[ss==1]
id_B<-CONT_DF$ID[ss==2]
id_C<-CONT_DF$ID[ss==3]

```
Run XWAS: 
```R
XWAS_results=xwas(df=CONT_DF,X=XVAR,cov = COV,mod = 'lm',IDA = id_A,removes = REM)
head(XWAS_results)
#       Estimate Std..Error   t.value     Pr...t.. nrow.stored.          fdr
# VAR_1 493.86983   25.00103 19.753983 2.975498e-73          982 6.189037e-71
# VAR_2 266.26276   32.41371  8.214511 6.770669e-16          982 1.279656e-13
# VAR_3  74.47888   34.42183  2.163710 3.073039e-02          982 1.000000e+00
# VAR_4  78.78474   34.90476  2.257134 2.422168e-02          982 1.000000e+00
# VAR_5  35.32911   34.42199  1.026353 3.049814e-01          982 1.000000e+00
# VAR_6 479.30374   25.70496 18.646355 1.517820e-66          982 3.141887e-64

#obtain significant X's
sigx=row.names(XWAS_results)[which(XWAS_results$fdr<0.05)]
sigx[12:length(sigx)]=substr(sigx[12:length(sigx)],1,nchar(sigx[12:length(sigx)])-1) #remove levels and only keep name of variable
sigx=unique(sigx)

sigx
#[1] "VAR_1"   "VAR_2"   "VAR_6"   "VAR_8"   "VAR_9"   "VAR_10"  "VAR_14"  "VAR_16"  "VAR_17"  "VAR_18"  "VAR_22C" "VAR_22"  "VAR_25" 

```
Visualize results from XWAS: 
```R
manhattan_xwas(xdff = XWAS_results,pval = 'fdr',thresh = 0.05) #plots p values on -log10 scale
plot_coeff_xwas(xdff = XWAS_results,pval = 'fdr',coeff='Estimate',thresh = 0.05) #plots coefficients of signficant results, set all=TRUE to plot all results
```
![image](https://user-images.githubusercontent.com/54297194/146267701-afa47654-6b01-4c86-bc43-42f747c27d38.png)

Run PXS (with only signficant exposures): 
```R
PXSS=PXS(df=CONT_DF,X=sigx,cov=COV,removes = REM,mod = 'lm',IDA = id_A,IDB = id_B,IDC = id_C,seed=5)

[1] "intiating PXS procedure with 13 variables"
[1] "excluding individuals..."
[1] "914 individuals remain"
[1] "transformed responsetab"
[1] "LASSO step initiating..."
[1] "cross validated LASSO complete"
[1] "the  min lamda  is: 0.021807607315826"
[1] "11 variables remain after LASSO"
[1] "excluding individuals..."
[1] "930 individuals remain"
[1] "8 remain after FS iteration 1"
8 remain after final FS iteration, they are:  VAR_22 VAR_25 VAR_1 VAR_2 VAR_6 VAR_8 VAR_10 VAR_17 
[1] "0 individuals removed due to factor having a new level"

nrow(PXSS) #number of individuals with PXS
# 2831  

head(PXSS)
#   ID    SEX AGE COV_Q_OTHER COV_C_OTHER      PC_1      PC_2       PC_3      PC_4        VAR_1         VAR_2        VAR_6         VAR_8
# 4572 FEMALE  47    7.738711  CATEGORY_5 1.2674114 0.9203929  2.6974477 0.9914804  0.003945463 -0.0007426558 -0.009524965 -0.0007623313
# 2754   MALE  44   15.349081  CATEGORY_5 0.5550317 0.9852589 -2.4922787 2.0085929 -0.005833392 -0.0070134080  0.018185588 -0.0192416176
# 2678   MALE  45   15.081943  CATEGORY_4 0.5632954 1.1129693 -1.5017750 2.6902221 -0.004395898 -0.0047589039  0.009265629 -0.0152872580
# 3064   MALE  50    5.657369  CATEGORY_3 0.1701464 0.9345670  0.5661696 1.1395434 -0.001120768  0.0061214494 -0.007926178  0.0099030477
# 3976   MALE  66    1.470289  CATEGORY_2 1.1001742 0.4358340  0.8061108 1.1374945  0.012189827  0.0047791377  0.003133365  0.0117178749
# 4364 FEMALE  48   16.530377  CATEGORY_1 0.9128850 1.3486736  0.5913076 0.3542070 -0.036801092 -0.0207197737 -0.011525390  0.0050126760

#  VAR_10       VAR_17       VAR_22 VAR_25 PHENO     pred       
# -0.02206128 -0.028053092      C      B    95    115.59746
# -0.04251531  0.004848575      D      C    77    81.76173
#  0.00294643  0.018504878      C      B    95    116.12381
#  0.03748063  0.012897271      E      B   116    112.76533
# -0.01354368 -0.017855771      G      C   102    93.17732
# -0.01347574 -0.001232650      D      C    79    77.85163
```
If you would like to consider interactions in calculate PXS, please use the PXSgl function instead:
```R
PXSinter=PXSgl(df=CONT_DF,X=sigx,cov=COV,removes = REM,mod = 'lm',IDA = id_A,IDB = id_B,IDC = id_C,seed=5,fold=5)
 
# "intiating group lasso PXS procedure with 13 variables"
# "excluding individuals..."
# "914 individuals remain"
# "cross validation with 5 folds"
# "the  min lamda  is: 0.0050374302730364"
# "recalibrating model in group B..."
# "excluding individuals..."
# "3761 individuals remain"

nrow(PXSinter) #number of individuals with PXS
# 2831

head(PXSinter)
#    ID    SEX AGE COV_Q_OTHER COV_C_OTHER      PC_1      PC_2       PC_3      PC_4        VAR_1         VAR_2        VAR_6         VAR_8
# 4572 FEMALE  47    7.738711  CATEGORY_5 1.2674114 0.9203929  2.6974477 0.9914804  0.003945463 -0.0007426558 -0.009524965 -0.0007623313
# 2754   MALE  44   15.349081  CATEGORY_5 0.5550317 0.9852589 -2.4922787 2.0085929 -0.005833392 -0.0070134080  0.018185588 -0.0192416176
# 2678   MALE  45   15.081943  CATEGORY_4 0.5632954 1.1129693 -1.5017750 2.6902221 -0.004395898 -0.0047589039  0.009265629 -0.0152872580
# 3064   MALE  50    5.657369  CATEGORY_3 0.1701464 0.9345670  0.5661696 1.1395434 -0.001120768  0.0061214494 -0.007926178  0.0099030477
# 3976   MALE  66    1.470289  CATEGORY_2 1.1001742 0.4358340  0.8061108 1.1374945  0.012189827  0.0047791377  0.003133365  0.0117178749
# 4364 FEMALE  48   16.530377  CATEGORY_1 0.9128850 1.3486736  0.5913076 0.3542070 -0.036801092 -0.0207197737 -0.011525390  0.0050126760
         VAR_9      VAR_10       VAR_14       VAR_16       VAR_17 VAR_22 VAR_25 PHENO      pred
#  0.014674896 -0.02206128  0.028193331 -0.003030358 -0.028053092      C      B    95 115.34940
#  0.010574598 -0.04251531 -0.004591859 -0.021717874  0.004848575      D      C    77  86.25764
# -0.037871902  0.00294643 -0.073432151 -0.015424781  0.018504878      C      B    95 124.16328
# -0.001291723  0.03748063  0.002486285  0.009385328  0.012897271      E      B   116 106.97016
#  0.018172056 -0.01354368  0.030513644 -0.002688886 -0.017855771      G      C   102  89.48601
# -0.010092473 -0.01347574 -0.010005814 -0.007991927 -0.001232650      D      C    79  79.16278



