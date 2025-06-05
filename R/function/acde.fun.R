library(OpenMx)
# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA
# Load Data
# this script is based on H.M CTD script https://hermine-maes.squarespace.com/

acde.fun <- function(
    mzData=NULL,
    dzData=NULL,
    selVars = c("T1_p","T2_p"),
    name = "ADCE",
    ...,
    
    # Phenotype mean
    M =  1,
    
    # Phenotype variances
    # by default fit an AE model
    ap.free = T, # Freely estimate a? (Focal Phenotype)
    cp.free = F, # Freely estimate c? (Focal Phenotype)
    dp.free = F, # Freely estimate d? (Focal Phenotype)
  
    ap.val  = .30, # Starting value for a
    cp.val  = .0, # Starting value for c
    dp.val  = .0, # Starting value for d
    ep.val  = .70  # Starting value for e
){
  # ----------------------------------------------------------------------------------------------------------------------
  # Create Algebra for expected Mean Matrices
  meanG <- mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=M, labels="mean", name="meanG" )
 
   # Create Matrices for Variance Components
  covA <- mxMatrix( type="Symm", nrow=1, ncol=1, free=ap.free, values=ap.val, label="VA11", name="VA" )
  covC <- mxMatrix( type="Symm", nrow=1, ncol=1, free=cp.free, values=cp.val, label="VC11", name="VC" )
  covD <- mxMatrix( type="Symm", nrow=1, ncol=1, free=dp.free, values=dp.val, label="VD11", name="VD" )
  covE <- mxMatrix( type="Symm", nrow=1, ncol=1, free=TRUE,    values=ep.val, label="VE11", name="VE" )
  
  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covP <- mxAlgebra( expression= VA+VD+VC+VE, name="V" )
  covMZ <- mxAlgebra( expression= VA+VD+VC, name="cMZ" )
  covDZ <- mxAlgebra( expression= 0.5%x%VA+ 0.25%x%VD + VC, name="cDZ" )
  expCovMZ <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
  expCovDZ <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )
  
  # Create Data Objects for Multiple Groups
  dataMZ <- mxData( observed=mzData, type="raw" )
  dataDZ <- mxData( observed=dzData, type="raw" )
  # Create Expectation Objects for Multiple Groups
  expMZ <- mxExpectationNormal( covariance="expCovMZ", means="meanG", dimnames=selVars )
  expDZ <- mxExpectationNormal( covariance="expCovDZ", means="meanG", dimnames=selVars )
  funML <- mxFitFunctionML()
  # Create Model Objects for Multiple Groups
  pars <- list( meanG, covA, covD, covC, covE, covP )
  modelMZ <- mxModel( pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
  modelDZ <- mxModel( pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
  multi <- mxFitFunctionMultigroup( c("MZ","DZ") )
  # Create Algebra for Unstandardized and Standardized Variance Components
  rowUS <- c('US')
  colUS <- rep(c('varP','varP_A','varP_D','varP_C','varP_E'),each=1)
  estUS <- mxAlgebra( expression=cbind(V,VA/V,VD/V,VC/V,VE/V), name="US", dimnames=list(rowUS,colUS) )
  # Create Confidence Interval Objects
  ciADE <- mxCI( "US[1,1:5]" )
  # Build Model with Confidence Intervals
  modelADE <- mxModel( "oneADEvc", pars, modelMZ, modelDZ, multi, estUS, ciADE )
}