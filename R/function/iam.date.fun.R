# iAM-DATE function
message("Version 30.05.2025")
# Written by Hans Fredrik Sunde & Giacomo Bignardi
# Adapted from the iAM-ACE by Hans Fredrik Sunde
# Partner of twins model
# (MZ and DZ twins, their spouses and one sibling)
# With both iAM and dAM
# With MZ group and DZ group
# mxOption(key="Default optimizer", value="NPSOL")

iam.date.fun <- function(
    # data fro the two groups
    mzData=NULL,
    dzData=NULL,
    
    # names of the variable in the two groups (exact order! names can be changed)
    selVars = c("S1_p","T1_p","FS1_p","T2_p","S2_p"),
    
    # name of the model
    name = "iAM_DATE",
    
    # Focal Phenotype
    a1.free = T, # Freely estimate a? (Focal Phenotype)
    c1.free = F, # Freely estimate c? (Focal Phenotype)
    t1.free = T, # Freely estimate t? (Focal Phenotype)
    d1.free = T, # Freely estimate d? (Focal Phenotype)
    n1.free = F, # Freely estimate n? (Focal Phenotype)
    e1.free = T, # Freely estimate e  (Focal Phenotype)
    a1.val  = .7, # Starting value for a
    c1.val  = .0, # Starting value for c
    t1.val  = .0, # Starting value for t
    d1.val  = .0, # Starting value for d
    n1.val  = .0, # Starting value for n
    e1.val  = .7, # Starting value for e
    
    # Assortative mating options
    mu.free = TRUE, # Freely estimate assortment strength (mu)
    mu.val = .3, # starting value for mu
    
    U.f = 1, # Intergenerational equilibrium (set to 0 if first generation of assortment, 1 if equilibirum)
    indirect.assortment = FALSE, # If FALSE, model will equate sorting factor and focal phenotype
    
    VarS_constraint = TRUE, # Will fix variance of the sorting factor to equal the focal phenotype. The model tends to do better if you instead constrain one of the paths (e.g., e1s)
    
    # Sorting Factor
    a1s.free = T, # Freely estimate scaling factor for VA on sorting factor? (genetic homogamy, additive)
    c1s.free = F, # Freely estimate scaling factor for VC on sorting factor? (social homogamy, sibling-shared)
    t1s.free = T, # Freely estimate scaling factor for VT on sorting factor? (social homogamy, twin-shared)
    d1s.free = T, # Freely estimate scaling factor for VD on sorting factor? (genetic homogamy, dominance)
    n1s.free = F, # Freely estimate scaling factor for VN on sorting factor? (genetic homogamy, epistasis)
    e1s.free = T, # Freely estimate scaling factor for VE on sorting factor?  (idiosyncratic homogamy) (Must be fixed unless others are fixed too)
    a1s.val  = NULL, # Starting value for scaling factors. If fixed: 0 to remove from sorting factor, 1 to constrain to no scaling.
    c1s.val  = NULL, # Starting value for scaling factors. If fixed: 0 to remove from sorting factor, 1 to constrain to no scaling.
    t1s.val  = NULL, # Starting value for scaling factors. If fixed: 0 to remove from sorting factor, 1 to constrain to no scaling.
    d1s.val  = NULL, # Starting value for scaling factors. If fixed: 0 to remove from sorting factor, 1 to constrain to no scaling.
    n1s.val  = NULL, # Starting value for scaling factors. If fixed: 0 to remove from sorting factor, 1 to constrain to no scaling.
    e1s.val  = NULL, # Starting value for scaling factors. If fixed: 0 to remove from sorting factor, 1 to constrain to no scaling.
    
    # Adjust for correlated mate-preferences? for now doesn't allow for mate-preferences adjustments
    cor.preferences = F,   # Freely estimated sibling-correlated mate preferences (inflates co-in-law correlations, can bias results unless accounted for)
    preferences.mz.dz = F, # Allow preference-correlation to differ across relatedness
    
    # Gene-environment correlation in parent-generation
    parental_rGE = 0, #for now the model doesn't allow rGE
    
    # Interaction deviation or non-additive rNA
    rI.f = 0, # default (not biologically likely but conservative)
    
    # Model options
    fitFunction = "ML" #for now only option avaible
) {
  
  timestamp()
  require(OpenMx)
  selVars_pa <- selVars[1:5]         # Names of parental variables
  ntv <- length(selVars)      # number of total variables
  
  # Data checks
  
  if (!cor.preferences) {preferences.mz.dz <- F} # Can only be free if correlated preferences are included
  
  
  # Free Parameters--------------------------------------------------------------
  paths <- list(
    
    # Misc set parameters
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE, values = U.f, labels = "U1", name = "U"), # Closeness to intergenerational equilibirum (AM)
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE, values = rI.f, labels = "riI", name = "rI"), # Closeness to intergenerational equilibirum (AM)
    # ACTE (parental generation, focal phenotype)
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = a1.free, values = a1.val, label = "est_a1", name = "a1"), # Genetic influences
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = c1.free, values = c1.val, label = "est_c1", name = "c1"), # (Sibling) shared environmental influences
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = t1.free, values = t1.val, label = "est_t1", name = "t1"), # Twin shared environmental influences
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = d1.free, values = d1.val, label = "est_d1", name = "d1"), # Dominant genetic influences
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = n1.free, values = n1.val, label = "est_n1", name = "n1"), # Epistatic genetic influences
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = e1.free, values = e1.val, label = "est_e1", name = "e1") # Unique environmental influences
  )
  
  # Assortative Mating Components --------------------------------------------------------------
  # Does checks and fills NULL-objects
  if (a1s.free & c1s.free & t1s.free & e1s.free & !VarS_constraint) {stop("Model cannot be run if all paths to sorting factor are free without constraining the variance of the sorting factor.
                                                                         \nFix one of the paths (e.g., 'e1s.free = FALSE'), or select 'VarS_constraint = TRUE'")}
  if (indirect.assortment) {
    # The following code adds labels and starting values for the paths to the sorting factor.
    # If the path is not freely estimated, it either constrains it to the value given in the function arguments
    # If no argument is given, it is constrained to equal the corresponding path to the phenotype
    if(!a1s.free & is.null(a1s.val)) {lab_a1s <- "est_a1"; a1s.free <- a1.free; a1s.val = a1.val} else {lab_a1s <- "est_a1s"; if (is.null(a1s.val)) {a1s.val = a1.val} } 
    if(!c1s.free & is.null(c1s.val)) {lab_c1s <- "est_c1"; c1s.free <- c1.free; c1s.val = c1.val} else {lab_c1s <- "est_c1s"; if (is.null(c1s.val)) {c1s.val = c1.val} } 
    if(!t1s.free & is.null(t1s.val)) {lab_t1s <- "est_t1"; t1s.free <- t1.free; t1s.val = t1.val} else {lab_t1s <- "est_t1s"; if (is.null(t1s.val)) {t1s.val = t1.val} } 
    if(!d1s.free & is.null(d1s.val)) {lab_d1s <- "est_d1"; d1s.free <- d1.free; d1s.val = d1.val} else {lab_d1s <- "est_d1s"; if (is.null(d1s.val)) {d1s.val = d1.val} } 
    if(!n1s.free & is.null(n1s.val)) {lab_n1s <- "est_n1"; n1s.free <- n1.free; n1s.val = n1.val} else {lab_n1s <- "est_n1s"; if (is.null(n1s.val)) {n1s.val = n1.val} } 
    if(!e1s.free & is.null(e1s.val)) {lab_e1s <- "est_e1"; e1s.free <- e1.free; e1s.val = e1.val} else {lab_e1s <- "est_e1s"; if (is.null(e1s.val)) {e1s.val = e1.val} } 
    
  } else {
  # If no indirect assortment, then all paths so sorting factor must equal the corresponding path to the focal phentotype  
    lab_a1s <- "est_a1"; a1s.free <- a1.free; a1s.val = a1.val
    lab_c1s <- "est_c1"; c1s.free <- c1.free; c1s.val = c1.val
    lab_t1s <- "est_t1"; t1s.free <- t1.free; t1s.val = t1.val
    lab_d1s <- "est_d1"; d1s.free <- d1.free; d1s.val = d1.val
    lab_n1s <- "est_n1"; n1s.free <- n1.free; n1s.val = n1.val
    lab_e1s <- "est_e1"; e1s.free <- e1.free; e1s.val = e1.val
    # VarS_constraint == FALSE
  }

  paths <- append(paths, list(
    
    # ACTE (parental generation, sorting factor)
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = a1s.free, values = a1s.val, label = lab_a1s, name = "a1s"),
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = c1s.free, values = c1s.val, label = lab_c1s, name = "c1s"),
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = t1s.free, values = t1s.val, label = lab_t1s, name = "t1s"),
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = d1s.free, values = d1s.val, label = lab_d1s, name = "d1s"),
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = n1s.free, values = n1s.val, label = lab_n1s, name = "n1s"),
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = e1s.free, values = e1s.val, label = lab_e1s, name = "e1s"),
    
    # Variance of sorting factor, and covariance between focal phenotype and sorting factor
    mxAlgebra(expression = a1s^2 + c1s^2 + e1s^2 + 2 * a1s * w * c1s + t1s^2 + d1s^2 + n1s^2, name = "VarS"),
    mxAlgebra(expression = a1 * a1s + c1 * c1s + e1 * e1s + c1 * w * a1s + a1 * w * c1s + t1 * t1s + d1 * d1s + n1 * n1s, name = "covPS"),
    mxAlgebra(expression = mu * (a1s + c1s * w)^2, name = "partner_rG"),
    
    # Estimates partner covariance on sorting factor, and rescales to copath
    mxMatrix(type = "Lower", nrow = 1, ncol = 1, free = mu.free, values = mu.val, labels = "mu1", name = "mu_param"),
    mxAlgebra(expression = mu_param / VarS^2, name = "mu"),
    
    # Gene-environment covariance
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE, values = parental_rGE, label = "w_est", name = "w")
  ))
  
  # Adds constraint if desired
  if (VarS_constraint) {VarS_cons <- mxConstraint(expression = VarS == var1, name = "VarS_constraint")} else {VarS_cons <- NULL}
  
  # Variance components--------------------------------------------------------------
  variances <- list(
    
    # Parental Generation (focal phenotype)
    mxAlgebra(expression = a1^2, name = "VA1"),
    mxAlgebra(expression = c1^2, name = "VC1"),
    mxAlgebra(expression = t1^2, name = "VT1"),
    mxAlgebra(expression = d1^2, name = "VD1"),
    mxAlgebra(expression = n1^2, name = "VN1"),
    mxAlgebra(expression = 2 * c1 * w * a1, name = "VrGE1"),
    mxAlgebra(expression = e1^2, name = "VE1"),
    
    # Parental Generation (sorting factor)
    mxAlgebra(expression = a1s^2, name = "VA1_s"),
    mxAlgebra(expression = c1s^2, name = "VC1_s"),
    mxAlgebra(expression = t1s^2, name = "VT1_s"),
    mxAlgebra(expression = d1s^2, name = "VD1_s"),
    mxAlgebra(expression = n1s^2, name = "VN1_s"),
    mxAlgebra(expression = 2 * c1s * w * a1s, name = "VrGE1_s"),
    mxAlgebra(expression = e1s^2, name = "VE1_s")
  )
  
  # Derived parameters--------------------------------------------------------------
  
  
  # Genetic correlation between siblings
  FS_rA <- mxAlgebra(expression = (1 + U * partner_rG) / 2, name = "FS_rA")
  FS_rD <- mxAlgebra(expression = 1 / 4, name = "FS_rD")
  FS_rN <- mxAlgebra(expression = rI, name = "FS_rN")
 
  # Algebra for expected variances--------------------------------------------------------------
  var1     <- mxAlgebra(expression= VA1 + VC1 + VrGE1 + VE1 + VT1 + VD1 + VN1, name="var1")
  
  # Algebra for expected covariances--------------------------------------------------------------
  # Between twins (or siblings or however the related pair in the parental generation is related)
  cov_MZ     <- mxAlgebra(expression =     1*VA1 + VC1 + VrGE1 +     1*VD1 +     1*VN1 + VT1 , name = "cov_MZ")  # Covariation between siblings
  cov_DZ     <- mxAlgebra(expression = FS_rA*VA1 + VC1 + VrGE1 + FS_rD*VD1 + FS_rN*VN1 + VT1 , name = "cov_DZ")  # Covariation between siblings
  cov_FS     <- mxAlgebra(expression = FS_rA*VA1 + VC1 + VrGE1 + FS_rD*VD1 + FS_rN*VN1       , name = "cov_FS")  # Covariation between siblings

  cov_In_MZ  <- mxAlgebra(expression = covPS*mu* (   1*a1*a1s  + c1*c1s + c1s*w*a1 + c1*w*a1s +     1*d1*d1s +     1*n1*n1s + t1*t1s)  , name = "cov_In_MZ") # Covariation between sibling-inlaw
  cov_In_DZ  <- mxAlgebra(expression = covPS*mu* (FS_rA*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + FS_rD*d1*d1s + FS_rN*n1*n1s + t1*t1s)  , name = "cov_In_DZ") # Covariation between sibling-inlaw
  cov_In_FS  <- mxAlgebra(expression = covPS*mu* (FS_rA*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + FS_rD*d1*d1s + FS_rN*n1*n1s         )  , name = "cov_In_FS") # Covariation between sibling-inlaw

  cov_CIn_MZ  <- mxAlgebra(expression = covPS^2 * (mu^2 *  (   1*a1s^2  + c1s^2 + 2*c1s*w*a1s + t1s^2 +     1*d1s^2 +     1*n1s^2)),  name = "cov_CIn_MZ") # Covariation between co-inlaws
  cov_CIn_DZ  <- mxAlgebra(expression = covPS^2 * (mu^2 *  (FS_rA*a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2 + FS_rD*d1s^2 + FS_rN*n1s^2)), name = "cov_CIn_DZ") # Covariation between co-inlaws

  cov_Mate <- mxAlgebra(expression = covPS*mu*covPS, name="cov_Mate") # Covariance between partners
  
  # Implied covariance matrices--------------------------------------------------------------
  expCov_MZ <- mxAlgebra(expression= rbind(cbind(var1,       cov_Mate,  cov_In_FS, cov_In_MZ, cov_CIn_MZ), 
                                           cbind(cov_Mate,   var1,      cov_FS,    cov_MZ,    cov_In_MZ), 
                                           cbind(cov_In_FS,  cov_FS,    var1,      cov_FS,    cov_In_FS), 
                                           cbind(cov_In_MZ,  cov_MZ,    cov_FS,    var1,      cov_Mate), 
                                           cbind(cov_CIn_MZ, cov_In_MZ, cov_In_FS, cov_Mate,  var1)), 
                         name="expCov_MZ")
  
  expCov_DZ <- mxAlgebra(expression= rbind(cbind(var1,        cov_Mate,  cov_In_FS, cov_In_DZ, cov_CIn_DZ), 
                                            cbind(cov_Mate,   var1,      cov_FS,    cov_DZ,    cov_In_DZ), 
                                            cbind(cov_In_FS,  cov_FS,    var1,      cov_FS,    cov_In_FS), 
                                            cbind(cov_In_DZ,  cov_DZ,    cov_FS,    var1,      cov_Mate), 
                                            cbind(cov_CIn_DZ, cov_In_DZ, cov_In_FS, cov_Mate,  var1)), 
                         name="expCov_DZ")
  
  # Create Algebra for expected Mean Matrix
  meanG <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = T, values = 0, labels = "mean", name = "meanG")
  exp_MZ <- mxExpectationNormal(covariance = "expCov_MZ", means = "meanG", dimnames = selVars)
  exp_DZ <- mxExpectationNormal(covariance = "expCov_DZ", means = "meanG", dimnames = selVars)
  
  # Model objetcs--------------------------------------------------------------
  # Create Data Objects
  groups <- c("MZ", "DZ")
  data_MZ <- mxData(observed = mzData, type = "raw")
  data_DZ <- mxData(observed = dzData, type = "raw")
  
  # Compile objects into lists (it eases handling)
  common_objects <- list(meanG, var1, cov_Mate, FS_rA, FS_rD, FS_rN)
  exp_cov_list <- list(expCov_MZ, expCov_DZ)
  
  # Group algebras
  algebras_MZ <- list(cov_MZ, cov_FS, cov_In_FS, cov_In_MZ, cov_CIn_MZ)
  algebras_DZ <- list(cov_DZ, cov_FS, cov_In_FS, cov_In_DZ, cov_CIn_DZ)
  algebras_all <- list(algebras_MZ, algebras_DZ)
  
  # function
  fitFun <- mxFitFunctionML()
  multi <- mxFitFunctionMultigroup(groups)
  
  # Group models
  model_MZ <- mxModel(paths, variances, common_objects, algebras_MZ, expCov_MZ, exp_MZ, data_MZ, fitFun, name = "MZ")
  model_DZ <- mxModel(paths, variances, common_objects, algebras_DZ, expCov_DZ, exp_DZ, data_DZ, fitFun, name = "DZ")
  model_list <- list(model_MZ, model_DZ)
  
  # CI objects--------------------------------------------------------------
  # Create Confidence Interval Objects
  est_var_Col <<- c(
    "varP", "varS",
    "varP_A", "varP_C", "varP_T", "varP_D", "varP_N", "varP_E",
    "varS_A", "varS_C", "varS_T", "varS_D", "varS_N", "varS_E",
    "Partner Correlation (P)", "Partner Correlation (S)", "Partner rG"
  )
  
  est_var <- mxAlgebra(
    expression = cbind(
      var1, VarS,
      VA1 / var1, VC1 / var1, VT1 / var1, VD1 / var1, VN1 / var1, VE1 / var1,
      VA1_s / VarS, VC1_s / VarS, VT1_s / VarS, VD1_s / VarS, VN1_s / VarS, VE1_s / VarS,
      cov_Mate / var1, mu * VarS, partner_rG
    ),
    name = "est_var",
    dimnames = list("est", est_var_Col)
  )
  
  est_path_Col <<- c(
    "mean",
    "mu",
    "a1s","c1s", "t1s", "d1s", "n1s", "e1s",
    "a1", "c1","t1","d1","n1","e1"
  )
  
  est_path <- mxAlgebra(
    expression = cbind(
      mean,
      mu_param,
      a1s, c1s, t1s, d1s, n1s, e1s,
      a1, c1, t1, d1, n1, e1
    ),
    name = "est_path",
    dimnames = list("est", est_path_Col)
  )
  
  
  CI <- mxCI(c("est_var", "est_path"))
  CI_names <<- c(est_var_Col)
  results_objects <- list(
    est_path,
    est_var
  )

    
    
  # Model build--------------------------------------------------------------
  model_iamCOTS <- mxModel(paths, multi, VarS_cons, common_objects, variances, algebras_all, model_list,exp_cov_list,CI,results_objects, name = name)
  
} 




