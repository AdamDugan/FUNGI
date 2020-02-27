


###################
## Load Packages ##
###################

library(gam)
library(fda)
library(nlme)
library(mgcv)

#pkgs = c("gam","fda","nlme")
#for(p in pkgs) if(p %in% rownames(installed.packages()) == FALSE) {install.packages(p, dependencies = TRUE)}
#for(p in pkgs) suppressPackageStartupMessages(library(p, quietly=TRUE, character.only=TRUE))
#rm(p, pkgs)




####################################
## Map the SNP Locations to [0,1] ##
####################################

Map_To_01 = function(SNP_Positions){
  return( (SNP_Positions - min(SNP_Positions)) / (SNP_Positions[ length(SNP_Positions) ] - min(SNP_Positions) ))
}




###########################################
## Remove Participants With Missing Data ##
###########################################

Remove_Missing_Cases = function(Genotype_Matrix = geno,
                                Phenotype_DataFrame = phe,
                                Protein_Variable_Names = c(),
                                Phenotype_ID_Variable_Name = "IID"){
  
  ## Input for troubleshooting
  # Genotype_Matrix = APO_Geno
  # Phenotype_DataFrame = adni
  # Protein_Variable_Names = c("ABETA","TAU")
  # Phenotype_ID_Variable_Name = "RID"
  # rm(APO_Geno, CLU_Geno, adni)
  
  ## Check the inputs
  if( !is.data.frame(Phenotype_DataFrame) ){ return( cat(Phenotype_DataFrame,"is not a data.frame!") ) }
  if( sum( Protein_Variable_Names %in% names(Phenotype_DataFrame) ) < length(Protein_Variable_Names) ){ return( cat(Protein_Variable_Names,"are not found in the Phenotype_DataFrame!") ) }
  if( !Phenotype_ID_Variable_Name %in% names(Phenotype_DataFrame) ){ return( cat(Phenotype_ID_Variable_Name,"is not found in the Phenotype_DataFrame!") ) }
  
  ## Identify inputs
  # geno = Genotype_Matrix
  # phe = Phenotype_DataFrame
  # proteins = phe[,Protein_Variable_Names]
  # rm(Genotype_Matrix, Phenotype_DataFrame, Protein_Variable_Names)
  
  ## How many cases are missing the protein?
  N_Missing_Proteins = unlist( lapply(X = Protein_Variable_Names, FUN = function(x){ return( sum(is.na(Phenotype_DataFrame[,x])) ) } ) )
  
  ## Identify which rows are missing at least 1 protein
  Missing_Any_Protein = apply(X = cbind(Phenotype_DataFrame[,Protein_Variable_Names]), MARGIN = 1,
                              FUN = function(x){ return( sum(is.na(x)) ) } )
  
  ## Remove cases missing the protein
  Phenotype_DataFrame = Phenotype_DataFrame[ Missing_Any_Protein < 1, ]
  
  ## Identify IDs from the genotype data
  IDs_Genotype = names( data.frame(Genotype_Matrix, check.names = FALSE) )
  
  ## Remove phenotype cases without genotype data
  Phenotype_DataFrame = Phenotype_DataFrame[ Phenotype_DataFrame[, Phenotype_ID_Variable_Name] %in% IDs_Genotype, ]
  
  ## Identify the IDs of the remaining subjects
  IDs_Phenotype = as.character( Phenotype_DataFrame[,Phenotype_ID_Variable_Name] )
  
  ## Identify which genotype columns to keep
  IDs_Genotype_ToKeep = IDs_Genotype[ IDs_Genotype %in% IDs_Phenotype ]
  
  ## Remove genotype data for cases missing phenotype data
  Genotype_Matrix = Genotype_Matrix[, IDs_Genotype_ToKeep ]
  
  ## Order the phenotype data by the genotype IDs
  Phenotype_DataFrame = Phenotype_DataFrame[ match(IDs_Genotype_ToKeep , Phenotype_DataFrame[,Phenotype_ID_Variable_Name] ), ]
  
  ## Do the number of phenotype rows match the number of genotype columns?
  if( !(nrow(Phenotype_DataFrame) == ncol(Genotype_Matrix)) ){ return( "Number of participants in the pheotype data does not match the number in the genotype data!" ) }
  else{
    ## Return a list containing the filtered matrices
    return( list(Proteins = Protein_Variable_Names,
                 N_Missing_Proteins = N_Missing_Proteins,
                 Genotype_Matrix = Genotype_Matrix,
                 Phenotype_DataFrame = Phenotype_DataFrame) )
  }
}





########################################
## Smooth the Genotype Data Using gam ##
########################################

Smooth_Genotypes_gam = function(Genotype_Matrix = geno,
                                SNP_Positions = SNP_Position_01,
                                Subject_Per_Column = TRUE){
  
  library(mgcv)
  library(gam)
  require(mgcv)
  
  ## Inputs for troubleshooting
  # Genotype_Matrix = geno_flipped
  # SNP_Positions = SNP_Position_01
  # Subject_Per_Column = TRUE
  # rm(geno_flipped, SNP_Position_01)
  
  
  ## Rename the inputs
  geno = Genotype_Matrix
  pos = SNP_Positions
  miss.subj = c()
  # rm(Genotype_Matrix, SNP_Positions)
  
  ## If the subjects are organized by row, then transpose the Genotype_Matrx
  if( Subject_Per_Column == FALSE ){ geno = t(geno) }
  
  ## Create an empty object to save the smoothed values
  result = matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
  
  ## Smooth the data for each participant using gam
  for(i in 1:ncol(geno)){
    
    ## See if any errors arise while fitting the model
    mod = tryCatch( gam(geno[,i] ~ s(pos, fx = FALSE, k=-1, bs='cr') ),
                    error = function(e) e,
                    warning = function(w) w)
    
    ## Check if there were too many missing values to fit the model
    if( is(mod, "error") ){
      cat("Participant", i, "has too many NA values to do anything meaningful", '\n')
      miss.subj = c(miss.subj, i) }
    else{
      
      ## Fit the model for a single participant
      mod = gam(geno[,i] ~ s(pos, fx = FALSE, k = -1, bs='cr') )
      
      ## Calcualte the fitted values for all SNPS (whether or not they were missing for this subject)
      result[,i] = predict.gam(mod, newdata = data.frame(pos = pos))
    }
  }
  
  ## Transpose the data if the genotype data were originally organized by row
  if( Subject_Per_Column == FALSE ){ return( list(Smoothed_Genotype = t(result), Too_Many_Missing_SNPs = miss.subj ) ) }
  else{ return( list(Smoothed_Genotype = result, Too_Many_Missing_SNPs = miss.subj )) }
}





#####################################################################################
## Calculate the Squared Eigen Values for the Empirical Variance-Covariance Matrix ##
#####################################################################################

## Inputs for troubleshooting:
# Knots = variant_positions_01
# SNP_Positons = variant_positions_01
# Genotype_Matrix = geno_flipped_smooth

Eigen_Values_EmpVarCov_fda = function(Knots = SNP_Position_01,
                                      SNP_Positons = SNP_Position_01,
                                      Genotype_Matrix = geno){
  
  library(fda)
  require(fda)
  
  ## Rename inputs
  pos = SNP_Positons
  geno = Genotype_Matrix
  
  ## Define values
  norder = 4
  nbasis = (length(Knots) + norder - 2)
  Lfdobj = 2
  
  ## Calculate the empirical variance-covariance matrix
  ybasis = create.bspline.basis(range(Knots), nbasis, norder, Knots)
  yfdPar = fdPar(ybasis, Lfdobj, lambda = 1e-13)
  yfd = smooth.basis(pos, geno, yfdPar)$fd
  cvar.fd = var.fd(yfd)
  E.fd = eval.bifd(pos, pos, cvar.fd)
  
  ## Calculate the eigen values for the empirical variance-covariance matrix
  eigen_val = eigen(x = E.fd)$values
  
  ## Return the squared eigen values
  return( eigen_val )
}
#rm(Knots, SNP_Positons, Genotype_Matrix)





###############################################
## Calculate the Adjusted Degrees of Freedom ##
###############################################

Adjusted_Degrees_of_Freedom = function(Number_of_Parameters_Full,
                                       Number_of_Parameters_Reduced,
                                       Eigen_Values,
                                       Number_of_Subjects){
  
  ## Rename inputs
  P = Number_of_Parameters_Full
  q = Number_of_Parameters_Reduced
  eigen_val = Eigen_Values
  eigen_val.sq = Eigen_Values^2
  N_Subjects = Number_of_Subjects
  
  ## Calcualte the adjusted degrees of freedom for the F distribution
  df_1 = (P-q)*((N_Subjects-1)*(sum(eigen_val)^2-2*sum(eigen_val.sq)/(N_Subjects-1)))/((N_Subjects-2)*(sum(eigen_val.sq)-sum(eigen_val)^2/(N_Subjects-2)))
  df_2 = (N_Subjects-P)*((N_Subjects-1)*(sum(eigen_val)^2-2*sum(eigen_val.sq)/(N_Subjects-1)))/((N_Subjects-2)*(sum(eigen_val.sq)-sum(eigen_val)^2/(N_Subjects-2)))
  
  ## Return the adjusted degrees of freedom
  return( c(df_1, df_2))
}





