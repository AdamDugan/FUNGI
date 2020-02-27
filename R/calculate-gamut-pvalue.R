


###############################################
## A Function to Calculate the GAMuT P-Value ##
###############################################

## Details may be found here:
## https://epstein-software.github.io/GAMuT/GAMuT-example.html

gamut_pvalue <- function(Phenotype_Matrix_Centered,
                         Genotype_Matrix){
  
  ## Load packages
  library(CompQuadForm)
  library(devtools)
  
  ## Load the GAMuT functions
  devtools::source_url(url = "https://raw.githubusercontent.com/epstein-software/GAMuT/master/GAMuT-functions.R")
  
  ## Form the phenotypic similarity matrix and corresponding eigenvalues
  proj_pheno <- proj_GAMuT_pheno(X = Phenotype_Matrix_Centered)
  Yc <- proj_pheno$Kc                # projection matrix
  lambda_Y <- proj_pheno$ev_Kc 
  
  ## Form the genotypic similarity matrix and corresponding eigenvalues
  MAF <- colMeans(Genotype_Matrix)/2                      # sample MAF of each variant in the sample
  beta_weight <- dbeta(MAF, 1, 25) / dbeta(0, 1, 25) # assume beta-distribution weights
  
  ## Then the weighted rare variants and the centered genotype matrix are calculated by
  G0 <- as.matrix(Genotype_Matrix) %*% diag(beta_weight) # Weighted rare variants
  G  <- as.matrix(scale(G0, center = TRUE, scale = FALSE))  # centered genotype matrix
  
  ## Next we call the function linear_GAMuT_geno() to construct the weighted linear kernel matrix for genotypes along with the eigenvalues:
  linear_geno <- linear_GAMuT_geno(X = G)
  Xc <- linear_geno$Lc                        # linear kernel similarity matrix
  lambda_X <- linear_geno$ev_Lc               # eigenvalues of Xc
  
  ## Construct the GAMuT and obtain the p-value
  GAMuT_pvalue <- TestGAMuT(Yc = Yc,
                            lambda_Y = lambda_Y,
                            Xc = Xc,
                            lambda_X = lambda_X)
  
  ## Return the p-value
  return(GAMuT_pvalue)
  
}


