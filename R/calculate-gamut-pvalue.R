


###############################################
## A Function to Calculate the GAMuT P-Value ##
###############################################

## Details may be found here:
## https://epstein-software.github.io/GAMuT/GAMuT-example.html

gamut_pvalue <- function(Phenotype_Matrix_Centered,
                         Genotype_Matrix){
  
  ## Load packages
  library(CompQuadForm)
  #library(devtools)
  
  
  ##############################
  ## Load the GAMuT Functions ##
  ##############################
  
  ## If devtools worked...
  #devtools::source_url(url = "https://raw.githubusercontent.com/epstein-software/GAMuT/master/GAMuT-functions.R")
  
  ## Manually load the functions from the Github page (last accessed on 2020-02-27):
  
  ## GAMuT-functions.R
##
## 2016-March-13
##
## this script contains functions used in the main program
## for GAMuT analysis
##
##----------------------------------------------------------------------
## descriptions of individual functions:
##----------------------------------------------------------------------
##
## * proj_GAMuT_pheno
##   constructs projection matrix and corresponding eigenvalues for phenotypes
##
## * linear_GAMuT_pheno
##   constructs linear kernel matrix and corresponding eigenvalues for phenotypes
##
## * linear_GAMuT_geno
##   constructs linear kernel matrix and corresponding eigenvalues for genotypes 
##
## * TestGAMuT
##   constructs GAMuT statistic and returns p-value


##----------------------------------------------------------------------
## phenotypic similarity:  projection matrix
##----------------------------------------------------------------------
proj_GAMuT_pheno <- function(X){
    out = vector("list", 2)
    names(out) = c("Kc", "ev_Kc")

    out$Kc = X %*% solve(t(X) %*% X) %*% t(X)   # projection matrix
    out$ev_Kc = rep(1, ncol(X))                 # find eigenvalues
    return(out)	
}

##----------------------------------------------------------------------
## phenotypic similarity:  linear kernel
##----------------------------------------------------------------------
linear_GAMuT_pheno <- function(X){
    out = vector("list", 2)
    names(out) = c("Kc", "ev_Kc")
    
    Kc = t(X) %*% X   # transposed kernel to find eigenvalues
    ev_Kc = eigen(Kc, symmetric=T, only.values=T)$values  
    
    out$Kc = X %*% t(X) # similiarity kernel for test
    out$ev_Kc = ev_Kc[ev_Kc > 1e-08]
    return(out)
}

##----------------------------------------------------------------------
## genotypic similarity:  linear kernel
##----------------------------------------------------------------------
linear_GAMuT_geno <- function(X){
    out = vector("list", 2)
    names(out) = c("Lc", "ev_Lc")
    
    Lc = t(X) %*% X   # transposed kernel to find eigenvalues
    ev_Lc = eigen(Lc, symmetric=T, only.values=T)$values  
    
    out$Lc = X %*% t(X) # similiarity kernel for test
    out$ev_Lc = ev_Lc[ev_Lc > 1e-08]
    return(out)	 
}


##----------------------------------------------------------------------
## constructing the GAMuT statistic and deriving p-value:
##----------------------------------------------------------------------
TestGAMuT <- function(Yc, lambda_Y, Xc, lambda_X) {

    ## test statistic:
    m = nrow(Yc) # number of subjects in study
    GAMuT = (1/m) * sum(sum(t(Yc) * Xc))  
    

    ## populate vector of all pairwise combination of eigenvalues
    ## from the phenotype and genotype similarity matrices:
    Z <- (as.matrix(lambda_Y)) %*% t(as.matrix(lambda_X))
    Zsort <- sort(Z, decreasing=T)

    ## derive p-value of GAMuT statistic:
    scoredavies = GAMuT*m^2
    results_score <- davies(scoredavies, Zsort)
    davies_pvalue <- (results_score$Qq)

    return(davies_pvalue)
} 
  
  
  ######################
  ## Run the Analysis ##
  ######################
  
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


