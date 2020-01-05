


#############################################
## Items to Help Troubleshoot the Function ##
#############################################

# ## Load the 1000 Genomes data
# ped <- read.table(file = "olgas-files/ped.ped",
#                   sep = "\t")
# info <- read.table(file = "olgas-files/info.info",
#                    header = TRUE)
# 
# ## Randomly select a starting position for the genetic region
# seq_start <- round( runif(n = 1, min = 0, max = (ncol(ped) - 500)) )
# 
# ## Calculate the ending point for the genetic region
# seq_end <- (seq_start + 500 - 1)
# 
# ## Extract the indexes for the selected SNPs
# variants_all = seq_start:seq_end
# 
# ## Randomly select individuals to include
# geno <- ped[ sample(x = 1:nrow(ped),
#                     size = 300,
#                     replace = TRUE), ]
# 
# ## Subset to the genetic region
# geno <- geno[, variants_all ]
# 
# ## Define some inputs
# genotype_matrix <- as.matrix(geno)
# variants_maf <- info$maf[ variants_all ]
# variants_index <- variants_all
# proportions_true <- c(0.10, 0.10, 0.10)
# effects_range_true <- data.frame(min = c(0.005, 0.005, 0.005),
#                                  max = c(0.025, 0.025, 0.025))
# effects_range_null <- c(0.000, 0.005)
# corr_matrix <- matrix(c(1, 0, 0,
#                         0, 1, 0,
#                         0, 0, 1),
#                       nrow = 3, ncol = 3, byrow = TRUE)
# dichotomize_probabilistically <- TRUE
# 
# ## Remove R objects
# rm(seq_start, seq_end, variants_all, geno)





########################################
## A Function to Simulate Binary Data ##
########################################

simulate_binary_phenotypes <- function(genotype_matrix,
                                       variants_maf,
                                       variants_index,
                                       proportions_true = c(0.10, 0.10, 0.10),
                                       effects_range_true = data.frame(min = c(0.005, 0.005, 0.005),
                                                                       max = c(0.025, 0.025, 0.025)),
                                       effects_range_null = c(0.000, 0.005),
                                       corr_matrix = matrix(c(1, 0, 0,
                                                               0, 1, 0,
                                                               0, 0, 1),
                                                             nrow = 3, ncol = 3, byrow = TRUE),
                                       dichotomize_probabilistically = TRUE){
  
  library(MASS)
  require(MASS)
  
  ########################
  ## Define Some Values ##
  ########################
  
  ## The number of variants
  n_variants <- ncol(genotype_matrix)
  
  ## The simulation details
  simulation_details <- data.frame(Variant_MAF = variants_maf,
                                   Variant_Index = variants_index)
  
  
  
  
  
  ####################################
  ## Simulate Continuous Phenotypes ##
  ####################################
  
  ## How many variants are truly associated?
  n_variants_true <- round(n_variants*proportions_true)
  
  ## Randomly select the truly associated variants
  variants_true <- apply(X = matrix(n_variants_true),
                         MARGIN = 1,
                         FUN = function(x){
                           tmp <- sample(x = 1:n_variants, size = x)
                           return( tmp[ order(tmp) ] ) } )
  
  ## Identify the non-causal variants
  variants_null <- apply(X = variants_true,
                         MARGIN = 2,
                         FUN = function(x){
                           tmp <- 1:n_variants
                           return( tmp[ !tmp %in% x ] ) } )
  
  ## Save the details
  simulation_details$Variants_Binary <- ""
  for(i in 1:nrow(simulation_details)){
    for(k in 1:ncol(variants_true)){
      if( i %in% variants_true[,k] ){ simulation_details$Variants_Binary[i] <- paste0(simulation_details$Variants_Binary[i], "Causal", k) }
    }
  }
  rm(i, k)
  
  ## Generate the causal effects
  effects_causal <- matrix(NA, ncol = length(proportions_true), nrow = max(n_variants_true))
  for(i in 1:ncol(effects_causal)){
    effects_causal[, i] <- runif(n = n_variants_true[i],
                                 min = effects_range_true$min[i],
                                 max = effects_range_true$max[i])
  }
  rm(i)
  
  ## Generate the null effects
  effects_null <- matrix(NA, ncol = length(proportions_true), nrow = (n_variants - min(n_variants_true)) )
  for(i in 1:ncol(effects_null)){
    effects_null[, i] <- runif(n = (n_variants - n_variants_true[i]),
                               min = effects_range_null[1],
                               max = effects_range_null[2])
  }
  rm(i)
  
  ## Combine the effects into a single matrix
  effects <- matrix(NA, nrow = n_variants, ncol = length(proportions_true))
  for(i in 1:ncol(effects)){
    effects[ variants_true[, i], i] <- effects_causal[, i]
    effects[ variants_null[, i], i] <- effects_null[, i]
  }
  rm(i, effects_causal, effects_null)
  
  ## Transform the proportions of variance explained
  effects_transformed <- matrix(NA, ncol = ncol(effects), nrow = nrow(effects))
  for(i in 1:nrow(effects)){
    for(k in 1:ncol(effects)){
      
      ## Identify the major and minor allele frequencies
      q <- variants_maf[i]
      p <- (1 - q)
      
      ## Calculate the transformed value
      effects_transformed[i, k] <- sqrt( ( effects[i, k] / (1 - effects[i, k] ) ) / (2*p*q) )
      
      ## Randomly make some effects negative
      effects_transformed[i, k] <- (effects_transformed[i, k] * sample(x = c(-1, 1), size = 1))
      
    }
  }
  rm(i, k, p, q)
  
  ## Generate some correlated errors
  errors <- MASS::mvrnorm(n = nrow(genotype_matrix),
                          mu = rep(0, length(proportions_true)),
                          Sigma = corr_matrix)
  
  ## Make the geno object a matrix
  geno_matrix <- as.matrix(genotype_matrix)
  
  ## Genreate the simulated data
  X <- (geno_matrix %*% effects_transformed) + errors
  
  ## Remove R objects
  rm(errors, geno_matrix)
  
  ## Combine the simualted data into a data.frame
  dat = data.frame(X_Cont_1 = X[, 1],
                   stringsAsFactors = FALSE)
  for(i in 2:ncol(X)){
    dat[, paste0("X_Cont_", i) ] <- X[, i]
  }
  rm(i)
  
  ## Remove row.names
  row.names(dat) = NULL
  
  ## Make the continuous phenotypes binary
  if( dichotomize_probabilistically ){
    for(k in 1:ncol(dat)){
      dat[, paste0("X", k) ] <- as.numeric( dat[, k] > runif(n = nrow(dat), min = min(dat[, k]), max = max(dat[, k])) )
    }
    rm(k)
  }
  else{
    for(k in 1:ncol(dat)){
      tmp <- runif(n = 1, min = min(dat[, k]), max = max(dat[, k]))
      dat[, paste0("X", k) ] <- as.numeric( dat[, k] > tmp )
    }
    rm(k, tmp)
  }

  ## Boxplots of the data
  # boxplot(dat$X_Cont_1 ~ dat$X1)
  # boxplot(dat$X_Cont_2 ~ dat$X2)
  # boxplot(dat$X_Cont_3 ~ dat$X3)
  
  ## Return the simulated data
  return(dat)
  
}