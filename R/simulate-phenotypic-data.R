


#####################################
## A Function to Simulate the Data ##
#####################################

simulate_phenotypes <- function(genotype_matrix,
                                variants_maf,
                                variants_index,
                                binary_proportion_true = 0.10,
                                binary_effect_range_true,
                                binary_effect_range_null,
                                binary_effect_on_continuous_range = c(0.10, 0.25),
                                continuous_correlation_matrix = matrix(c(1, 0.3, 0.3, 1),
                                                                       ncol = 2,
                                                                       byrow = TRUE),
                                continuous_propotions_true = c(0.05, 0.05),
                                continuous_effects_range_true = c(0.005, 0.015),
                                continuous_effects_range_null = c(0, 0.0025)){
  
  library(MASS)
  require(MASS)
  
  ########################
  ## Define Some Values ##
  ########################
  
  ## The number of variants
  n_variants <- ncol(genotype_matrix)
  
  ## Save the simualted SNP information
  # simulation_details <- data.frame(Variant = info$SNP[ variants_all ],
  #                                  Variant_Position = info$posi[ variants_all ],
  #                                  Variant_MAF = info$maf[ variants_all ],
  #                                  Variant_Index = variants_all)
  simulation_details <- data.frame(Variant_MAF = variants_maf,
                                   Variant_Index = variants_index)
  
  
  
  
  
  ########################################################
  ## Simulate a Continuous Phenotype and Dichotomize It ##
  ########################################################
  
  ## The proportion of variants that are truly associated with the binary phenotype
  proportion_true <- binary_proportion_true
  
  ## The range of effect sizes for these truly associated variants
  effect_true <- binary_effect_range_true
  
  ## The range of effect sizes for the other variants
  effect_null <- binary_effect_range_null
  
  ## How many variants are truly associated?
  n_variants_true <- round(n_variants*proportion_true)
  
  ## Randomly select the truly associated variants
  variants_true <- sample(x = 1:n_variants, size = n_variants_true)
  variants_true <- variants_true[ order(variants_true) ]
  
  ## Identify the other unrelated variants
  variants_null <- 1:n_variants
  variants_null <- variants_null[ !variants_null %in% variants_true ]
  
  ## Save the details
  simulation_details$Variants_Binary <- "Null"
  simulation_details$Variants_Binary[ variants_true ] <- "True"
  
  ## Randomly generate the effects
  tmp_true <- runif(n = n_variants_true, min = effect_true[1], max = effect_true[2])
  tmp_null <- runif(n = (n_variants - n_variants_true), min = effect_null[1], max = effect_null[2])
  
  ## Combine the effects into a single vector
  effects <- rep(NA, n_variants)
  effects[ variants_true ] <- tmp_true
  effects[ variants_null ] <- tmp_null
  
  ## Remove R objects
  rm(tmp_true, tmp_null)
  
  ## Save the results
  simulation_details$Effects_Binary <- effects
  
  ## Transform the proportions of variance explained
  effects_transformed <- NA
  for(i in 1:ncol(genotype_matrix)){
    
    ## Identify the major and minor allele frequencies
    q <- variants_maf[i]
    p <- (1 - q)
    
    ## Calculate the transformed value
    effects_transformed[i] <- sqrt( ( effects[i] / (1 - effects[i] ) ) / (2*p*q) )
  }
  rm(i, p, q)
  
  ## Randomly make some effects negative
  effects_transformed <- (effects_transformed * sample(x = c(-1, 1), size = length(effects_transformed), replace = TRUE))
  
  ## Generate independent errors
  tmp_errors <- rnorm(n = nrow(genotype_matrix))
  
  ## Simulate a continuous outcome
  z <- as.vector( t(matrix(effects_transformed)) %*% t(genotype_matrix) + tmp_errors )
  
  ## Remove R objects
  rm(tmp_errors)
  
  # ## Save the transformed effects
  # simulation_details$Effects_Transformed_Binary <- effects_transformed
  
  ## Make the continuous phenotype binary
  z_binary <- as.numeric( z > runif(n = nrow(geno), min = min(z), max = max(z)) )
  
  
  
  
  
  #############################################
  ## Create Correlated Continuous Phenotypes ##
  #############################################
  
  ## The number of correlated phenotypes
  n_corr_phen <- 2
  
  ## The proportion of variants that are truly associated with the binary phenotype
  proportion_true <- continuous_propotions_true
  
  ## The correlation matrix for these phenotypes
  corr_matrix <- continuous_correlation_matrix
  
  ## The range of effect sizes for these truly associated variants
  effect_true <- continuous_effects_range_true
  
  ## The range of effect sizes for the other variants
  effect_null <- continuous_effects_range_null
  
  ## Generate some correlated errors
  errors <- MASS::mvrnorm(n = nrow(genotype_matrix),
                          mu = rep(0, n_corr_phen),
                          Sigma = corr_matrix)
  
  ## How many variants are truly associated?
  n_variants_true <- round(n_variants*proportion_true)
  
  ## Randomly select the truly associated variants
  variants_true <- apply(X = matrix(n_variants_true),
                         MARGIN = 1,
                         FUN = function(x){
                           tmp = sample(x = 1:n_variants, size = x)
                           tmp = tmp[ order(tmp) ]
                           return(tmp)
                         })
  
  ## Identify the other unrelated variants
  variants_null <- matrix( rep(1:n_variants, times = n_corr_phen),
                           ncol = n_corr_phen)
  variants_null <- apply(X = matrix(1:n_corr_phen),
                         MARGIN = 1,
                         FUN = function(x){
                           tmp = variants_null[, x]
                           tmp[ tmp %in% variants_true[, x] ] = NA
                           return(tmp)
                         })
  
  # ## Save the variant information
  # vars <- paste0("Variant_True_X", 1:n_corr_phen)
  # for(i in 1:length(vars)){
  #   
  #   ## Bring over the null IDs
  #   tmp <- variants_null[, i]
  #   
  #   ## Change the lables
  #   tmp[ !is.na(tmp) ] <- "Null Variant"
  #   tmp[ is.na(tmp) ] <- "True Variant"
  #   
  #   ## Save the results
  #   dat_variants[, vars[i] ] <- tmp
  # }
  # rm(i, vars, tmp)
  
  ## Save the results
  tmp_vars <- paste0("Variants_X", 1:n_corr_phen)
  for(i in 1:length(tmp_vars)){
    tmp <- variants_null[, i]
    tmp[ !is.na(tmp) ] <- "Null"
    tmp[ is.na(tmp) ] <- "True"
    simulation_details[, tmp_vars[i] ] <- tmp
  }
  rm(i, tmp)
  
  ## Randomly generate the effects
  tmp_true <- apply(X = matrix(n_variants_true),
                    MARGIN = 1,
                    FUN = function(x){ return( runif(n = x, min = effect_true[1], max = effect_true[2]) ) } )
  tmp_null <- apply(X = matrix(n_variants_true),
                    MARGIN = 1,
                    FUN = function(x){ runif(n = (n_variants - x), min = effect_null[1], max = effect_null[2]) } )
  
  ## Combine the effects into a matrix
  effects <- variants_null
  for(i in 1:n_corr_phen){
    effects[ !is.na(effects[, i]), i] <- tmp_null[, i]
    effects[ is.na(effects[, i]), i] <- tmp_true[, i]
  }
  rm(i)
  
  ## Remove R objects
  rm(tmp_true, tmp_null)
  
  ## Save the results
  tmp_vars <- paste0("Effects_X", 1:n_corr_phen)
  for(i in 1:length(tmp_vars)){
    simulation_details[, tmp_vars[i] ] <- effects[, i]
  }
  rm(i, tmp_vars)
  
  ## Transform the proportions of variance explained
  effects_transformed <- matrix(NA, nrow = nrow(effects), ncol = ncol(effects))
  tmp_vars <- paste0("Effects_Transformed_X", 1:n_corr_phen)
  tmp_vars_orig <- paste0("Effects_X", 1:n_corr_phen)
  for(i in 1:nrow(simulation_details)){
    
    ## Identify the major and minor allele frequencies
    p <- simulation_details$Variant_MAF[i]
    q <- (1 - p)
    
    ## Transform the effects
    for(k in 1:length(tmp_vars)){
      effects_transformed[i, k] <- sqrt( ( simulation_details[i, tmp_vars_orig[k] ] / (1 - simulation_details[i, tmp_vars_orig[k] ] ) ) / (2*p*q) )
    }
  }
  rm(i, p, q, tmp_vars_orig)
  
  ## Randomly make some effects negative
  effects_transformed <- apply(X = effects_transformed,
                               MARGIN = 2,
                               FUN = function(x){ return( x * sample(x = c(-1, 1), size = length(x), replace = TRUE) ) } )
  
  ## Save the results
  for(i in 1:length(tmp_vars)){
    simulation_details[, tmp_vars[i] ] <- effects_transformed[, i]
  }
  rm(i, tmp_vars)
  
  ## Randomly generate proportions of variance explained
  effects_binary <- runif(n = n_corr_phen,
                          min = binary_effect_on_continuous_range[1],
                          max = binary_effect_on_continuous_range[2])
  
  ## Randomly make some effects negative
  effects_binary <- ( effects_binary * sample(x = c(-1, 1), size = length(effects_binary), replace = TRUE) )
  
  ## Make into a matrix
  effects_binary <- matrix(c(effects_binary[1], 0, 0, effects_binary[2]),
                           ncol = n_corr_phen,
                           byrow = TRUE)
  
  ## Make the binary phenotype a matrix
  z_binary_matrix <- matrix(rep(z_binary, times = n_corr_phen),
                            ncol = n_corr_phen)
  
  ## Make the geno object a matrix
  geno_matrix <- as.matrix(genotype_matrix)
  
  ## Genreate the simulated data
  X <- (geno_matrix %*% effects_transformed) + (z_binary_matrix %*% effects_binary) + errors

  ## Remove R objects
  rm(errors, z_binary_matrix, geno_matrix)
  
  ## Combine the simualte data into a data.frame
  dat <- data.frame(Z_Cont = z,
                    Z_Binary = z_binary,
                    X1 = X[, 1],
                    X2 = X[, 2],
                    stringsAsFactors = FALSE)
  
  ## Remove row.names
  row.names(dat) = NULL
  
  ## Return the simulated data
  return(dat)
  
}