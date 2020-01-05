


########################################################
## Fit an lme() Model and Calculate the Fitted Values ##
########################################################

compute_lme_fitted_values <- function(Y_Vector = Y_vec,
                                      X_Matrix = X,
                                      Random_Basis = random_basis,
                                      Variant_Positions_01 = variant_positions_01){
  
  ## Packages
  library(nlme)
  require(nlme)
  
  ## Create the fixed part of the model
  X_full <- kronecker(X_Matrix, cbind(1, Variant_Positions_01, Variant_Positions_01^2))
  
  ## Calculate the random part of the model
  Z_full <- kronecker(X_Matrix, Random_Basis)
  
  ## Fit the model
  group <- rep(1, length(Y_Vector))
  mod <- try( lme(Y_Vector ~ X_full - 1,
                  random = list(group = pdIdent(~ Z_full - 1) ) ) )
  if( class(mod) == "try-error" ){
    mod <- lme(Y_Vector ~ X_full - 1,
               random = list(group = pdIdent(~ Z_full - 1) ),
               control = lmeControl(opt = "optim") )
    print("Used the 'optim' method")
  }
  
  ## Calculate the fitted values
  return( as.vector( cbind(X_full, Z_full) %*% t( coef(mod) ) ) )
}


