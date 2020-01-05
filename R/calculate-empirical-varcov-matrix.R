


########################################################
## Calculate the Empitical Variance-Covariance Matrix ##
########################################################

## Inputs for troubleshooting
# Genotype_Matrix = geno_flipped_smooth
# Variant_Positions_01 = variant_positions_01

compute_empirical_varcov_matrix <- function(Genotype_Matrix = geno_flipped_smooth,
                                    Subject_Per_Column = TRUE,
                                    Variant_Positions_01 = variant_positions_01){
  
  ## Load packages
  library(fda)
  require(fda)
  
  ## Rename the input
  tmp_pos = Variant_Positions_01
  tmp_Knots = Variant_Positions_01
  
  ## Transpose the genotype matrix if necessary
  if( !Subject_Per_Column ){ Genotype_Matrix = t(Genotype_Matrix) }
  
  ## Define fda values
  norder = 4
  nbasis = (length(tmp_Knots) + norder - 2)
  Lfdobj = 2
  
  ## Calculate the empirical variance-covariance matrix
  ybasis = create.bspline.basis(rangeval = range(tmp_Knots),
                                nbasis = nbasis,
                                norder = norder,
                                breaks = tmp_Knots)
  yfdPar = fdPar(fdobj = ybasis,
                 Lfdobj = Lfdobj,
                 lambda = 1e-13)
  yfd = smooth.basis(argvals = tmp_pos,
                     y = Genotype_Matrix,
                     fdParobj = yfdPar)$fd
  cvar.fd = var.fd(yfd)
  E.fd = eval.bifd(tmp_pos, tmp_pos, cvar.fd)
  
  ## Extract the eigenvalue decomposition of the correlation matrix
  # S = cov2cor(E.fd)
  # ee = eigen(S)
  # eivec = ee$vectors
  # eigva = ee$values
  
  ## Return the matrix
  if( !Subject_Per_Column ){ return( t(E.fd) ) }
  else{ return(E.fd) }
  
}


