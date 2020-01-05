


###################################
## Calculate the New F-Statistic ##
###################################

## Inputs for troubleshooting
# RSS_Full <- rss_Full
# RSS_Reduced <- rss_ZBinary_X2
# DF1 <- (4 - 3)
# DF2 <- (50 - 4)
# Covariance_Matrix_Eigen_Values <- eigen_values_corr

calculate_new_f_statistic <- function(RSS_Full,
                                      RSS_Reduced,
                                      DF1,
                                      DF2,
                                      Covariance_Matrix_Eigen_Values){
  
  library(CompQuadForm)
  require(CompQuadForm)
  library(testit)
  require(testit)
  
  ## Calculate the numerator and denominator
  X1 = (RSS_Reduced -  RSS_Full)
  X2 = RSS_Full
  
  ## Define the first and second accuracy parameters to try
  acc_parm_1 <- 1e-20
  acc_parm_2 <- 1e-10
  
  ## Calculate the empirical CDF values (equivalent to p-values)
  p1 = has_warning( davies(X1,
                           lambda = Covariance_Matrix_Eigen_Values,
                           h = rep(DF1, length(Covariance_Matrix_Eigen_Values)),
                           acc = acc_parm_1,
                           lim = 2e6)$Qq )
  p2 = has_warning( davies(X2,
                           lambda = Covariance_Matrix_Eigen_Values,
                           h = rep(DF2, length(Covariance_Matrix_Eigen_Values)),
                           acc = acc_parm_1,
                           lim = 2e6)$Qq )
  
  ## If either produces a warning, then increase the accuracy parameter
  if( p1 ){
    p1 = davies(X1,
                lambda = Covariance_Matrix_Eigen_Values,
                h = rep(DF1, length(Covariance_Matrix_Eigen_Values)),
                acc = acc_parm_2,
                lim = 2e6)$Qq
  }
  else{
    p1 = davies(X1,
                lambda = Covariance_Matrix_Eigen_Values,
                h = rep(DF1, length(Covariance_Matrix_Eigen_Values)),
                acc = acc_parm_1,
                lim = 2e6)$Qq
  }
  
  if( p2 ){
    p2 = davies(X2,
                lambda = Covariance_Matrix_Eigen_Values,
                h = rep(DF2, length(Covariance_Matrix_Eigen_Values)),
                acc = acc_parm_2,
                lim = 2e6)$Qq
  }
  else{
    p2 = davies(X2,
                lambda = Covariance_Matrix_Eigen_Values,
                h = rep(DF2, length(Covariance_Matrix_Eigen_Values)),
                acc = acc_parm_1,
                lim = 2e6)$Qq
  }
  
  ## Map these eCDF values to chi-square distributions
  c1 = qchisq(1-p1, df = DF1)
  c2 = qchisq(1-p2, df = DF2)
  
  ## Return the results
  return( list(
    FStatistic = (c1 / DF1) / (c2 / DF2),
    FStatistic_DF = c(DF1, DF2),
    Pvalue = (1 - pf(((c1 / DF1) / (c2 / DF2)),
                     df1 = DF1,
                     df2 = DF2))
  ))
}


