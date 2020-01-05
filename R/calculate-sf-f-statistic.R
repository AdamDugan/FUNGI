


#######################################################
## Calculate the Shen-Faraway Functional F-Statistic ##
#######################################################

calculate_sf_f_statistic <- function(RSS_Full,
                                     RSS_Reduced,
                                     DF1,
                                     DF2,
                                     VarCov_Matrix_Eigen_Values,
                                     Adjust_DF = TRUE,
                                     Round_To_Closest_Integer = FALSE){
  
  ## Adjust the degrees of freedom if needed
  if( Adjust_DF ){
    
    ## Square the eigen values
    eigen_values_squared <- VarCov_Matrix_Eigen_Values^2
    
    ## Calculate the degrees of freedom adjustment factor
    df_af <- ( ((sum(VarCov_Matrix_Eigen_Values))^2) / (sum(eigen_values_squared)))
    
    ## Calcualte the adjusted degrees of freedom for the F distribution
    DF1 <- df_af*DF1
    DF2 <- df_af*DF2
    
  }
  
  ## Round the degrees of freedom if necessary
  if( Round_To_Closest_Integer ){
    DF1 <- round(DF1)
    DF2 <- round(DF2)
  }
  
  ## Calculate the F-statistic
  f_statistic <- ( (RSS_Reduced -  RSS_Full) / DF1) / (RSS_Full / DF2)
  
  ## Return the results
  return(
    list( FStatistic = f_statistic,
          FStatistic_DF = c(DF1, DF2),
          Pvalue = (1 - pf(f_statistic, df1 = DF1, df2 = DF2) ) )
    )
}


