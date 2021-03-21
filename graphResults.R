
melted_df_Plot <- function(CDE,
                           z1_grid,
                           z2_grid){
  
  df_result = data.frame('CDE' = as.numeric(CDE),
                         'z1_grid' = rep(z1_grid, nrow(CDE)),
                         'z2_grid' = rep(z2_grid, each = ncol(CDE)))
  
  return(df_result)
}

contourGraph <- function(melted_cormat,
                        melted_cormat_real = NULL,
                        x_lim,
                        y_lim,
                        tol = 1E-6){
  
  melted_cormat$CDE[melted_cormat$CDE < tol] = 0
  
  
  plot <- ggplot() +
    geom_contour(data = melted_cormat,
                 aes(x = z1_grid, y = z2_grid, z = CDE, color = 'Estimated')) 
  
  if(!is.null(melted_cormat_real)){
  plot <- plot + geom_contour(data = melted_cormat_real,
                 aes(x = z1_grid, y = z2_grid, z = CDE, color = 'Real'))
  }
  
  plot <- plot + xlim(x_lim[1], x_lim[2]) + ylim(y_lim[1], y_lim[2])
    
  return(plot)
}
