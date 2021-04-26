######################### Smoothing function using gam #########

Landscape.gam.Smoothing <- function(gene,
                                    data,
                                    newgrid,
                                    scale.exp){
  
  df <- data[,c("RC.DV.axis1","RC.DV.axis2", gene)]
  colnames(df) <- c("RC.DV.axis1", "RC.DV.axis2", "z")
  
  # Fit a gam regression 
  mod <- gam(z ~ s(RC.DV.axis1, RC.DV.axis2),
             family= gaussian(),
             data = df,
             method = "REML") #GCV.CP
  
  
  # Predict the new landscape
  smoothed.landscape <- predict(mod,
                                newdata= newgrid,
                                se=TRUE)
  if (scale.exp == T) {
    # Smooth the fitted values
    fit <- scale(smoothed.landscape$fit)
  } else {
    fit <- smoothed.landscape$fit
  }
  
  
  
  return(fit)
}


Gene.smooth.landscape.gam <- function(data, scale.exp = T){
  
  # New coordinates grid
  x = seq(min(data$RC.DV.axis1), max(data$RC.DV.axis1), length.out = 200)
  y = seq(min(data$RC.DV.axis2), max(data$RC.DV.axis2), length.out = 200)
  newgrid <- expand.grid(RC.DV.axis1 = x, RC.DV.axis2 = y, KEEP.OUT.ATTRS = FALSE)
  
  fitlist <- mapply(FUN = Landscape.gam.Smoothing, gene = colnames(data)[3:ncol(data)],
                    MoreArgs = list(data = data, scale.exp =scale.exp, newgrid = newgrid),
                    SIMPLIFY = FALSE)
  
  newgrid <- cbind(newgrid, as.data.frame(fitlist))
  
  return(newgrid)
  
  
}


######## Function to generate pheatmap of the gene landscape #####
Landscape <- function(gene,
                      data,
                      col.pal){
  df <- data[,c("RC.DV.axis1","RC.DV.axis2", gene)]
  colnames(df) <- c("RC.DV.axis1", "RC.DV.axis2", "z")
  data_wide <- tidyr::spread(df, RC.DV.axis1, z)
  rownames(data_wide) <- round(data_wide$RC.DV.axis2,3)
  data_wide <- data_wide[,-1]
  colnames(data_wide) <- round(as.numeric(colnames(data_wide)),3)
  
  if (!col.pal == "viridis") {
    p <- ggplotify::as.ggplot(pheatmap::pheatmap(data_wide[200:1,],
                                                 cluster_rows= F,
                                                 cluster_cols = F,
                                                 scale = "none",
                                                 show_colnames = F,
                                                 show_rownames = F,
                                                 border_color = NA,
                                                 legend = F,
                                                 color = colorRampPalette(rev(brewer.pal(n = 10, name = col.pal)))(100),
                                                 breaks = seq(-2.5,2.5, length.out = 100),
                                                 main = gene,
                                                 silent = T)) 
    
    
  } else{
    p <- ggplotify::as.ggplot(pheatmap::pheatmap(data_wide[200:1,],
                                                 cluster_rows= F,
                                                 cluster_cols = F,
                                                 scale = "none",
                                                 show_colnames = F,
                                                 show_rownames = F,
                                                 border_color = NA,
                                                 legend = F,
                                                 color = viridis::viridis(100),
                                                 breaks = seq(-2.5,2.5, length.out = 100),
                                                 main = gene,
                                                 silent = T)) 
    
  }
  
  
  return(p)
  
}


######## Function to plot several gene landscapes ########
Plot.landscape <- function(data,
                           genes,
                           ncols,
                           col.pal){
  
  pList <- mapply(FUN = Landscape, gene = genes,
                  MoreArgs = list(data= data, col.pal = col.pal),
                  SIMPLIFY = FALSE)
  
  print(x = cowplot::plot_grid(plotlist = pList, ncol = ncols))
} 
