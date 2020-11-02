##########################################################################
############### Specific gene plotting function ##########################
##########################################################################
Gene.trend <- function(Dataset,gene, color.by = "", Use.scale.data = F){
  data <- as.data.frame(Dataset@meta.data$Phase)
  data$PseudotimeScore <- Dataset@meta.data$Pseudotime
  
  if(Use.scale.data == T){
    data$Gene <- Dataset@scale.data[gene,]
  } else {data$Gene <- Dataset@data[gene,]
  }
  
  if(color.by == "Cluster"){
    data$Cluster <- Dataset@meta.data$Transfered.Ident
    p <- ggplot(data=data, aes(x=PseudotimeScore, y=Gene, color=Cluster)) +
         geom_point(size = .5) +
         ylim(0,max(data$Gene)) +
         geom_smooth(method="loess", n= 30, color="red", fill="grey") +
         ggtitle(gene)
  } else {
    data$lineages <- Dataset@meta.data$Lineage
    p <- ggplot(data=data, aes(x=PseudotimeScore, y=Gene, color=lineages)) +
         geom_point(size = .5) +
         scale_color_manual(values= c("#cc391b","#026c9a")) +
         ylim(0,max(data$Gene)) + 
         geom_smooth(method="loess", n= 30,  aes(colour=lineages), fill="grey") +
         ggtitle(gene)
  }
  return(p)
}
##########################################################################
############### Define function to plot several genes ####################
##########################################################################

Plot.Genes.trend <- function(Dataset, genes, color.by, Use.scale.data){
  pList <- mapply(FUN = Gene.trend, gene = genes,
                  MoreArgs = list(Dataset = Dataset, color.by = color.by, Use.scale.data),
                  SIMPLIFY = FALSE)
  print(x = cowplot::plot_grid(plotlist = pList, ncol = 2))
} 


