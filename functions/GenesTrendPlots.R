##################################################################################
##################################################################################
##################################################################################

Plot.gene.trend <- function(Dataset, 
                            gene,
                            x.intercept = NULL,
                            Use.scale.data = F,
                            Axis= "DorsoVentral.Score"){
  
  if(Axis != "DorsoVentral.Score"){
    data <- data.frame(Cluster = paste0("Cluster",as.character(Dataset@ident)),
                       RostroCaudal.score =  Dataset@meta.data$RostroCaudal.score)
    
    if(!Use.scale.data){
      data$Gene <- Dataset@data[gene,]
    } else {data$Gene <- Dataset@scale.data[gene,]
    }
    
    p <- ggplot(data=data, aes(x=-RostroCaudal.score, y=Gene, color=Cluster)) +
      geom_point() +
      scale_color_manual(values= c("#68b041", "#e3c148", "#b7d174", "#e46b6b")) + 
      geom_smooth(method="loess", n= 30, color="red", fill="grey") + 
      ggtitle(gene) +
      ylim(0,NA) +
      theme(legend.position="none")
    
    if(!is.null(x.intercept)){
      p <- p + geom_vline(xintercept = x.intercept, colour = "red", linetype = 2)
    }
    
    return(p)
    
  } else {
    data <- data.frame(Cluster = paste0("Cluster",as.character(Dataset@ident)),
                       DorsoVentral.Score =  Dataset@meta.data$DorsoVentral.Score)
    
    if(!Use.scale.data){
      data$Gene <- Dataset@data[gene,]
    } else {data$Gene <- Dataset@scale.data[gene,]
    }
    
    p <- ggplot(data=data, aes(x=DorsoVentral.Score, y=Gene, color=Cluster)) +
      geom_point() +
      scale_color_manual(values= c("#68b041", "#e3c148", "#b7d174", "#83c3b8", "#009fda", "#3e69ac", "#e46b6b")) + 
      geom_smooth(method="loess", n= 30, color="red", fill="grey") + 
      ggtitle(gene) +
      ylim(0,NA) +
      theme(legend.position="none")
    
    if(!is.null(x.intercept)){
      p <- p + geom_vline(xintercept = x.intercept, colour = "red", linetype = 2)
    }
    
    return(p)
  }
  
}

##################################################################################
##################################################################################
##################################################################################

Plot.Genes.trend <- function(Dataset,
                             genes,
                             Axis= "DorsoVentral.Score",
                             x.intercept= NULL,
                             Use.scale.data = FALSE){
  pList <- mapply(FUN = Plot.gene.trend, gene = genes,
                  MoreArgs = list(Dataset = Dataset,
                                  Axis= Axis,
                                  x.intercept= x.intercept,
                                  Use.scale.data = Use.scale.data),
                  SIMPLIFY = FALSE)
  print(x = cowplot::plot_grid(plotlist = pList, ncol = 2))
} 


Plot.Cluster.trend <- function(Datasrt, #The Seurat object
                               Which.cluster, #number of the cluster to plot
                               clust.list, #list from the clustering results
                               group.by, # smooth by all genes by lineages or gene by gene
                               span, #Interger value parameter of the smoother argument of geom_smooth
                               Smooth.method, #Smoothing method (function) use by geom_smooth e.g. "auto", "lm", "glm", "gam", "loess"
                               Use.scale.data #Use scaled expression matrix
){
  
  clust.vec <- clust.list
  Cluster.Genes <- names(clust.vec[clust.vec == Which.cluster])
  
  if(!Use.scale.data){
    Expression <- as.data.frame(t(as.matrix(Datasrt@data[Cluster.Genes,])))
  } else {Expression <- as.data.frame(t(as.matrix(Datasrt@scale.data[Cluster.Genes,])))
  }
  
  Expression$Cluster = paste0("Cluster",as.character(Datasrt@ident))
  Expression$DorsoVentral.Score =  Datasrt@meta.data$DorsoVentral.Score
  Melt.Data <- melt(Expression, id=c("DorsoVentral.Score", "Cluster")) ; head(Melt.Data)
  
  if (group.by == "Cluster"){
    p <- ggplot(data=Melt.Data,
                aes(x=DorsoVentral.Score, y=value, colour=Cluster)) +
      geom_smooth(method= Smooth.method,
                  span = span, n= 30, fill="grey", aes(colour=lineages))+
      ggtitle(paste0("Cluster ", Which.cluster, " (", length(Cluster.Genes), " genes)"))
  } else if (group.by == "genes") {
    p <- ggplot(data=Melt.Data,
                aes(x=DorsoVentral.Score, y=value, colour=variable)) +
      geom_smooth(method= Smooth.method,
                  span = span, n= 30, fill="grey", aes(colour=variable)) +
      ggtitle(paste0("Cluster ", Which.cluster, " (", length(Cluster.Genes), " genes)"))
  } else if (group.by == "global"){
    p <- ggplot(data=Melt.Data,
                aes(x=DorsoVentral.Score, y=value)) +
      geom_smooth(method= Smooth.method,
                  span = span, n= 30, fill="grey") +
      ggtitle(paste0("Cluster ", Which.cluster, " (", length(Cluster.Genes), " genes)"))
  } else { print("group by must be 'Cluster' or 'genes' ")
  }
  
  return(p)
}

Clusters.trend <- function(Datasrt,
                           Which.cluster, #number of the clusters to plot
                           clust.list, #list from the clustering results
                           group.by, #smooth all genes by lineages
                           span, #Interger value parameter of the smoother argument of geom_smooth
                           Smooth.method, #Smoothing method (function) use by geom_smooth e.g. "auto", "lm", "glm", "gam", "loess"
                           Use.scale.data #Use scaled expression matrix
){
  pList <- mapply(FUN = Plot.Cluster.trend,
                  Which.cluster = Which.cluster,
                  MoreArgs = list(Datasrt = Datasrt,
                                  clust.list = clust.list,
                                  group.by = group.by,
                                  span = span,
                                  Smooth.method = Smooth.method,
                                  Use.scale.data = Use.scale.data),
                  SIMPLIFY = FALSE)
  print(x = cowplot::plot_grid(plotlist = pList, ncol = 2))
}
