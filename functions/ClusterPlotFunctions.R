
########### Cluster.dotplot ############################
Cluster.dotplot <- function(Dataset,
                            Marker.genes,
                            min.expression,
                            percent.min,
                            maxdot.size){
  #Extract and transpose the matrix of marker genes expression
  data.to.plot <- data.frame(t(as.matrix(Dataset@data[Marker.genes,])))
  
  data.to.plot$Cell <- rownames(data.to.plot)
  data.to.plot$id <- Dataset@ident
  
  #Reshape the dataframe
  data.to.plot <- data.to.plot %>% gather(key = Marker.genes, value = expression, -c(Cell, id)) 
  
  #For each genes in each cluster calculate mean expression and percent cell with norm expression > 0
  data.to.plot <- data.to.plot %>%
		  group_by(id, Marker.genes) %>% 
    		  summarize(avg.exp = mean(expm1(x = expression)), pct.exp = length(x = expression[expression > min.expression]) / length(x = expression) ) 
  
  data.to.plot <- data.to.plot %>% ungroup() %>%
    group_by(Marker.genes) %>% 
    mutate(avg.exp.scale = scale(x = avg.exp)) %>%
    mutate(avg.exp.scale = MinMax(data = avg.exp.scale, max = 2, min = -2)) # add column with scaled expression values from -2 to 2
  
  data.to.plot$genes.plot <- factor(x = data.to.plot$Marker.genes, levels = rev(x = Marker.genes)) #Put gene names as factor 
  
  data.to.plot$pct.exp[data.to.plot$pct.exp < percent.min] <- NA #Set to Na if less than percent.min of cells expresse the gene
  data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
  
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
	    geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) + # modify the colors by if want to color by domains or cluster ident
	    scale_size_area(max_size= maxdot.size) + # Scale the radius of the dot from 0 to 6 
	    scale_x_discrete(position = "top") +
	    theme(axis.text.x = element_text(angle = 90, vjust = 1), axis.title.y = element_blank()) +
	    xlab("") + ylab("") +
	    scale_colour_gradientn(colours = brewer.pal(11,"RdPu"))
} 
