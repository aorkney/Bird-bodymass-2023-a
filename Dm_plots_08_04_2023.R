# The purpose of this script is to produce a combination of plots
# that demonstrate how integration within the head, wing and
# between the wing and trunk, change with increasing body mass.
# This will be achieved by computing major axes of covariance of 
# skeletal proportions within a phylogenetic framework, 
# before determining whether heteroskedasticity (scatter) of the individual
# taxa depends significantly upon body mass. 
#
# This code was written by Dr. Andrew Orkney, parsed by Prof. Brandon Hedrick, 
# and draws upon existing functions from dependent R packages and online contributions from Stackexchange.

library( geomorph ) # 4.0.5
library( ape ) # 5.7-1
library( nlme ) # 3.1-162
library( ggplot2 ) # 3.4.1
library( ggdendro ) # 0.1.23
library( dendextend ) # 1.17.1
library( phytools ) # 1.5-1
# Load required packages

setwd('C:/Users/Lab/Documents/AOrkney/bird_scripts_and_data_15_02_2023')
# Set work directory appropriately (you will need to adjust this if you
# are replicating the work undertaken here)

# Load centroid sizes for study birds. 
load('Csize.22.10.2022.RData')
# In our study, we defined this object as an array called 'GPA.Csize'
# Be aware that the taxon names do not yet match the names in the phylogeny we will use.
# This is because some of the collected taxa were not represented in the phylogeny of Prum et al., 2015, 
# and they were subsequently matched to their nearest relatives. 

# Load phylogenetic tree of birds
load('tree.22.10.2022.RData')
# (Pruned from Prum et al., 2015)
# The object is a phylogeny called 'pruned.tree'.

# Load the tree names required to match birds to the closest genus on the tree.
load('tree.names.22.10.2022.RData')
# The object is a vector called 'tree_names'.

# Load bird masses
load('masses.22.10.2022.Rdata')
# The object is a named numeric vector called 'masses'.
# The names will not yet match those on the phylogeny.

# Let's now ensure names match the phylogeny. 
names(masses) <- tree_names
for( i in 1:length(GPA.Csize) ){
	names(GPA.Csize[[i]]) <- tree_names 
}
# Done


# Data are loaded and names have been set

den.tree <- as.dendrogram(force.ultrametric(pruned.tree))
seg.tree<-dendro_data(den.tree)$segments
# Prepare the phylogeny for phylogram plot production
# The phylogeny differed slightly from an ultrametric topology.

# https://stackoverflow.com/questions/21474388/colorize-clusters-in-dendogram-with-ggplot2
# I have borrowed code that was discussed in this online forum, in order to construct a plot of
# the avian phylogeny.

library(zoo) # 1.8-12
library(dplyr) # 1.1.1
# These packages will be required to plot the phylogeny.

cut <- 20
# We will colour 20 sub-clades of the avian phylogeny, to help orientate the reader.
dendr <- dendro_data(as.dendrogram(force.ultrametric(pruned.tree))) 
clust <- cutree(force.ultrametric(pruned.tree), k= cut)
clust.df <- data.frame(label = names(clust), cluster = clust)
# A dataframe with the constituent clustering information has been produced.

height <- unique(dendr$segments$y)[order(unique(dendr$segments$y), decreasing = TRUE)]
cut.height <- mean(c(height[cut], height[cut-1]))
dendr$segments$line <- ifelse(dendr$segments$y == dendr$segments$yend &
   dendr$segments$y > cut.height, 1, 2)
dendr$segments$line <- ifelse(dendr$segments$yend  > cut.height, 1, dendr$segments$line)
# These lines divide the dendrogram into the roots and the branches that comprise the 20 sub-clades.


dendr$segments$cluster <- c(-1, diff(dendr$segments$line))
change <- which(dendr$segments$cluster == 1)
for (i in 1:cut) dendr$segments$cluster[change[i]] = i + 1
dendr$segments$cluster <-  ifelse(dendr$segments$line == 1, 1, 
             ifelse(dendr$segments$cluster == 0, NA, dendr$segments$cluster))
dendr$segments$cluster <- na.locf(dendr$segments$cluster) 
# These lines ascribe numbers to the sub-clades


clust.df$label <- factor(clust.df$label, levels = levels(dendr$labels$label))
clust.df <- arrange(clust.df, label)
clust.df$cluster <- factor((clust.df$cluster), levels = unique(clust.df$cluster), labels = (1:cut) + 1)
dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")
# These lines assure consistent numbering between segment$cluster and label$cluster

n.rle <- rle(dendr$segments$cluster)
N <- cumsum(n.rle$lengths)
N <- N[seq(1, length(N), 2)] + 1
N.df <- dendr$segments[N, ]
N.df$cluster <- N.df$cluster - 1
# These lines compute positions for sub-clade labels (we will not actually need these)

createAngleHJustCols <- function(labeldf) {        
    nn <- length(labeldf$y)
    halfn <- floor(nn/2)
    firsthalf <- rev(90 + seq(0,360, length.out = nn))
    secondhalf <- rev(-90 + seq(0,360, length.out = nn))
    angle <- numeric(nn)
    angle[1:halfn] <- firsthalf[1:halfn]
    angle[(halfn+1):nn] <- secondhalf[(halfn+1):nn]

    hjust <- numeric(nn)
    hjust[1:halfn] <- 0
    hjust[(halfn+1):nn] <- 1

    return(list(angle = angle, hjust = hjust))
}
# This is a function that I borrowed from Stackexchange: https://stackoverflow.com/questions/38034663/rotate-labels-for-ggplot-dendrogram
# The function rotates labels on fan phylogenies. 

lab.dat<-dendro_data(den.tree)$labels
lab.dat$label<- gsub('_.*','',lab.dat$label)
# Extract text to label genera 

scale<-0.1
# This defines the scale of the root to tip length of the phylogeny.
phylo.plot <- 
ggplot()+
geom_segment(
data=segment(dendr),
  aes(x = x, y = y *scale, xend = xend, yend = yend *scale, size=factor(line), colour=factor(cluster)), 
      lineend = "square", show.legend = FALSE) + 
  	scale_colour_manual(values = c("grey60", rev(hcl.colors(cut,palett='Hawaii')))) +
	scale_size_manual(values = c(.1, 1.1)) +
	geom_text(data= lab.dat, aes(x=x, y=y, label=gsub('_.*','',label)), angle=(createAngleHJustCols(lab.dat)[['angle']]), hjust=(createAngleHJustCols(lab.dat)[['hjust']]),size = 6 / .pt )+
	scale_y_reverse(expand = c(0.2, 0)) + 
	coord_polar()+
	expand_limits(x=0)+
	theme(legend.position='bottom',
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		axis.ticks.x=element_blank(),
       )

# This is a coloured cut of the phylogeny with 20 different subclades coloured
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

# Now that we have the clusters, colours and taxon labels we want to be able to produce subplots in R as described at the start of the code.



# Define a sundry function to interrogate data objects
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

# Define a function to compute landmark constellation centroid sizes, corrected for allometric scaling.
# We want to account for 'allometry' (size-dependent scaling between species), because it could itself be a major driver of trait integration.
get.residual.Csize <- function( array, masses, phylogeny, taxa, bones ){ # Input arguments of shape-data array, mass vector, phylogeny, taxa and bones of interest
	allometry.Csize <- list() # Dumby variable defined to receive allometric models
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Phylogeny pruned to taxa of interest
	for(i in 1:length(array[bones]) ){ # For each bone of interest 
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ bones[i] ]][taxa] ), species = taxa ) # Organise the data into a dataframe
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( 1, phy=newphy, form= ~ species  ), data=df ) # Compute allometric model
	}
	names(allometry.Csize) <- names(array[bones]) # Check correspondance of allometric models with bone names that describe them 
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c ) # Extract allometric model residuals 
	return(residual_Csize) # Return residuals 
}

# Define a function to compute the major axis of trait covariance between two sets of skeletal proportions. 
m.axis.resid <- function( array, masses, bone1, bone2, phylogeny, taxa ){ # This function takes shape data, a mass vector, bone names, phylogeny and taxa of interest inputs
	bones<-c(bone1,bone2) # Define a pairwise combination of bones 
	r.Csize <- get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = phylogeny, taxa= taxa,bones= bones ) # Compute residual centroid size
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Prune phylogeny to taxa of interest 
	two.block <- phylo.integration(r.Csize[[bone1]],r.Csize[[bone2]],phy=newphy) # Compute 2 block partial least squares (integration) 
	residuals<- prcomp(cbind(two.block$XScores[,1],two.block$YScores[,1]))$x[,2] # Find residuals from the major axis of covaraince 
	list<-list()
	list[[1]]<-two.block
	list[[2]]<-residuals
	return(list) # Return residuals 
}


bone1<-'carpometacarpus'
bone2<-'humerus'

output<-m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, names(masses) )
# Find the residuals from the major axis relating bones 1 and 2.
model <- lm(output[[2]] ~log10(masses) )
# Compute a model of the residuals depending on body mass.
car::ncvTest(model, ~log10(masses))$p # 3.1.1
# Apply a Breusch-Pagan test to determine whether heteroskedasticity (scatter) depends on body mass. 
# If it does, then this shows that trait integration varies significantly with body mass.
# We will need to visualise the relationship to determine its sense, however. 

library(caTools) # 1.18.2
# We will now begin to visualise the relationship.
# We need to load this package to use running quantile functions for plotting sigma contour intervals. 

# Compute 1 and 2 sigma confidence intervals for the data distribution
# A running window of k=30 is applied to comptue quantiles 
ribbons <- as.data.frame(cbind( rep(0,length(masses)),
log10(masses)[order(masses)],
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Done 

library(ggrepel) # 0.9.3
# This package can be applied to minimise label overlap.

 axis.title <- 15
 axis.text <- 10
 point.size <- 3
 line.width <- 1/2
 genus.text<-2.5

df <- as.data.frame(cbind(log10(masses),output[[2]], clust[match(names(masses),names(clust))]))

outside <- list()
for(i in 1:dim(ribbons)[1]){
	if(df$V2[order(masses)][i] > ribbons$V3[i] | df$V2[order(masses)][i] < ribbons$V6[i] ){
		outside[[i]] <- match(rownames(ribbons)[i], rownames(df))
	}
}
# Find the bird taxa which fall outside of the 2 sigma confidence intervals

wing.plot <- ggplot()+ # Make plot
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+ # Axial limits
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+ # Confidence envelope 1
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+ # Confidence envelope 2 
geom_point( aes( x=V1, y=V2, fill=factor(V3) ), size=point.size, stroke=1, data=df, col='black', pch=21, show.legend=F)+
  	scale_fill_manual(values = c( rev(hcl.colors(cut,palett='Hawaii')))) +
geom_text_repel( data= df[unlist(outside),], aes( x=V1, y=V2 ), label=gsub('_.*','',rownames(df[unlist(outside),])),size = 6 / .pt )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')) , title='carpometacarpus-humerus *')+ # Axial labels 
# Sundry thematic instructions that regulate plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks.x=element_line(size=2),
      axis.ticks.y=element_blank(),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"),
	title=element_text(size=16))
# This is a plot which illustrates deviation from the major axis of trait covariance between the carpometacarpus and humerus, and how it is 
# structured over body mass. 

# I will desist commenting for the subsequent combinations:

bone1<-'skull'
bone2<-'mandible'

output<-m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, names(masses) )
model <- lm(output[[2]] ~log10(masses) )
car::ncvTest(model, ~log10(masses))$p 3.1.1


# Compute 1 and 2 sigma confidence intervals for the data distribution
ribbons <- as.data.frame(cbind( rep(0,length(masses)),
log10(masses)[order(masses)],
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Done 

df <- as.data.frame(cbind(log10(masses),output[[2]], clust[match(names(masses),names(clust))]))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))

outside <- list()
for(i in 1:dim(ribbons)[1]){
	if(df$V2[order(masses)][i] > ribbons$V3[i] | df$V2[order(masses)][i] < ribbons$V6[i] ){
		outside[[i]] <- match(rownames(ribbons)[i], rownames(df))
	}
}
df[unlist(outside),]
# Find the bird taxa which fall outside of the 2 sigma confidence intervals


head.plot <- ggplot()+ # Make plot
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+ # Axial limits
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+ # Confidence envelope 1
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+ # Confidence envelope 2 
geom_point( aes( x=V1, y=V2, fill=factor(V3) ), size=point.size, stroke=1, data=df, col='black', pch=21, show.legend=F)+
  	scale_fill_manual(values = c( rev(hcl.colors(cut,palett='Hawaii')))) +
geom_text_repel( data= df[unlist(outside),], aes( x=V1, y=V2 ), label=gsub('_.*','',rownames(df[unlist(outside),])),size = 6 / .pt )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')), title='cranium-mandible')+ # Axial labels 
# Sundry thematic instructions that regulate plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks.x=element_line(size=2),
      axis.ticks.y=element_blank(),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"),
	title=element_text(size=16))




bone2<-'scapula'
bone1<-'sternum'

output<-m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, names(masses) )
model <- lm(output[[2]] ~log10(masses) )
car::ncvTest(model, ~log10(masses))$p 3.1.1



# Compute 1 and 2 sigma confidence intervals for the data distribution
ribbons <- as.data.frame(cbind( rep(0,length(masses)),
log10(masses)[order(masses)],
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Done 

df <- as.data.frame(cbind(log10(masses),output[[2]], clust[match(names(masses),names(clust))]))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))

outside <- list()
for(i in 1:dim(ribbons)[1]){
	if(df$V2[order(masses)][i] > ribbons$V3[i] | df$V2[order(masses)][i] < ribbons$V6[i] ){
		outside[[i]] <- match(rownames(ribbons)[i], rownames(df))
	}
}
df[unlist(outside),]
# Find the bird taxa which fall outside of the 2 sigma confidence intervals


trunk.plot <- ggplot()+ # Make plot
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+ # Axial limits
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+ # Confidence envelope 1
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+ # Confidence envelope 2 
geom_point( aes( x=V1, y=V2, fill=factor(V3) ), size=point.size, stroke=1, data=df, col='black', pch=21, show.legend=F)+
  	scale_fill_manual(values = c( rev(hcl.colors(cut,palett='Hawaii')))) +
geom_text_repel( data= df[unlist(outside),], aes( x=V1, y=V2 ), label=gsub('_.*','',rownames(df[unlist(outside),])),size = 6 / .pt )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')), title='scapula-sternum *')+ # Axial labels 
# Sundry thematic instructions that regulate plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks.x=element_line(size=2),
      axis.ticks.y=element_blank(),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"),
	title=element_text(size=16))




bone1<-'carpometacarpus'
bone2<-'sternum'

output<-m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, names(masses) )
model <- lm(output[[2]] ~log10(masses) )
car::ncvTest(model, ~log10(masses))$p 3.1.1


# Compute 1 and 2 sigma confidence intervals for the data distribution
ribbons <- as.data.frame(cbind( rep(0,length(masses)),
log10(masses)[order(masses)],
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Done 

df <- as.data.frame(cbind(log10(masses),output[[2]], clust[match(names(masses),names(clust))]))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))


outside <- list()
for(i in 1:dim(ribbons)[1]){
	if(df$V2[order(masses)][i] > ribbons$V3[i] | df$V2[order(masses)][i] < ribbons$V6[i] ){
		outside[[i]] <- match(rownames(ribbons)[i], rownames(df))
	}
}
df[unlist(outside),]
# Find the bird taxa which fall outside of the 2 sigma confidence intervals


wing.trunk.plot <- ggplot()+ # Make plot
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+ # Axial limits
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+ # Confidence envelope 1
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+ # Confidence envelope 2 
geom_point( aes( x=V1, y=V2, fill=factor(V3) ), size=point.size, stroke=1, data=df, col='black', pch=21, show.legend=F)+
  	scale_fill_manual(values = c( rev(hcl.colors(cut,palett='Hawaii')))) +
geom_text_repel( data= df[unlist(outside),], aes( x=V1, y=V2 ), label=gsub('_.*','',rownames(df[unlist(outside),])),size = 6 / .pt )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')), title='carpometacarpus-sternum *')+ # Axial labels 
# Sundry thematic instructions that regulate plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks.x=element_line(size=2),
      axis.ticks.y=element_blank(),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"),
	title=element_text(size=16))


library(ggpubr) # 0.6.0 
library(grid) # base
# These packages contain functions that we can use to combine subplots:

quadrate <- ggarrange(ncol=2,nrow=2,wing.plot,head.plot,trunk.plot,wing.trunk.plot, labels=c('b','c','d','e'),font.label=list(size=25) , hjust=.7, vjust= 3)

annotated<- annotate_figure(quadrate, left = textGrob(expression(paste( '',italic('D'[m]),'' )), rot = 90, vjust = 0.5,  gp = gpar(fontsize=20)),
bottom= textGrob(expression(paste(log[10],'(',italic(mass),')',' [g]')), gp = gpar(fontsize=20)))

dev.new(width=18,height=9,unit='cm')

ggarrange(phylo.plot,annotated, labels=c('a'),font.label=list(size=25), hjust=-3, vjust= 3 )

setwd('C:/Users/Lab/Documents/Aves_2023')
# Set the work directory to the location where you wish to save the combined plots.
ggsave(filename='Dm_plot_08_04_2023.pdf')
# Save the plots 


