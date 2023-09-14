# This script is part of a project that seeks to determine how variation in body mass structures patterns of 
# interspecific evolutionary covariance in skeletal proportions over the avian skeleton. 
# We should expect locomotory structures, such as wings, to become more integrated as mechanical stress increases, 
# which might have implications for the pathways of evolutionary change that are available to avian subclades exploring
# different ranges of body masses.
#
# The purpose of this script is to produce a plot
# That illustrates how avian body mass is distributed across phylogeny. 
# (subplot a) and perform tests for phylogenetic signal (Blomberg's K and Pagel's Lambda)
#
# In addition to this, raw variance structure for the different
# bone centroid sizes of the avian skeleton need to be plotted 
# as a function of mass.
# (subplot b)
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


dendr <- dendro_data(as.dendrogram(force.ultrametric(pruned.tree))) 
basic <- as.dendrogram(force.ultrametric(pruned.tree))
# Prepare the phylogeny for phylogram plot production
# The phylogeny differed slightly from an ultrametric topology.

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

den.tree <- as.dendrogram(force.ultrametric(pruned.tree))
lab.dat<-dendro_data(den.tree)$labels
lab.dat$label<- gsub('_.*','',lab.dat$label)
# Extract text to label genera 

df.plot <- as.data.frame( log10(masses[match(pruned.tree$tip,names(masses))]) )
colnames(df.plot)<-'mass'
# Prepare a dataframe of log-transformed bird masses for plotting


fit<-fastAnc(pruned.tree,log10(masses)[pruned.tree$tip],vars=TRUE,CI=TRUE)
# Estimate the ancestral character states for each node
# We will colour each line according to the node colour 


n.col <- rep(NA, length(dendr$segments[,1]))
n.col[which(dendr$segments$x==dendr$segments$xend)] <- fit$ace[match(pruned.tree$edge[,1],names(fit$ace))]
n.col[which(dendr$segments$x==dendr$segments$xend)-1] <- fit$ace[match(pruned.tree$edge[,1],names(fit$ace))]

dendr.col<- cbind(dendr$segments, n.col  )
colnames(dendr.col)[5]<-'col'
dendr.col <- data.frame(dendr.col)
# Set information to colour tip and branch colours (these will indicate body mass variety over phylogeny)

BlomK <- phylosig(tree=pruned.tree, x= log10(masses)[pruned.tree$tip], method='K',test=T )[[1]]
# Blomberg's K statistic is 1.07
lambda <- phylosig(tree=pruned.tree, x= log10(masses)[pruned.tree$tip], method='lambda',test=T )[[1]]
# Pagel's lambda is 0.87
# These results are compatible with a Brownian walk model of body mass evolution across the studied birds. 


scale<-.1
# This defines the scale of the root to tip length of the phylogeny.
plotA <- 
ggplot()+
geom_col(data = df.plot,aes(x=1:dim(df.plot)[1], y=rep(-3.6,dim(df.plot)[1])), fill='black', linewidth=1, colour='black' )+
geom_col(data = df.plot,aes(x=1:dim(df.plot)[1], y=rep(-3.3,dim(df.plot)[1]), fill=mass, col=mass), linewidth=1 )+
scale_fill_gradient2(low='blue',high='red',midpoint=mean(df.plot$mass))+
scale_color_gradient2(low='blue',high='red',midpoint=mean(df.plot$mass))+
geom_segment(data=dendr.col, aes(x = x, y = y *scale, xend = xend, yend = yend *scale ),linewidth=2, col='black')+
geom_segment(data=dendr.col, aes(x = x, y = y *scale, xend = xend, yend = yend *scale, col=col ),linewidth=1)+
scale_color_gradient2(low='blue',high='red',midpoint=mean(df.plot$mass))+
	geom_text(data= lab.dat, aes(x=x, y=y, label=label), angle=createAngleHJustCols(lab.dat)[['angle']], hjust=createAngleHJustCols(lab.dat)[['hjust']],size = 5 / .pt )+
labs(fill=expression(paste(log[10],'(',italic(mass),')',' [g]')), col= expression(paste(log[10],'(',italic(mass),')',' [g]')))+
	scale_y_reverse(expand = c(0.2, 0)) + 
 coord_polar()+
	geom_text( aes(x= c(75), y= c(-5.5)  ) , label= as.expression(bquote(K%~~%.(round(BlomK,digits=2))~','~lambda%~~%.(round(lambda,digits=2)))), size=20 /.pt )+
theme(legend.position=c(0.5,0.88),
	legend.direction='horizontal',
	legend.title=element_text(size= 20 ),
	legend.text=element_text(size=15),
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
		#plot.margin=margin(rep(-10,4)),
		legend.box.spacing = unit(0, "pt"),
		plot.margin=unit(c(-0.05,-0.10,-0.05,-0.10), "null"),
       )
# This is a coloured phylogeny plot of body mass over the study bird taxa.


# A K > 1 indicates that variance tends to be partitioned among clades more strongly than expected 
# under Brownian motion; clades have therefore diverged more than expected 
# By contrast, the pagel's lambda is less than expected under Brownian motion, though. 
#
# In sum, the combined interpretation of these two indices probably means phylogenetic comparative methods,
# that use a variance-covariance matrix assuming Brownian motion to account for phylogeny, are probably
# representing this process fairly in birds. 

# The next task is to assess whether there is 
# heteroskedasticity in the variance structure of the 
# residual centroid sizes which will be the subject of integration tests.

# Define a sundry function to interract with data objects such as GPA.coords
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

# Define a function to compute the residual skeletal element centroid cizes, after the affects of allometry are removed.
# We want to account for 'allometry' (size-dependent scaling between species), because it could itself be a major driver of trait integration.
get.residual.Csize <- function( array, masses, phylogeny, taxa ){ # The function takes a shape array, mass vector, phylogeny and list of taxa
	species<-taxa
	allometry.Csize <- list() # Dumby variable to receive allometric models
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Phylogeny pruned to desired taxa
	for(i in 1:length(array) ){ # For each skeletal element 
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ i ]][taxa] ), species = taxa  ) # Define a data frame
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( 1, phy=newphy, form= ~ species ), data=df ) # Compute an allometric model
	}
	names(allometry.Csize) <- names(array) # Ensure names of allometric models match the skeletal elements they describe 
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c ) # Extract model residuals
	return(residual_Csize) # Return residuals
} # Conclude function
# This function uses phylogenetic Generalised Least Squares to compute the allometrically adjusted centroid sizes of skeletal elements.
# The user must provide the array of original bird skeletal element centroid sizes, bodymasses, phylogeny and taxa of interest. 

# Compute residual centroid sizes by applying function:
residual.Csize <- get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = pruned.tree, taxa=names(masses) )
# Which taxa do we wish to investigate? 
# In this example, we will investigate all birds. 

# An array of residual centroid sizes is now available for investigation. 
# 
# We will compute a proxy of variance over a running window of 30 taxa organised by mass, 
# and then plot this against mass.
# The caTools package has a function 'runquantile' that can be used to produce the 
# the interquartile range of a running bin of 30 birds 
# This will need to be plotted as four subplots- one for bones within each module, 
# and each will require their own legend to distinguish line types. 
# They can share a common axis for mass, but they will have different axes
# for the range of residuals 

library(caTools) # 1.18.2

inter.quart <- function(input,k){
	return(runquantile(x=input[order(masses)],k=k,probs=.75,endrule='quantile')-
	runquantile(x=input[order(masses)],k=k,probs=.25,endrule='quantile'))
}
# Compute the interquartile range (range between first and third quartile of variance in skeletal proportions)
# over body mass.
# k = number of values included in each interquartile range
# input = ordered masses from lowest to highest

df<- cbind(rep(masses[order(masses)],length(residual.Csize)), 
rep(names(residual.Csize), each=length(masses)),
unlist(lapply(residual.Csize, inter.quart,k=30))
)

df<-as.data.frame(df)
colnames(df)<-c('mass','bone','range')
df$mass<-as.numeric(df$mass)
df$range<-as.numeric(df$range)
# Compiled results into a dataframe.

# Produce plots of variance structure with body mass for individual modules in the avian body plan.
bi<-
ggplot( data= df[which(df$bone=='skull' | df$bone=='mandible'),] )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=0.5 )+
geom_line( aes(x=log10(mass),y=range, col=bone, linetype=bone), linewidth=1 )+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
bii<-
ggplot( data= df[which(df$bone=='humerus' | df$bone=='ulna' | df$bone=='radius' | df$bone=='carpometacarpus'),] )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=0.5 )+
geom_line( aes(x=log10(mass),y=range, col=bone, linetype=bone), linewidth=1 )+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
biii<-
ggplot( data= df[which(df$bone=='scapula' | df$bone=='coracoid' | df$bone=='sternum' | df$bone=='synsacrum'),] )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=0.5 )+
geom_line( aes(x=log10(mass),y=range, col=bone, linetype=bone), linewidth=1 )+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
biv<-
ggplot( data= df[which(df$bone=='femur' | df$bone=='tibiotarsus' | df$bone=='tarsometatarsus'), ] )+
geom_line( aes(x=log10(mass),y=range, col=bone), linewidth=0.5 )+
geom_line( aes(x=log10(mass),y=range, col=bone, linetype=bone), linewidth=1 )+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())


library(ggpubr) # 0.6.0
library(grid) # base
# These packages contain functions that we can use to combine subplots:

quadrate <- ggarrange(bii,bi,biii,biv,ncol=2, nrow=2, labels=c('Wing','Head','Trunk','Leg'), hjust=-1.5)
annotated<- annotate_figure(quadrate, left = textGrob(expression(paste('Interquartile range ',log[10],'(',italic(rCsize),')',' ')), rot = 90, vjust = 0.5,  gp = gpar(cex = 1.8)),
bottom= textGrob(expression(paste(log[10],'(',italic(mass),')',' [g]')), gp = gpar(cex = 1.8)))

dev.new(width=18,height=9,unit=cm)
ggarrange(plotA, annotated, labels=c('a','b'),font.label=list(size=25))

setwd('C:/Users/Lab/Documents/Aves_2023')
# Set the work directory to the location where you wish to save the combined plots.
ggsave(filename='diagnostic_plot_08_09_2023.pdf')
# Save the plots 

