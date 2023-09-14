# The purpose of this script is to determine whether 
# the masses among avian subclades predict species richness (a proxy for origination rate).
# In short; do small-bodied avian subclades speciate more quickly?
#
# This code was written by Dr. Andrew Orkney, parsed by Prof. Brandon Hedrick, 
# and draws upon existing functions from dependent R packages

library(caper) # 1.0.2 
# Load required analytical package. 

# Set the work directory
setwd('C:/Users/Lab/Documents/AOrkney/bird_scripts_and_data_15_02_2023')
# The user will need to customise this

# Load centroid sizes for study birds. 
load('Csize.22.10.2022.RData')
# In our study, we defined this object as an array called 'GPA.Csize'
# Be aware the names do not yet match the names in the phylogeny we will use.

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

bird.data <- cbind(names(masses),log10(masses), rep(1, length(masses)))
bird.data <- data.frame(bird.data)
bird.data[,2] <- as.numeric(bird.data[,2])
bird.data[,3] <- as.numeric(bird.data[,3])
colnames(bird.data)[1]<- 'binomial' 
colnames(bird.data)[2]<- 'masses'
colnames(bird.data)[3] <- 'species.rich'
# Compile a dataframe. 
# All tips in the phylogeny represent lineages with a single species, 
# so species richness at all tips is 1. 

bird.comparative <- comparative.data(pruned.tree, bird.data, binomial, na.omit=T)
# Organise the dataframe ahead of analysis

birdBodySize <- macrocaic(species.rich ~ masses, data=bird.comparative)
summary(birdBodySize)
# There is strong evidence that bird lineage origination rate depends on body mass 
# The coefficient is negative, therefore subclades dominated by small birds are more likely to be species rich 


