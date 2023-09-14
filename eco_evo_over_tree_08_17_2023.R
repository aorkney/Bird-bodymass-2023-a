# The purpose of this script is to compute evolutionary models of 
# Brownian variance of the exploration of ecological trait space.
#
# We will explicitly test whether birds with lower body masses 
# evolve at a different rate. 
#
# This code was written by Dr. Andrew Orkney, parsed by Prof. Brandon Hedrick, 
# and draws upon existing functions from dependent R packages

library(geomorph) # 4.0.5 
library(cluster) # 2.1.4
library(ape) # 5.7-1
# Load requisite packages 

# Set the work directory
setwd('C:/Users/Lab/Documents/AOrkney/bird_scripts_and_data_15_02_2023')
# The user will need to customise this

# Load some requisite datasets
load('tree.22.10.2022.RData')
load('tree.names.22.10.2022.RData')
load('masses.22.10.2022.RData')

# Ensure names match the phylogeny. 
names(masses) <- tree_names

# Define upper and lower quartiles of avian body mass 
lower.quartile <- names(which(masses <= quantile(masses)[2])) 
# The 25% smallest birds
upper.quartile <- names(which(masses >= quantile(masses)[4]))
# The 25% largest birds

# The user could also set a bespoke threshold as follows
#lower.quartile <- names(which(masses <= 300))
#upper.quartile <- names(which(masses > 300))

short.tree <- keep.tip(pruned.tree,c(lower.quartile,upper.quartile))
# The phylogeny has been pruned to only the specific taxa we need. 

# The analyses in this script will require categorial flight style information that describes the bird species' ecologies;
# These are represented as a matrix of binary scores (they are not mutually exclusive)
# that will need to be transformed into a continuous representation of ecological space. 

# Load the flight style data
flight<- read.csv('flight_masses_22_10_2022_plus_A.csv') 
flight.rownames<-flight$X
flight<-flight[,-c(1,2,14)]
rownames(flight)<-flight.rownames
flight<-flight[short.tree$tip,]
flightdist <- daisy(flight, metric='gower', stand=F, type=list(asymm=c(1:11)))
 # Compute distance matrix
flightdist[which(is.na(flightdist)==T)]<-1
 # Replace Na values with 1 
flightpcoa <- cmdscale(flightdist,eig=T)
 # Compute a Principal Coordinate Analysis to make the data continuously distributed 
k <- max(which(flightpcoa$eig/sum(flightpcoa$eig)>=0.05))
 # Only the lowest axes of variation each explaining >5% of variance were retained.
full.flightpcoa <- cmdscale(flightdist, k=k, eig=T) # Perform a truncation 
# Done 

gp <- rep('0',length(short.tree$tip))
names(gp) <- short.tree$tip
gp[lower.quartile]<-'lower'
gp[upper.quartile]<-'upper'
# Define a group vector for birds belonging to the upper and lower quartiles of body mass distribution

compare.evol.rates(A=full.flightpcoa$points[short.tree$tip,], phy=short.tree,
  method="simulation",gp=gp,iter=999)
# The compare.evol.rates function will test whether Brownian covariances differ across the partition. 

# The interested user may wish to
# repeat these tests using residual carpometacarpus length, rather than flight style variety, 
# as a representation of ecological space (because residual carpometacarpus length is a continuous variable
# with a known strong relationship to ecological variety across birds. 


