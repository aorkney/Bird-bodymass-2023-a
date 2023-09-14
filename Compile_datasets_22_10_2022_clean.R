# The purpose of this script is to generate
# a suite of data to be used in subsequent analyses for the SICBY
# conference abstract submission, and potentially future research. 

require( geomorph )

load("D:/Documents/Alex_birds/NEW_analysis/2021_Bjarnason_125_SIdata/Supp 3 Datasets for Bjarnason & Benson 2021 revised submission 26 January 2021/Bird_landmarks_clean_26January2021_Bjarnason&Benson taxon set.RData")

landmark.info<-read.csv("D:\\Documents\\Alex_birds\\Bird_landmark_info_23Nov2019.csv" , row.names = 1)

# We are now loading the Avian Phylogeny of Prum et al., (2015)
setwd("D:/Documents/Alex_birds/")
library(ape)
PrumTree <- read.nexus(file = "Avian-TimeTree.tre")

setwd('D:/Documents/Alex_birds/NEW_analysis/')
metadata <- read.csv( "Eco_meta_data.csv" , row.names = 1 )

# Load metadata (this contains information about bird masses and ecological properties)

SL.counts <- "min"

elements <- c( "skull" , "mandible" , "scapula", "coracoid", "sternum" , "humerus", "ulna", "radius", "carpometacarpus", "synsacrum" , "femur", "tibiotarsus", "tarsometatarsus" )

name_list<-list()
for(i in 1:length(elements)){
	name_list[[i]]<-dimnames(compiled.bird.landmarks.list[[ elements[ i ] ]][[ SL.counts ]][[ "output.landmarks.array" ]] )[[ 3 ]]
}
AA<-table(unlist(name_list))
birds <- names( AA )[ AA == length( elements ) ]
birds <- birds[ birds != "Falco_sparverius_UMMZ_154452" ]

GPA.results <- list()
	for( i in 1:length( elements ) )	{
		element <- elements [ i ]
			GPA.results[[ i ]] <- gpagen( compiled.bird.landmarks.list[[ element ]][[ SL.counts ]][[ "output.landmarks.array" ]][ , , birds ] ,  curves = compiled.bird.landmarks.list[[ element ]][[ SL.counts ]][[ "sliders" ]] , ProcD = T )
		}
names( GPA.results ) <- elements

birds[7]<-'Anseranas semipalmata NHMUK 1852.7.22.1'
birds[10]<-'Archilochus_colubris_FMNH_484738'
birds[12]<-'Choriotis_kori_UMMZ_210444'
birds[13]<-'Arenaria_interpres_FMNH_313992'
birds[30]<-'Pipra_erythrocephala_UMMZ_157210'
birds[31]<-'Chaetura_brachyrua_UMMZ_157689'
birds[32]<-'Charadris_vociferus_FMNH_470173'
birds[37]<-'Larus_novae-hollandiae_UMZC_274.C_'
birds[44]<-'Columbina_minuta_FMNH_289160'
birds[48]<-'Chloropipio_holochlora_FMNH_288169'
birds[66]<-'Podica_senegalensis_UMZC_209A_'
birds[75]<-'Leptotila_rufazilla_FMNH_318659'
birds[84]<-'Myiobius_varbatus_FMNH_386812'
birds[105]<-'Poecile_atricapillis_FMNH_504323'
birds[107]<-'Probosciger_aterriumus_NHMUK_S2006.15.5'
birds[117]<-'Rhamphastos_ambiguus_NHMUK_S2002.21'
birds[124]<-'Rupicola_peruviana_UMMZ_119210'

# We have performed the rotation to align the skeletons. We have rejected birds that do not appear in the Meta data 
# All of the rotated constellations host 149 birds


get.item <- function( X , item ) { X[[ item ]] }	
GPA.coords <- lapply( GPA.results , get.item , item = "coords" )
GPA.Csize <- lapply( GPA.results , get.item , item = "Csize" )

# We have organised the data into a single object

tree_genera <- as.character( metadata[ birds , "Prum_genus" ] )
tree_names <- PrumTree$tip.label[ match( tree_genera , sub( "_.*" , "\\1" , PrumTree$tip.label ) ) ]
pruned.tree <- drop.tip( PrumTree , PrumTree$tip.label[ !PrumTree$tip.label %in% tree_names ] )

masses <- apply( metadata[ birds , c( "Mass_F_.hbw_alive." , "Mass_M_.hbw_alive." ) ] , 1 , mean )
#names( masses ) <- tree_names

# Now it is time to save the data
setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')

save(masses, file='masses.22.10.2022.RData')
save(pruned.tree,file='tree.22.10.2022.RData')
save(GPA.Csize,file='Csize.22.10.2022.RData')
save(GPA.coords,file='coords.22.10.2022.RData')
save(tree_names,file='tree.names.22.10.2022.RData')

# That should be all the data that is needed for the moment. 

# What about a version of the shape data that is corrected for allometry?

allometry_list <- list()
allometry.residuals <- list()

library(abind)

for(i in 1:length(GPA.coords)){
	shape.temp<-GPA.coords[[i]]
	dimnames(shape.temp)[[3]]<-tree_names
	gdf.temp<-geomorph.data.frame(shape=shape.temp,size=log10(masses),phy=pruned.tree)
	allometry_list[[i]]<-procD.pgls(shape~size,phy=pruned.tree,data=gdf.temp)
	names(allometry_list)[[i]]<-elements[i]
	for(j in 1:nrow(allometry_list[[i]]$pgls.residuals)){
		if(j==1){
			allometry.residuals[[i]]<-matrix(allometry_list[[elements[i]]]$pgls.residuals[j,],ncol=3,byrow=T)+matrix(t(GPA.results[[elements[i]]]$consensus),ncol=3,byrow=T)
			rownames(allometry.residuals[[i]])<-dimnames(shape.temp)[[1]]
			colnames(allometry.residuals[[i]])<-c('x','y','z')
		} else {
		allometry.residuals[[i]]<-abind(allometry.residuals[[i]],matrix(allometry_list[[elements[i]]]$pgls.residuals[j,],ncol=3,byrow=T)+matrix(t(GPA.results[[elements[i]]]$consensus),ncol=3,byrow=T),along=3)
 		}
	}
}
names(allometry.residuals)<-names(allometry_list)<-names(GPA.coords)

for(i in 1:dim(GPA.coords[['skull']])[3]){
plot3d(GPA.coords[['skull']][,,i],aspect='iso',col='red')
Sys.sleep(2)
plot3d(allometry.residuals[['skull']][,,i],aspect='iso',col='blue')
Sys.sleep(2)
}

for(i in 1:length(GPA.coords)){
	dimnames(allometry.residuals[[i]])<-dimnames(GPA.coords[[i]])
}

save(allometry.residuals,file='allom_coords.29.10.2022.RData')


