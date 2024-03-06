# Molemorph 1
# SimPhy import
# This replaces the scripts I & II (simulation of species and gene trees) 
# as long as they dont work properly (potentially forever)
# These Simulations have now been outsourced to SimPhy (https://github.com/adamallo/SimPhy)
# This Script exclusively deals with the import of the simulations from SimPhy and is totally redundant if you have
# species and gene trees that are according to the documentation preferred by read.nexus() 

rm(list = ls())
direc <- "C:/Users/timbr/Paleobiology/Gene_Morpho/"
setwd(paste0(direc,"PhySim/")) # your directory

# libraries
library(ape)

# Reading of the SimPhy species trees. This basically just moves them, removes a level folders and changes newick to nexus

tre_nums <- sprintf("%02d", 1:50)
tre_es <- paste0(tre_nums, "/s_tree.trees")

sTrees <- list()
for (i in 1:length(tre_nums)){
  
  sTrees[[i]] <- read.tree(file = tre_es[i])
  
}

for (i in 1:length(sTrees)) {
  write.nexus(sTrees[i], file = paste0(direc,"data/speciestrees/species_tree",i,".nex") )
}

###################################################################################################

# same as above just with the gene trees

tre_nums <- sprintf("%02d", 1:50)
tre_es <- paste0(tre_nums, "/g_trees1.trees")

gTrees <- list()
for (i in 1:length(tre_nums)){
  
  gTrees[[i]] <- read.tree(file = tre_es[i])
  
}

for (i in 1:length(gTrees)) {
  write.nexus(gTrees[i], file = paste0("C:/Users/timbr/Paleobiology/Gene_Morpho/data/genetrees/gene_tree",i,".nex") )
}
