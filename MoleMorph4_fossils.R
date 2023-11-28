# MoleMorph IV
# This script is used to simulate fossils along a phylogeny
rm(list = ls())  # clean local environment
setwd("C:/Users/timbr/Paleobiology/Gene_Morpho/")                # path  

#install.packages("devtools","ape")
library(ape)
library(devtools)
#install_github("fossilsim/fossilsim")
library(FossilSim)

# load in trees
sTrees <- list()

for (i in 1:50) {
  t <- paste0("data/speciestrees/species_tree",i,".nex")
  sTrees[[i]] <- read.nexus(t)
}

# Set_seed for reproducibility
set.seed(6180)

# Set up Parameters
sample_r <- 1.5
rho <- 1

# Simulate fossils on a tree
fossil <- list()
for (i in 1:length(sTrees)) {
  fossil[[i]] <- sim.fossils.poisson(sample_r, tree = sTrees[[i]])
}

# Combine Phylogenies & Fossils
fTrees <- list()

for (i in 1:length(fossil)){
  tempFtree <- SAtree.from.fossils(sTrees[[i]],fossil[[i]])
  fTrees[[i]]<- sampled.tree.from.combined(tempFtree[[1]], rho = rho)
}

# Plotting 
plot(fTrees[[1]])

# Print trees to file
for (i in 1:length(fTrees)) {
  write.nexus(fTrees[i], file = paste0("data/fossiltrees/fossil_tree",i,".nex") )
}



