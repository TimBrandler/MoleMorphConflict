# MoleMorph V
# This script is used to simulate morphological characters along a phylogeny
# A species Tree (with or w/o) sampled ancestors
rm(list = ls())  # clean local environment
setwd("C:/Users/timbr/Paleobiology/Gene_Morpho/")                # path  

#install.packages("geiger", "ouch")
library(geiger)
library(ouch)
#library(devtools)
#install_github("liamrevell/phytools")
# read in trees
fTrees <- list()

for (i in 1:50) {
  t <- paste0("data/fossiltrees/fossil_tree",i,".nex")
  fTrees[[i]] <- read.nexus(t)
}


library(ape)
library(phytools)


# Function for transition matrix with n number of states

Q <- function(s, c){
  temp <-  matrix(1/s, nrow = s, ncol = s) - diag(1, s)
  rep(list(temp), c)
}



# Simulate Morphological Characters
# Loop through the first 20 tip labels and remove the "_integer" part
# This makes it easier to identify the extant species (for you but more importantly during the xml file creation)  

for (j in 1:length(fTrees)) {
  for (i in 1:20) {
    fTrees[[j]][["tip.label"]][i] <- gsub("_\\d+$", "", fTrees[[j]][["tip.label"]][i])
  }
}

# the sim.char function we are using does not have a substitution rate parameter (or effectively fixed at 1)
# so we scale the branch lengths beforehand
scaled <- list()

for (i in 1:length(fTrees)) {
scaled[[i]] <- fTrees[[i]]
scaled[[i]][["edge.length"]] <- fTrees[[i]][["edge.length"]]*0.2

}

## This just creates vectors of characters with different max numbers of states
morpho <- list()

for (h in 1:3) {
  if (h == 1) {
    for (i in 1:length(fTrees)){
      morpho[[i]] <- sim.char(scaled[[i]], par = Q(2, 50), model = "discrete")[,,1]
      }
  } else if (h == 2) {
    for (i in 1:length(fTrees)){
      morpho[[i]] <- cbind(morpho[[i]], sim.char(scaled[[i]], par = Q(3, 150), model = "discrete")[,,1])
    }
  } else if (h == 3) {
   for (i in 1:length(fTrees)){
      morpho[[i]] <- cbind(morpho[[i]], sim.char(scaled[[i]], par = Q(4, 300), model = "discrete")[,,1])
    }
  }
}
# for a more complex setup you can use different scale values for the different partitions (bcs like this we effectively have no partitions)


# print to file 
for (i in 1:length(morpho)) {
  write.nexus.data(morpho[[i]], paste0("data/characters/chars",i,".nex"), format = "standard", datablock = T, interleaved = F)
}


