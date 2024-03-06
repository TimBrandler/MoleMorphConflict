# MoleMorph IV
# As long as script I and II are not working and being substituted by script Ii this effectively is the third script
# This script is used to simulate and sample fossils along a phylogeny
# requirements: A species tree of class "phylo"
# Note: depending on which assumptions you are making you need to go to different sections
# 1. Natural: Species and Gene sorting is subject to a stochastic process, start at the beginning
# 2. Meeting MSC-FBD assumptions: Morphology evolves along the gene tree, go to l.66

rm(list = ls())  # clean local environment
direc <- "C:/Users/timbr/Paleobiology/Gene_Morpho/"   # your path
setwd(direc)                  

#install.packages("devtools","ape")
library(ape)
#library(devtools)
#install_github("fossilsim/fossilsim")
library(FossilSim)

# 'Natural' assumptions
# load in trees
sTrees <- list()

for (i in 1:50) { # Note that the "50" is just the number of simulations that were run and no necessity   
  t <- paste0("data/speciestrees/species_tree",i,".nex")
  sTrees[[i]] <- read.nexus(t)
}

# Set up Parameters
sample_r <- 0.2      # fossil sampling rate, we will use this in a poisson sampling process
rho <- 1             # recent sampling, we assume perfect sampling in the the recent

# Simulate fossils on a tree
fossil <- list()

for (i in 1:length(sTrees)) {
  fossil[[i]] <- sim.fossils.poisson(sample_r, tree = sTrees[[i]])
}

# Combine Phylogenies & Fossils
fTrees <- list()
tempFtree <- list()

for (i in 1:length(fossil)){
  tempFtree[[i]] <- SAtree.from.fossils(sTrees[[i]],fossil[[i]])
  
  fTrees[[i]] <- sampled.tree.from.combined(tempFtree[[i]][[1]], rho = rho)
}

# Print trees to file
for (i in 1:length(fTrees)) {
  write.nexus(fTrees[i], file = paste0("data/fossiltrees/fossil_tree",i,".nex") )
}

dates <- list()

for (i in 1:length(fossil)){

  d <- matrix(data = c(tempFtree[[i]]$fossils[,6],round(fossil[[i]][[3]], 2)), nrow = length(fossil[[i]][[3]]), ncol = 2, dimnames = list(NULL, c("taxon", "age")))
  dates[[i]] <- as.data.frame(d)
  dates[[i]][["taxon"]] <- paste0(dates[[i]][["taxon"]])

  cat(write.table(dates[[i]], paste0("data/dates/fossils_",i,"_Stndrd.dat"), sep = "\t", row.names = F, col.names = F, quote = F ))
}


# Morphology evolves based on the gene tree
# This is bascially the same steps as above just with a different starting point and saving location
rm(list = ls())  # clean local environment
direc <- "C:/Users/timbr/Paleobiology/Gene_Morpho/"   # your path
setwd(direc) 

# MSC assumptions
# load in trees
mscTrees <- list()

for (i in 1:50) { # Note that the "50" is just the number of simulations that were run and no necessity   
  t <- paste0("data/genetrees/gene_tree",i,".nex")
  mscTrees[[i]] <- read.nexus(t)
}

# Set up Parameters
sample_r <- 0.2      # fossil sampling rate, we will use this in a poisson sampling process
rho <- 1             # recent sampling, we assume perfect sampling in the the recent

# Simulate fossils on a tree
fozzil <- list()

for (i in 1:length(sTrees)) {
  fozzil[[i]] <- sim.fossils.poisson(sample_r, tree = sTrees[[i]])
}

# Combine Phylogenies & Fossils
mfcTrees <- list()
tempFtree <- list()

for (i in 1:length(fossil)){
  tempFtree[[i]] <- SAtree.from.fossils(sTrees[[i]],fozzil[[i]])
  
  mfcTrees[[i]] <- sampled.tree.from.combined(tempFtree[[i]][[1]], rho = rho)
}

# Print trees to file
for (i in 1:length(mfcTrees)) {
  write.nexus(fTrees[i], file = paste0("data/fossiltrees/mfc_tree",i,".nex") )
}

dates <- list()

for (i in 1:length(fossil)){
  
  d <- matrix(data = c(tempFtree[[i]]$fossils[,6],round(fossil[[i]][[3]], 2)), nrow = length(fozzil[[i]][[3]]), ncol = 2, dimnames = list(NULL, c("taxon", "age")))
  dates[[i]] <- as.data.frame(d)
  dates[[i]][["taxon"]] <- paste0(dates[[i]][["taxon"]])
  
  cat(write.table(dates[[i]], paste0("data/dates/fossils_",i,"_MSC.dat"), sep = "\t", row.names = F, col.names = F, quote = F ))
}
