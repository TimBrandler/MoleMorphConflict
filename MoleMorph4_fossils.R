# MoleMorph IV
# This script is used to simulate fossils along a phylogeny
# requirements: A species tree of class "phylo"
rm(list = ls())  # clean local environment
setwd("C:/Users/timbr/Paleobiology/Gene_Morpho/")                # path  

#install.packages("devtools","ape")
library(ape)
library(devtools)
#install_github("fossilsim/fossilsim")
library(FossilSim)
library(dplyr)

# load in trees
sTrees <- list()

for (i in 1:50) {
  t <- paste0("data/speciestrees/species_tree",i,".nex")
  sTrees[[i]] <- read.nexus(t)
}

# Set seed for reproducibility
set.seed(6180)

# Set up Parameters
sample_r <- 0.5
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

# Print trees to file
for (i in 1:length(fTrees)) {
  write.nexus(fTrees[i], file = paste0("data/fossiltrees/fossil_tree",i,".nex") )
}


# Extract and export fossil ages to file

dates <- list()

for (i in 1:length(fossil)){

  d <- matrix(data = round(c(fossil[[i]][[1]],fossil[[i]][[3]]), 2), nrow = length(fossil[[i]][[3]]), ncol = 2, dimnames = list(NULL, c("taxon", "age")))
  d <- as.data.frame(d)
  dates[[i]] <- d %>%
    group_by(taxon) %>%
    filter(age == max(age)) %>%
    ungroup()

  cat(write.table(dates[[i]], paste0("data/dates/fossils_",i,".dat"), sep = "\t", row.names = F, col.names = F ))
  }


# Plotting 
plot(fTrees[[1]])

## Shiny App
#library(FossilSimShiny)
#launchFossilSimShiny()
