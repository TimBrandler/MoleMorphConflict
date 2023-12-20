# MoleMorph V
# This script is used to simulate morphological characters along a phylogeny
# A species Tree (with or w/o) sampled ancestors
rm(list = ls())  # clean local environment
setwd("C:/Users/timbr/Paleobiology/Gene_Morpho/")                # path  

#install.packages("geiger", "ouch")
library(geiger)
library(ouch)

# read in trees
fTrees <- list()

for (i in 1:50) {
  t <- paste0("data/fossiltrees/fossil_tree",i,".nex")
  fTrees[[i]] <- read.nexus(t)
}

# Set number of characters, number of states
n_chars <- 300
n_states <- 3

# Function for transition matrix with n number of states
Q <- function(s, c){
  temp <-  matrix(1/s, nrow = s, ncol = s) - diag(1, s)
  rep(list(temp), c)
}

#QQ <- 

#q <- rep(list(rbind(c(-2/3, 1/3, 1/3), c(1/3, -2/3, 1/3),c(1/3, 1/3, -2/3))),n_chars)

# Simulate Morphological Characters
morpho <- list()

for ( i in 1:length(fTrees)){
morpho[[i]] <- sim.char(fTrees[[i]], par = Q(n_states, n_chars), model="discrete" )[,,1]
}

# morph <- matrix(, nrow = length(morpho))

#for (i in 1:length(morpho)) {
#  for (j in 1:length(morpho[[i]][,1])) {
#    morph[[i]][j] <- list(morpho[[i]][j,]) 
#  }  
#}

# print to file 

for (i in 1:length(morpho)) {
  write.dna(morpho[[i]], paste0("data/characters/phylo",i,".nex"), format = "fasta", nbcol = -1, colsep = " ", colw = n_chars)
}
