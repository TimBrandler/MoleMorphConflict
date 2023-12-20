# MoleMorph III
# This script is used to simulate molecular alignements for gene trees
# requirements: A genetree(s) of class "phylo"
rm(list = ls())  # clean local environment
setwd("C:/Users/timbr/Paleobiology/Gene_Morpho/")                # path     

#install.packages("phyclust")
library(phyclust)

# set seed
set.seed(314159)

# read in 
gTrees <- list()

for (i in 1:50){             # this is obviously not clean but solutions come later
  t <- paste0("data/genetrees/gene_tree",i,".nwk")
  gTrees[[i]] <-read.tree(t)
}

# Alignment Sim Setup
model <- "-mHKY -t0.6 -a0.25 -g5 -l2500"  # -m controls model, -t controls transistion/transversion, -a controls


for (i in 1:length(gTrees)){
  t <- gTrees[[i]]
  seqgen(opts = model, newick.tree = write.tree(t),temp.file = paste0("data/sequences/SimSeq",i,".phylip")) # this directly exports the alignments to your folder
}


# Or this
#mean.rate = rgamma(1, shape = 2, scale = 1/2)
#lstdev <- 0.1

#simulated.rates <- rlnorm(internal_branches, log(mean.rate), lstdev)

# generate new branch lengths
#new.blens <- tree$edge.length * simulated.rates /100

# assign them to the tree
#new.tree <- tree
#new.tree$edge.length <- new.blens

#seqs <- seqgen(opts = paste(model, " -l", alignment.length, sep = ""), newick.tree = write.tree(new.tree))[-1]



