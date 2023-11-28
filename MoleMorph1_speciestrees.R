# MoleMorph I
# This script is used to simulate species trees under a BD - process
#rm(list = ls())  # clean local environment
setwd("C:/Users/timbr/Paleobiology/Gene_Morpho/")                # path     

#install.packages("treeducken")
library(treeducken)

#set a seed for reproducibility
set.seed(1123581321)


# Set BD paramaters
lambda <- 1 # Speciation
mu <- 0.3 # Extinction
Tips <- 20 # Size
Sims <- 50 # Different trees / Simulations

# Simulate using the treeducken package
tree_list <- sim_stBD(sbr = lambda,               
                      sdr = mu,
                      numbsim = Sims,
                      n_tips = Tips)

# Plot the trees (if u feel like it)
for (i in 1:length(tree_list)) {
  plot(tree_list[[i]])
}

# Print trees to files
for (i in 1:length(tree_list)) {
  write.nexus(tree_list[i], file = paste0("data/speciestrees/species_tree",i,".nex") )
}

