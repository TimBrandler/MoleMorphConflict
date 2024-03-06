# MoleMorph I
# This script is used to simulate species trees under a BD - process
# NOTE: This script is inactive due to treeducken not being maintained !!
# Start with MoleMorph1i_SimPhy_import !!

rm(list = ls())  # clean local environment
setwd("C:/Users/timbr/Paleobiology/Gene_Morpho/")                # path     

#install.packages("treeducken", "taxize")
library(treeducken)

#set a seed for reproducibility
set.seed(1123581321)


# Set BD paramaters
# Carnivora set based on Lee Hsiang Liow & John A. Finarelli, 2014:
# lambda <- 0.244 #+/- 0.066
# mu <- 0.197 # +/- 0.06

# 

lambda <- 0.25 # Speciation
mu <- 0.0 # Extinction
Tips <- 20 # Size # 20 gets quite extensive with l/m of ~1.2 (up to ~100 Tips) 
Sims <- 50 # Different trees / Simulations

# Simulate using the treeducken package
tree_list <- sim_stBD(sbr = lambda,               
                      sdr = mu,
                      numbsim = Sims,
                      n_tips = Tips)


for (i in 1:length(tree_list)) {
  if (length(tree_list[[i]]) < 5) {
    mt <- print(i)
    tree_list[i] <- tree_list[-i]
  }
}

num_entries <- 300

# Generate the sequences
sequence_generator <- function(n) {
  paste0(letters[(n - 1) %% 26 + 1], letters[((n - 1) %/% 26) %% 26 + 1], letters[((n - 1) %/% (26^2)) %% 26 + 1])
}

shorts <- data.frame(ID = sapply(seq_len(num_entries), sequence_generator))

# Get the number of tips in the tree

for (i in 1:length(tree_list)) {
  num_tips <- length(tree_list[[i]]$tip.label)
  # Assign the generated tip labels to the tree
  tree_list[[i]][["tip.label"]] <- shorts[1:num_tips,]
}


#for (i in 1:length(tree_list)){
#  tree_list[[i]][["tip.label"]] <- seq(1:length(tree_list[[i]][["tip.label"]]))
#}

# Print trees to files
for (i in 1:length(tree_list)) {
  write.nexus(tree_list[i], file = paste0("data/speciestrees/species_tree",i,".nex") )
}

# Plot the trees (if u feel like it)
for (i in 1:length(tree_list)) {
  plot(tree_list[[i]])
}


temp <- sim.bd.taxa(10, 50, 0.25, 0)
tree_list <- list()
for (i in 1:50){
  write.tree(temp[[i]], file = paste0("data/test/species_tree",i,".nwk") )
}
