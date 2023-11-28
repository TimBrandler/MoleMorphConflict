# MoleMorph II
# This script is used to simulate gene trees along species trees
rm(list = ls())  # clean local environment
setwd("C:/Users/timbr/Paleobiology/Gene_Morpho/")                # path     

#install.packages("treeducken")
library(treeducken)

#set a seed for reproducibility
set.seed(1123581321)


# read in 
sTrees <- list()

for (i in 1:50){             # this is obviously not clean but solur´tions come later
  t <- paste0("data/speciestrees/species_tree",i,".nex")
  sTrees[[i]] <-read.nexus(t)
}

# Set parameters
N <- 10000
mut_r <- 3
gen_t <- 1e-6
n_sampled <- 1 


# Simulate Gene trees along species phylogeny
gTrees <- list()

for (i in 1:length(sTrees)){
gTrees[[i]] <- sim_msc(sTrees[[i]], 
                    ne = N,
                    mutation_rate = mut_r,
                    generation_time = gen_t,
                    num_sampled_individuals = n_sampled,
                    num_genes = n_sampled)
}


# Plotting
plot(gTrees[[36]][[1]][["container.tree"]])  # same as sTrees[[36]]
plot(gTrees[[36]][[1]][["gene.trees"]][[1]])

# Print to .nex file
for (i in 1:length(gTrees)) {
  write.tree(gTrees[[i]][[1]][["gene.trees"]][[1]],  # Q: One [1] is number of locus, what is the other [1?]
              file = paste0("data/genetrees/gene_tree",i,".nwk") )
}
