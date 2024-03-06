# MoleMorph III
# NOTE: This script is inactive due to treeducken not being maintained!!!
# Continue with MoleMorph1i.SimPhy_import
# This script is used to simulate gene trees along species trees which can either be pure extant or +extinct but do not allow for sampled ancestors
# requirements: A species tree object of class "phylo"
rm(list = ls())  # clean local environment
setwd("C:/Users/timbr/Paleobiology/Gene_Morpho/")                # path     

#install.packages("treeducken")
library(treeducken)
library(devtools)

#set a seed for reproducibility
set.seed(162)


# read in 
sTrees <- list()

for (i in 1:50){             
  t <- paste0("data/speciestrees/species_tree",i,".nex")
  sTrees[[i]] <-read.nexus(t)
}

# Set parameters
lambda <- 0.012
mu <- 0.005
trans_r <- 0
n_loci <- 1
trans_t <- "cladewise" # options are "random" & "cladewise", if lateral gene transfer = 0 it is ignored
Ne <- 1000
mut_r <- 3e-3
gen_t <- 1e-6
n_sampled <- 1 

test <- sim_ltBD(sTrees[[1]],lambda, mu,trans_r,n_loci, trans_t)
plot(test[[1]])
plot(sTrees[[2]])

# Simulate locus trees along species phylogeny
lTrees <- list()

for (i in 1:length(sTrees)){
  lTrees[[i]] <- sim_ltBD(sTrees[[i]],
                          lambda,
                          mu,
                          trans_r,
                          n_loci,
                          trans_t)
  
}

# Simulate Gene trees along species phylogeny
#gTrees <- list()

#for (i in 1:length(lTrees)){
#  gTrees[[i]] <- sim_msc(lTrees[[i]][[1]], 
#                    effective_pop_size = Ne,
#                   mutation_rate = mut_r,
#                    generation_time = gen_t,
#                    num_reps = n_sampled)
#}

gTrees <- list()

for (i in 1:length(lTrees)){
gTrees[[i]] <- drop_extinct(lTrees[[i]][[1]])
}

# Print to .nex file
for (i in 1:length(gTrees)) {
  tryCatch(
    {
      write.tree(gTrees[[i]], file = paste0("data/genetrees/gene_tree", i, ".nwk"))
    },
    error = function(e) {
      # Print the error message (optional)
      cat(i,"Error:", conditionMessage(e), "\n")
    }
  )
}


##
# Simulate Gene trees along species phylogeny
gTrees <- list()

for (i in 1:length(sTrees)){
  gTrees[[i]] <- sim_msc(sTrees[[i]], 
                    ne = Ne,
                   mutation_rate = mut_r,
                    generation_time = gen_t,
                    num_sampled_individuals = n_sampled,
                   num_genes = 2)
}

# Print to .nex file
for (i in 1:length(gTrees)) {
  write.tree(gTrees[[i]][[1]][["gene.trees"]][[1]],  # Q: One [1] is number of locus, what is the other [1?]
              file = paste0("data/genetrees/gene_tree",i,".nwk") )
}

# Plotting
plot(gTrees[[32]][[1]][["container.tree"]])  # same as sTrees[[i]] so same as the species/fossil Tree
plot(gTrees[[32]][[1]][["gene.trees"]][[1]])
par(mfrow=c(2,1))



