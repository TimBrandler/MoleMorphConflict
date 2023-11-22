## Hello again 
rm(list = ls())  # clean local environment
setwd("C:/Users/timbr/Paleobiology/Gene_Morpho/")                      # path


############# Get "Phybase" #####################
#require(devtools)                                        
#devtools::install_github("lliu1871/phybase")               # do this or the one below # or neither
#devtools::install("PhyloProj/R.packages/phybase-master")  # this requires to download the phybase.zip manually
#library(phybase)                                           
#devtools::install_github("fossilsim/fossilsim")
library(FossilSim)
library(treeducken)
library(phyclust)


#?sim_sptree_bdp  # sbr, sdr, numbsim, n_tips, (gsa_stop_mult = 10??

# Set_seed for reproducibility
set.seed(50)

# Set Parameters for tree simulation

species_lambda <- 1              #
species_mu <- 0.2
No_of_trees <- 50
No_of_tips <- 20

# Simulate Trees under Birth-Death process model

tree_list <- sim_stBD(sbr = species_lambda,                # sim_stBD is the better of the two
                      sdr = species_mu,
                      numbsim = No_of_trees,
                      n_tips = No_of_tips)

# Plot trees (Extant + Extinct)
#par(mfrow=c(1,1))
#x11()

#for (i in 1:length(tree_list)) {
  plot(tree_list[[i]])
}


# Drop Extinct lineages - if dropped now the pipeline will
# -have the assumption that all fossils are sampled ancestors of extant lineages
#################################################################
                                                                #
exta_trees <- list()                                            #
                                                                #
for (i in 1:length(tree_list)) {                                #
  exta_trees[[i]] <- drop_extinct(tree_list[[i]])               #
}                                                               #
                                                                #
# Plot trees (Extant)                                           #
                                                                #
for (i in 1:length(exta_trees)) {                               #
  plot(exta_trees[[i]])                                        #
}                                                               #
#################################################################


### Simulate locus trees

#?sim_ltBD
# gbr  - gene birth rate
# gdr  - gene death rate 
# lgtr - gene transfer rate 
# num_loci - number of locus trees to simulate 
# transfer_type - "cladewise" or "random" - 
# random -> same odds for every contemporaneous lineage
# cladewise -> higher odds the more closely related
 
#par(mfrow=c(1,1))

# to do: check logical/interesting rate relationships 
# ?help uses lower than species tree related rates but within the same order of magnitude
# lower rates make trees more uniform (duuuuh)


gene_lambda <- 0.3
gene_mu <- 0.1
gene_gamma <- 0.1
no_loc <- 5

loc_tree <- list()

# Option A simulate loci trees along a pure extant phylogeny ####
# !! IT DOES NOT WORK TO SIMULATE GENE TREES ALONG THIS!!!      #
for (i in 1:length(exta_trees)){                                #
                                                                #
loc_tree[[i]] <- sim_ltBD(species_tree = exta_trees[[i]],
                     gbr = gene_lambda,
                     gdr = gene_mu,
                     lgtr = gene_gamma,
                     num_loci = no_loc,
                     transfer_type = "cladewise")

}

# option B simulate loci trees along a tree with extinct lineages #
                                                                  #
for (i in 1:length(tree_list)){                                   #
                                                                  #
  loc_tree[[i]] <- sim_ltBD(species_tree = tree_list[[i]],        #
                            gbr = gene_lambda,                    #
                            gdr = gene_mu,                        #
                            lgtr = gene_gamma,                    #
                            num_loci = no_loc,                    #
                            transfer_type = "cladewise")          #
                                                                  #
}                                                                 #
###################################################################

#Plotting (might cause session to abort)
#plot(loc_tree[[1]][[1]])


### Simulate coalescent gene trees 
#effective_pop_size	- the effective population size
#generation_time	- unit time per generation (default 1 year per generation)
#mutation_rate	- number of mutations per unit time
#num_reps	- number of coalescent simulations per locus
?sim_mlc

N <-100 #exp(6)
mut_r <- 3
gen_t <- 1
num_reps <- 5

gene_loc <- sim_mlc(loc_tree[[1]][[1]],
                              effective_pop_size = N,
                    generation_time = gen_t,
                    mutation_rate = mut_r,
                              num_reps = 5)
plot(loc_tree[[1]][[1]])
plot(gene_loc$parent_tree[[3]])

gene_tree_list <- list()



for (i in 1:length(loc_tree)){

  for (j in 1:length(loc_tree[[i]])) {
    
                 temp <- sim_mlc(loc_tree[[i]][[j]],
                         effective_pop_size = N,
                         num_reps = num_reps)
         gene_tree_list <- c(gene_tree_list, temp)
                
  }
  }


counter <- ((i)/length((loc_tree)))*100
print(counter)





#plot(gene_loc_tree[["parent_tree"]][[1]])
#par(mfrow = c(2,1))
#plot(tree_list[[1]])
#plot(ext_tree1)
#plot(loc_tree[[2]])
#pdf("gene_trees1.pdf")
#plot(gene_loc_tree[["parent_tree"]][[1]])
#plot(gene_loc_tree[["parent_tree"]][[2]])




# new.tree = your tree

model <- "-mHKY -t5 -a0.25 -g5" 
alignment.length <- 300

seqs41 <- phyclust::seqgen(opts = paste(model, " -l", alignment.length, sep = ""), newick.tree = ape::write.tree(genetreext4[[1]][["gene.trees"]][[1]]))[-1] 

View(seqs41)
seqs41[-1]

# the following will simulate sequences under a relaxed clock model

# sample a mean substitution rate for each locus
# under the gamma distribution
mean.rate = rgamma(1, shape = 2, scale = 1/2)

# simulate variable rates under a log normal distribution
lstdev <- 0.1
simulated.rates <- rlnorm(internal_branches, log(mean.rate), lstdev)



# generate new branch lengths
new.blens <- tree$edge.length * simulated.rates /100

# assign them to the tree
new.tree <- tree
new.tree$edge.length <- new.blens

seqs <- phyclust::seqgen(opts = paste(model, " -l", alignment.length, sep = ""), newick.tree = ape::write.tree(new.tree))[-1]

# the following will simulate fossils
sampl_rate <- 1.5 # fossil sampling rate
fossils <- FossilSim::sim.fossils.poisson(sampl_rate, tree = tree_list[[4]])

# the following is for simulating morpho data for extant and fossil samples

# combines the tree and fossils data frame objects
ftree <- FossilSim::SAtree.from.fossils(tree_list[[4]],fossils)
tree <- FossilSim::sampled.tree.from.combined(ftree[[1]], rho = 0)
plot(tree)
# q matrix for a binary character
?geiger
q <- list(rbind(c(-.5, .5), c(.5, -.5)))
morpho1 <- geiger::sim.char(tree_list[[1]], q, model="discrete") #[,,1]
morpho2 <- geiger::sim.char(tree_list[[2]], q, model="discrete") #[,,1]
morpho3 <- geiger::sim.char(tree_list[[3]], q, model="discrete") #[,,1]
morpho4 <- geiger::sim.char(tree_list[[4]], q, model="discrete") #[,,1]

morpho1[,,1]
morpho2[,,1]
morpho3[,,1]
morpho4[,,1]

plot(tree_list[[1]])
plot(tree_list[[2]])
plot(tree_list[[3]])
plot(tree_list[[4]])



library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

N <- 10 # effective population size  

# fill up the probility transition matrix
P <- matrix(NA, ncol = 2*N + 1, 2*N + 1)
P_df <- data.frame()
for(i in 0:(2*N)){
  for(j in 0:(2*N)){
    # add to matrix for examples
    P[i+1, j+1] <- dbinom(j, 2*N, (i / (2*N))) 
    
    # add to data.frame for viz below
    P_df <- bind_rows(P_df, data.frame(i = i, j = j, p = P[i+1, j+1]))
  }
}

# plot it up!
p <- ggplot(P_df, aes(x = i, y = j, fill = p)) +   
  geom_tile() + scale_fill_viridis(option="D")
p



#### GENE TREE BASED ON SPECIES TREE  ###################
# species tree - phylo-obj
# ne - pop size
# num_sampled_individuals - within each lineage
# num_genes - to simulate within each locus
# rescale - TRUE/FALSE rescale tree to coalescent units
# mutation_rate - per generation
# generation_time - No of time units per generation 
#?sim_msc

Pop <- 1000000
mut_r <- 1e-9
gen_t <- 1e-6
n_sampled <- 3
n_genes <- 5

N <-10 #exp(6)
mut_r <- exp(-6)
gen_t <- 25


x11()
genetree <- sim_msc(ext_tree4, 
                    ne = Pop,
                    mutation_rate = mut_r,
                    generation_time = gen_t,
                    num_sampled_individuals = n_sampled,
                    num_genes = n_genes)

genetreext4 <- sim_msc(tree_list[[4]], 
                       ne = Pop,
                       mutation_rate = mut_r,
                       generation_time = gen_t,
                       num_sampled_individuals = n_sampled,
                       num_genes = n_genes)


plot(tree_list[[1]])
plot(ext_tree1)
plot(genetree[[1]][["gene.trees"]][[1]])
plot(genetreexta[[1]][["gene.trees"]][[2]])