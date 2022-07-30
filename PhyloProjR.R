## Hello again 
rm(list = ls())  # clean local environment !OPTIONAL!
setwd("C:/Users/timbr/Paleobiology")                      # set your own path!


############# Get "Phybase" #####################
require(devtools)                                          # Why did I ever bother using library()?
#devtools::install_github("lliu1871/phybase")               # do this or the one below # or neither
#devtools::install("PhyloProj/R.packages/phybase-master")  # this requires to download the phybase.zip manually
library(phybase)                                           # xD
devtools::install_github("fossilsim/fossilsim")
library(FossilSim)
require(treeducken)
require(phyclust)


#?sim_sptree_bdp  # sbr, sdr, numbsim, n_tips, (gsa_stop_mult = 10???????????????)

# Set_seed for reproducibility
set.seed(50)

# Set Parameters for tree simulation

species_birth_rate <- 4              # Think about all of them again
species_death_rate <- 2.5
No_of_trees <- 4
No_of_tips <- 20 # 20?

# Simulate Trees under Birth-Death process model
?sim_stBD

tree_list <- sim_stBD(sbr = species_birth_rate,                # sim_stBD is the better of the two
                      sdr = species_death_rate,
                      numbsim = No_of_trees,
                      n_tips = No_of_tips)

# Plot trees (Extant + Extinct)

par(mfrow = c(2,2))                                  ## There...........
plot(tree_list[[1]])
plot(tree_list[[2]])
plot(tree_list[[3]])
plot(tree_list[[4]])
par(mfrow = c(1,1))                                  ## ..... and back again :)

# Drop Extinct lineages

ext_tree1 <- drop_extinct(tree_list[[1]])
ext_tree2 <- drop_extinct(tree_list[[2]])
ext_tree3 <- drop_extinct(tree_list[[3]])
ext_tree4 <- drop_extinct(tree_list[[4]])

# Plot trees (Extant)

par(mfrow = c(2,2))                           ## There............
plot(ext_tree1)
plot(ext_tree2)
plot(ext_tree3)
plot(ext_tree4)
par(mfrow = c(1,1))                           ## ........ and back again (JRR Tolkien btw:)

######## this is the start of the beginning so far ....

### Simulate locus trees

#?sim_ltBD
# gbr  - gene birth rate
# gdr  - gene death rate 
# lgtr - gene transfer rate 
# num_loci - number of locus trees to simulate 
# transfer_type - "cladewise" or "random"

gene_b <- 0.1
gene_d <- 0.04
transfer_r <- 0.02
no_loc <- 4

loc_tree <- sim_ltBD(species_tree = ext_tree1,
         gbr = gene_b,
         gdr = gene_d,
         lgtr = transfer_r,
         num_loci = no_loc,
         transfer_type = "random")

plot(loc_tree[[2]])


### Simulate coalescent gene trees 
?sim_mlc
#locus_tree - a locus tree from 'sim_ltBD' of class 'phy'
#effective_pop_size	- the effective population size

#generation_time	- unit time per generation (default 1 year per generation)

#mutation_rate	- number of mutations per unit time

#num_reps	- number of coalescent simulations per locus

gene_loc_tree <- sim_mlc(loc_tree[[2]],
                         effective_pop_size = 1e6,
                         num_reps = 20)

plot(gene_loc_tree[["parent_tree"]][[1]])
par(mfrow = c(2,1))
plot(tree_list[[1]])
plot(ext_tree1)
plot(loc_tree[[2]])
pdf("gene_trees1.pdf")
plot(gene_loc_tree[["parent_tree"]][[1]])
plot(gene_loc_tree[["parent_tree"]][[2]])
dev.off()


#### GENE TREE BASED ON SPECIES TREE  ###################
# species tree - phylo-obj
# ne - pop size
# num_sampled_individuals - within each lineage
# num_genes - to simulate within each locus
# rescale - TRUE/FALSE rescale tree to coalescent units
# mutation_rate - per generation
# generation_time - No of time units per generation 
?sim_msc

Pop <- 1000000
mut_r <- 1e-9
gen_t <- 1e-6
n_sampled <- 3
n_genes <- 5

x11()
genetree <- sim_msc(ext_tree4, 
                    ne = Pop,
                    mutation_rate = mut_r,
                    generation_time = gen_t,
                    num_sampled_individuals = n_sampled,
                    num_genes = n_genes)
plot(genetree[[1]][["gene.trees"]][[1]])

genetreext4 <- sim_msc(tree_list[[4]], 
                    ne = Pop,
                    mutation_rate = mut_r,
                    generation_time = gen_t,
                    num_sampled_individuals = n_sampled,
                    num_genes = n_genes)

plot(genetreext4[[1]][["gene.trees"]][[1]])


plot(tree_list[[1]])
plot(ext_tree1)
plot(genetree[[1]][["gene.trees"]][[1]])
plot(genetreexta[[1]][["gene.trees"]][[2]])

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


