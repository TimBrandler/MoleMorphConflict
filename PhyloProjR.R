## Hello again 
rm(list = ls())  # clean local environment !OPTIONAL!
setwd("C:/Users/timbr/Paleobiology")                      # set your own path!


############# Get "Phybase" #####################
require(devtools)                                          # Why did I ever bother using library()?
devtools::install_github("lliu1871/phybase")               # do this or the one below 
#devtools::install("PhyloProj/R.packages/phybase-master")  # this requires to download the phybase.zip manually
library(phybase)                                           # xD
require(treeducken)


#?sim_sptree_bdp  # sbr, sdr, numbsim, n_tips, (gsa_stop_mult = 10???????????????)

# Set Parameters for tree simulation

species_birth_rate <- 4              # Think about all of them again
species_death_rate <- 2.5
No_of_trees <- 4
No_of_tips <- 20 # 20?

# Simulate Trees under Birth_Death process model

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
