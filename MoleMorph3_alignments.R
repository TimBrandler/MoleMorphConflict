# MoleMorph III
# As long as script I and II are not working and being substituted by script Ii this effectively is the second script
# This script is used to simulate molecular alignements for gene trees
# requirements: Genetree(s) of class "phylo", nexus format preferred in this script, other formats should work but you need  to change the code for that
rm(list = ls())  # clean local environment
direc <- "C:/Users/timbr/Paleobiology/Gene_Morpho/"  # your path to directory
setwd(direc)                     

#install.packages("phyclust", "phangorn")
library(phyclust)
library(phangorn)

# read in 

gTrees <- list()

for (i in 1:50) {
  t <- paste0("data/genetrees/gene_tree", i, ".nex")
  
  if (file.exists(t)) {           # this check was implemented when I was playing with different simulation Setups
    gTrees[[i]] <- read.nexus(t)  # currently it should not be possible to have extra extinctions in the gene trees, therefore this check is unnecessary but it also does not hurt             
  } else {
    cat("File", t, "does not exist.\n")
  }
}
# the scripts continues below this inactivated section (l.77)
############################################################################################################################################################
## same as above: this was a check for when gene death was set to be >0; currently useless 
#for (i in 1:length(gTrees)) {       
#  tryCatch(
#    {
#      gTrees[[i]][["node.label"]] <- NULL
#    },
#    error = function(e) {
#      cat("Error at index", i, ":", conditionMessage(e), "\n")
#    }
#  )
#}

### The following was also made for a different setup where gene birth was >0: not currently in use
## gene birth generates orthologs in the target locus of present genes, therefore species are left with a non-uniform number of genes
## These Orthologs are what is identified by the second integer 
## The following filters the orthologs 
#tax <- list()
#strings <- list()
#windfall <- list()
#for (j in 1:length(gTrees)) {
#  tax[[j]] <- gTrees[[j]][["tip.label"]]
#
#  strings[[j]] <- list()
#  windfall[[j]] <- list()
#  
#  for (i in 0:100) {
#    strings[[j]][[1+i]] <- tax[[j]][grep(paste0("_", i, "_0$"), tax[[j]])]
#  
#    if (length(strings[[j]][[1 + i]]) == 0) {
#      break
#    }
#    windfall[[j]][[1+i]] <- tax[[j]][-grep(paste0("_", i, "_0$"), tax[[j]])]
#  }
#} 

## Clean empty entries
#for (j in 1:length(gTrees)) {
#  if (length(strings[[j]][[length(strings[[j]])]]) == 0) {
#    strings[[j]] <- strings[[j]][-length(strings[[j]])]
#    windfall[[j]] <- windfall[[j]][-length(windfall[[j]])]
#  }
#}

## Clean entries with only one entry
#for (j in 1:length(gTrees)) {
#  to_remove_indices <- sapply(strings[[j]], function(x) length(x) %in% c(0, 1))
#  strings[[j]] <- strings[[j]][!to_remove_indices]
#}
##########################################################################################################################################################################


# Alignment Sim Setup
# Code below specifies the model of molecular evolution for in-depth info see https://github.com/rambaut/Seq-Gen/blob/master/manual.md
model <- "-mHKY -t1.2 -a0.25 -g4 -l2500 -f0.3,0.2,0.2,0.3 -or"  
# -m controls model, -t controls transistion/transversion which is not equal to but 1/2*kappa, 
# -a controls gamma rate parameter shape, -g controls rate category count 
# -l controls sequence length, -f controls nucleotide frequency, -o controls output format


# print to file

for (i in 1:length(gTrees)){
  t <- gTrees[[i]]
  seqgen(opts = model, newick.tree = write.tree(t),temp.file = paste0("data/sequences/SimSeq",i,".phylip")) 
}





