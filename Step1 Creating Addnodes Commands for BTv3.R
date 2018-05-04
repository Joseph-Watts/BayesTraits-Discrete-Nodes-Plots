#' Creating figures with ancestral state estimations and phlogenetic uncertanty plotted at the nodes in BayesTraitsV3 Discrete (Step 1 of 2)
#' 
#' This is the first of two scripts needed.
#' The purpose of this script is to create a txt file with Addnodes commands for each node on a consensus tree. These commands can then be entered into BayestraitsV3.
#' 
#' The second step, provided in another script, processes the outputs from BayesTraits and creates a figure from them. 
#' 
#' For more information, please see: https://github.com/Joseph-Watts/BayesTraits-Discrete-Nodes-Plots
#'---------------------------------------------

rm(list = ls())
setwd("D:/Code to Share/BayesTraits-Discrete-Nodes-Plots/Demo")
library(ape)

# Read in a consensus tree. The nodes and taxa on this tree must be ordered how they will appear in the final figure. FigTree provides an easy way to view and order nodes. 
cons <- read.nexus("Primates_Consensus.trees")

# Tests to see whether the tree is binary, should be true.
is.binary.tree(cons)

# Create a list of all subtrees of the cons tree. There is a tree for each node.
subt = subtrees(cons) 

# Create a matrix which will serve as the input commands for BayesTraits and numbered node names
BT = matrix(nrow = (cons$Nnode * 2), ncol = (2 + (length(cons$tip.label))))

# Include the addnode command at the start of each row and 
for(i in 1:cons$Nnode){
  BT[i*2-1,1] = paste0("AddTag ","T",i)
  BT[i*2,1] = paste0("AddNode ",i," T",i)
}

# For each tree in the subt list this loops enters the descendant taxa into a row of BT. Each row in this matrix serves to reconstruct a single node on the tree. 
for(j in 1:cons$Nnode)  {
  for(i in 1: length(subt[[j]]$tip.label)){
    BT[j*2-1,i+1] = subt[[j]]$tip.label[i]
  }
}

# Write commands for BayesTraits. This requires some additional lines to create a command file that can be used for BayesTraits.
write.table(BT, file = "Primates_AddNodes.txt", quote = F, sep = "\t", na = "", col.names = F, row.names = F)

