#' Creating figures with ancestral state estimations and phlogenetic uncertanty plotted at the nodes in BayesTraitsV3 Discrete (Step 2 of 2)
#' 
#' This is the second of two scripts needed.
#' The purpose of this script is to read in the outputs from BayesTraits and create two figures, one for each trait. 
#' 
#' For more information, please see: https://github.com/Joseph-Watts/BayesTraits-Discrete-Nodes-Plots
#'---------------------------------------------

rm(list = ls())
setwd("D:/Code to Share/BayesTraits-Discrete-Nodes-Plots/Demo")
library(ape)

# Read in a consensus tree
tree = read.nexus("Primates_Consensus.trees")

# Read BT ouput. Skip all the lines that are part of the header (length will vary depending on settings and taxa).
BT_output = read.table(file = "Primates.txt.Log.txt", sep = "\t", header = T, na = "--", skip = 168)

first = which(colnames(BT_output)=="X1...P.0.0.")
last = (tree$Nnode * 4) + first - 1
BT_output = BT_output[,c(first:last)]

# Reading in taxa names for figure
fig_taxa_names = read.table(file = "Primate_Display_Names_Index.txt", header = T, sep = "\t")
fig_taxa_names = fig_taxa_names[order(fig_taxa_names[,1]), ] # ordering the rows alphabetically by the first column

# Reading in the trait data used in the bayestraits analysis. The first collumn is taxa names, second trait A, and third trait B.
traits = read.table(file = "Primates.txt", header = F, sep = "\t", na.strings = "-" )

row.names(traits) = traits[ ,1]
traits = traits[order(row.names(traits)), ] # ordering the rows in the traits file alphabetically
colnames(traits) = c("char_taxa", "trait_A", "trait_B")

traits = merge(traits, fig_taxa_names, by.x = "char_taxa", by.y = "Tree_Name")

## Ordering traits by the tip.labels in the consensus tree
index = 1:nrow(traits)
tree_taxa = tree$tip.label
tree_taxa_index = data.frame(tree_taxa, index)
row.names(tree_taxa_index) = tree$tip.label
tree_taxa_index = tree_taxa_index[order(row.names(tree_taxa_index)), ]

combined = data.frame(tree_taxa_index, traits)

table(combined$"tree_taxa" == combined$char_taxa)

# -------------------------------------------------------------------------------------------- 
# At this point it needs to be checked that there are no false responses.
# If there are false responses then the taxa in the tree file do not match the taxa in the traits file and this needs to be fixed.
# --------------------------------------------------------------------------------------------  

combined = combined[order(combined$index), ]
traits_ordered = combined[ ,4:6]

# If the row names in the trait dataframe are the same as those in the tree tip lable then change them all to the Society names. 
if(sum(!row.names(traits_ordered) == tree$tip.label) <= 0){
  tree$tip.label = as.character(traits_ordered$Display_Name)
  row.names(traits_ordered) = traits_ordered$Display_Name
  traits_ordered = traits_ordered[,1:2]
}


# Construct two matricies, one for each trait. 
# In each matrix there are 3 collumns, 0, 1, NA. The first column is the liklihood of the trait being in state 0 at the given node, the second collumn is the liklihood of the trait being in state 1, and the third collumn is the percentage of the trees in which that node is absent.
# For trait A (the first (e.g. (A,B)))

st0 = NULL; st1 = NULL; st_NA = NULL
d = 1; k = 2; nodes = tree$Nnode
n = nrow(BT_output)
for(i in 1:nodes){
  d = i * 4 - 3
  av1 = sum(as.numeric(BT_output[,d]), na.rm = TRUE) / n
  av2 = sum(as.numeric(BT_output[,d + 1]), na.rm = TRUE) / n
  st0[i] = av1 + av2
  k = i * 4 - 1
  av1 = sum(BT_output[,k], na.rm = TRUE) / n
  av2 = sum(BT_output[,k + 1], na.rm = TRUE) / n
  st1[i] = av1 + av2
  st_NA[i] = 1 - st0[i] - st1[i]
}
A_p_nodes = cbind(st0, st1, st_NA) 
A_p_nodes = A_p_nodes[1:nrow(A_p_nodes), ] 

#For trait B (the second (e.g. (A,B)))

st0 = NULL; st1 = NULL; st_NA = NULL
d = 1; k = 2; nodes = tree$Nnode
n = nrow(BT_output)
for(i in 1:nodes){
  d = i * 4 - 3
  av1 = sum(as.numeric(BT_output[,d]), na.rm = TRUE) / n
  av2 = sum(as.numeric(BT_output[,d + 2]), na.rm = TRUE) / n
  st0[i] = av1 + av2
  k = i * 4 - 2
  av1 = sum(as.numeric(BT_output[,k]), na.rm = TRUE) / n
  av2 = sum(as.numeric(BT_output[,k + 2]), na.rm = TRUE) / n
  st1[i] = av1 + av2
  st_NA[i] = 1 - st0[i] - st1[i]
}
B_p_nodes = cbind(st0, st1, st_NA)
B_p_nodes = B_p_nodes[1:nrow(B_p_nodes), ] 

# ------------------------
# Plotting trait A on the consensus tree with pie charts at tips and nodes
# ------------------------

# Specifying the colours to use
A_present_colour = "olivedrab"
Absent_colour = "white"
NA_colour = "gray"

Trait_A_colours = ifelse(is.na(traits_ordered$trait_A), NA_colour,
                     ifelse(traits_ordered$trait_A == 0, Absent_colour,
                            ifelse(traits_ordered$trait_A == 1, A_present_colour,NA)
                            )
                     )

pdf(file = "Primates_Trait_A.pdf", width = 6, height = 10, compress = F)

plot.phylo(tree, edge.width = .5, label.offset = 0.01, cex= .5, adj = 0, direction = "rightwards")

tiplabels(pch = 21, bg = Trait_A_colours, adj = .5, cex = 1.0)

# Reads in the csv file that was created at the start, csv file must have been edited before it is read back in.
# A matrix is constructed from the csv file and then pie charts are constructed at the nodes. 

A_p_nodes_2 = ifelse(A_p_nodes<0, 0, A_p_nodes)
nodelabels(pie=A_p_nodes_2, cex= 0.5, edge.width = .3, piecol=c(Absent_colour, A_present_colour, NA_colour))

dev.off()

# ------------------------
# Plotting trait B on the consensus tree with pie charts at tips and nodes
# ------------------------

# Specifying the colours to use
B_present_colour = "brown3"
Absent_colour = "white"
NA_colour = "gray"

Trait_B_colours = ifelse(is.na(traits_ordered$trait_B), NA_colour,
                         ifelse(traits_ordered$trait_B == 0, Absent_colour,
                                ifelse(traits_ordered$trait_B == 1, B_present_colour,NA)
                         )
)

pdf(file = "Primates_Trait_B.pdf", width = 6, height = 10, compress = F)

plot.phylo(tree, edge.width = .5, label.offset = 0.01, cex= .5, adj = 0, direction = "rightwards")

tiplabels(pch = 21, bg = Trait_B_colours, adj = .5, cex = 1.0)

# Reads in the csv file that was created at the start, csv file must have been edited before it is read back in.
# A matrix is constructed from the csv file and then pie charts are constructed at the nodes. 

B_p_nodes_2 = ifelse(B_p_nodes<0, 0, B_p_nodes)
nodelabels(pie=B_p_nodes_2, cex= 0.5, edge.width = .3, piecol=c(Absent_colour, B_present_colour, NA_colour))

dev.off()

# Note: If you change the data, you will need to adjust the sizing values of the commands and the canvas.
