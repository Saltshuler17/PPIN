# Protein-Protein Interaction Networking
# Creating the figure legends components
# By Sam Altshuler
# Originated on: 05/09/20



library(dplyr)
library(stringr)
library(data.table)

# Create a unified legend for all figures that are presented together
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Functions -----
# Function to each gene in the time points to all of the GO terms
# they are associated with (term A -> gene 1, 2, 3), this is collecting
# information already gathered in the network data table, but allows for
# more pruning of the GO terms
map2term <- function(allgenes, TP){
  # create a data frame with all of the GO terms in one column
  # and all of the genes associated with that term in the second
  masterdf <- data.frame("Parent" = TP, stringsAsFactors = F)
  # create these terms that can be used to differentiate the 
  # different ontologies from the preformed network data
  xterm <- paste("GoID", "Source", TP[[1]], sep = "_")
  nterm <- paste("GoID", "Node", TP[[1]], sep = "_")
  lil <- vector(length = 0)
  t <- 1
  # loop through every GO term (parent level == second level) for 
  # this specific ontology (denoted by the first value in the list
  # of GO terms, either BP, MF, or CC)
  for (i in 1:length(TP)){
    masterdf$Genes[[i]] <- vector("list", 1)
    # create a list to store all the genes for each term
    glist <- vector(length = 0)
    c <- 1
    # loop through all genes, if the go term is associated with that
    # gene, the gene is added to the list. Check both source and Node columns
    for (j in 1:length(allgenes$Source)){
      if (TP[[i]] %in% allgenes[[xterm]][[j]]){
        glist[[c]] <- as.character(allgenes$SourceNCU[[j]])
        c <- c + 1
      }
      if (TP[[i]] %in% allgenes[[nterm]][[j]]){
        glist[[c]] <- as.character(allgenes$NodeNCU[[j]])
        c <- c + 1
      }
    }
    # get rid of duplicates if any show up
    glist <- as.list(unique(glist))
    
    # ARBITRARY NUMBER THIS IS VARIABLE, can be 10 or 20 or whatever
    # gets rid of all terms that don't have "enough"
    # genes associated with them across all of the data
    if (length(glist) >= 20) {
      masterdf$Genes[[i]] <- glist
      masterdf$Leng[[i]] <- length(glist)
    } else {
      masterdf$Leng[[i]] <- length(glist)
      lil[[t]] <- i
      t <- t +1
    }
  }
  # delete all terms that have less than the arbitrary # chosen
  masterdf <- masterdf[-c(lil),]
  # return the complete data frame that contains all GO terms (per ontology)
  # and all genes in all time points (+LL) associated with each term
  return(masterdf)
}

# Function to map from the data frame of terms back to the network data frame
# (from Terms A -> gene 1, 2, 3 to gene 1 -> term a, b, c)
remap2gene <- function(gene_frame, masterdf){
  # gene_frame = the network for a specific time point
  # use the name of the data frame with the GO terms to know
  # which ontology you are saving into 
  goname <- deparse(substitute(masterdf))
  xterm <- paste("GoID", "Source", goname, sep = "_")
  nterm <- paste("GoID", "Node", goname, sep = "_")
  # Loop through all interactions in the time point,
  # both source and node genes
  for (k in 1:length(gene_frame$SourceNCU)){
    c <- 1
    h <- 1
    gene_frame[[xterm]][[k]] <- vector("list",1)
    gene_frame[[nterm]][[k]] <- vector("list",1)
    Gsource <- vector(length = 0)
    Gnode <- vector(length = 0)
    # loop through all GO terms, if the GO term contains
    # the gene the loop is on, that term is added to a list
    for (i in 1:length(masterdf$Parent)){
      #print(i)
      if (gene_frame$SourceNCU[[k]] %in% masterdf$Genes[[i]]){
        Gsource[[c]] <- as.character(masterdf$Parent[[i]])
        c <- c + 1
      }
      if (gene_frame$NodeNCU[[k]] %in% masterdf$Genes[[i]]){
        Gnode[[h]] <- as.character(masterdf$Parent[[i]])
        h <- h + 1
      }
    }
    
    # add the GO terms to either the source or the node
    # and if no GO terms were detected, that gene is labeled "Unknown"
    
    Gsource <- as.list(unique(Gsource))
    if (length(Gsource) > 0) {
      gene_frame[[xterm]][[k]] <- Gsource
    } else {
      gene_frame[[xterm]][[k]] <- c("Unknown")
    }
    
    Gnode <- as.list(unique(Gnode))
    if (length(Gnode) > 0) {
      gene_frame[[nterm]][[k]] <- Gnode
    } else {
      gene_frame[[nterm]][[k]] <- c("Unknown")
    }
  }
  # return the updated data frame for each timepoint
  return(gene_frame)
}

# Main -----
# Load the network data for each time point and save as a new variable
load(file.path("enrichment_encore_DD4_Update.RData"))
TP_MF1 <- TP_MF
TP_BP1 <- TP_BP
TP_CC1 <- TP_CC
Gene4 <- gene_interest

load(file.path("enrichment_encore_DD8_Update.RData"))
TP_MF2 <- TP_MF
TP_BP2 <- TP_BP
TP_CC2 <- TP_CC
Gene8 <- gene_interest

load(file.path("enrichment_encore_DD12_Update.RData"))
TP_MF3 <- TP_MF
TP_BP3 <- TP_BP
TP_CC3 <- TP_CC
Gene12 <- gene_interest

load(file.path("enrichment_encore_DD16_Update.RData"))
TP_MF4 <- TP_MF
TP_BP4 <- TP_BP
TP_CC4 <- TP_CC
Gene16 <- gene_interest

load(file.path("enrichment_encore_DD20_Update.RData"))
TP_MF5 <- TP_MF
TP_BP5 <- TP_BP
TP_CC5 <- TP_CC
Gene20 <- gene_interest

load(file.path("enrichment_encore_DD24_Update.RData"))
TP_MF6 <- TP_MF
TP_BP6 <- TP_BP
TP_CC6 <- TP_CC
Gene24 <- gene_interest

load(file.path("enrichment_encore_DDLL_Update_2.RData"))
TP_MF7 <- TP_MF
TP_BP7 <- TP_BP
TP_CC7 <- TP_CC
GeneLL <- gene_interest

# Compile the GO terms into one big list for all time points
TP_M <- unlist(unique(c("MF", #TP_MFTP, 
                 TP_MF1, TP_MF2, TP_MF3, TP_MF4, TP_MF5, TP_MF6, TP_MF7)))
TP_B <- unlist(unique(c("BP", #TP_BPTP, 
                 TP_BP1, TP_BP2, TP_BP3, TP_BP4, TP_BP5, TP_BP6, TP_BP7)))
TP_C <- unlist(unique(c("CC", #TP_CCTP, 
                 TP_CC1, TP_CC2, TP_CC3, TP_CC4, TP_CC5, TP_CC6, TP_CC7)))


# Compile the networks into one big dataframe
allgenes <- unique(rbind(Gene4, Gene8, Gene12, Gene16, Gene20, Gene24, GeneLL))

# Map each protein across all timepoints to their associated GO term
# (returning to the Term 1 -> protein a, b, c form)
BP <- map2term(allgenes, TP_B)
MF <- map2term(allgenes, TP_M)
CC <- map2term(allgenes, TP_C)

# Re-Map the terms to their associated proteins for each time point
# however only terms with at least 20 proteins across all time points
# are included now.

# NOTE: this changes the network to include only terms 
# with a chosen number genes or more, please use the network
# from the previous script (example: enrichment_encore_DD24.Rdata)
# for a complete network with all of the GO terms (second level only)

GeneTP <- remap2gene(GeneTP, BP)
GeneTP <- remap2gene(GeneTP, MF)
GeneTP <- remap2gene(GeneTP, CC)

Gene4 <- remap2gene(Gene4, BP)
Gene4 <- remap2gene(Gene4, MF)
Gene4 <- remap2gene(Gene4, CC)

Gene8 <- remap2gene(Gene8, BP)
Gene8 <- remap2gene(Gene8, MF)
Gene8 <- remap2gene(Gene8, CC)

Gene12 <- remap2gene(Gene12, BP)
Gene12 <- remap2gene(Gene12, MF)
Gene12 <- remap2gene(Gene12, CC)

Gene16 <- remap2gene(Gene16, BP)
Gene16 <- remap2gene(Gene16, MF)
Gene16 <- remap2gene(Gene16, CC)

Gene20 <- remap2gene(Gene20, BP)
Gene20 <- remap2gene(Gene20, MF)
Gene20 <- remap2gene(Gene20, CC)

Gene24 <- remap2gene(Gene24, BP)
Gene24 <- remap2gene(Gene24, MF)
Gene24 <- remap2gene(Gene24, CC)

GeneLL <- remap2gene(GeneLL, BP)
GeneLL <- remap2gene(GeneLL, MF)
GeneLL <- remap2gene(GeneLL, CC)

# List of all the GO terms across all time points
TP_B <- unique(as.list(BP$Parent))
TP_M <- unique(as.list(MF$Parent))
TP_C <- unique(as.list(CC$Parent))


# Save data -----
# save the legend data
save(file = file.path("ALLTP_Old_Legend.RData"), list = c("TP_B", "TP_M", "TP_C"))

# save each time point to be used for the network graphing
# rename each network to the same thing to make the network graphing easier
gene_interest <- Gene4
save(file = file.path("DD4.RData"), list = c("gene_interest"))
gene_interest <- Gene8
save(file = file.path("DD8.RData"), list = c("gene_interest"))
gene_interest <- Gene12
save(file = file.path("DD12.RData"), list = c("gene_interest"))
gene_interest <- Gene16
save(file = file.path("DD16.RData"), list = c("gene_interest"))
gene_interest <- Gene20
save(file = file.path("DD20.RData"), list = c("gene_interest"))
gene_interest <- Gene24
save(file = file.path("DD24.RData"), list = c("gene_interest"))
gene_interest <- GeneLL
save(file = file.path("DDLL.RData"), list = c("gene_interest"))
gene_interest <- GeneTP
save(file = file.path("ALLTP.RData"), list = c("gene_interest"))

# Proceed to the "Network_Graphing.R" file to create the final figures
