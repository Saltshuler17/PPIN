
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
map2term20 <- function(allgenes, TP){
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
    if (length(glist) > 0) {
      masterdf$Genes[[i]] <- glist
      masterdf$Leng[[i]] <- length(glist)
    } else {
      masterdf$Leng[[i]] <- length(glist)
      lil[[t]] <- i
      t <- t +1
    }
    masterdf$Type[[i]] <- TP[[1]]
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
# load(file.path("enrichment_encore_ALLTP.RData"))
# GeneTP <- gene_interest
# 
# load(file.path("Union_Intersect.RData"))


load(file.path("DD4.RData"))
TP_MF1 <- TP_MF
TP_BP1 <- TP_BP
TP_CC1 <- TP_CC
Gene4 <- gene_interest

load(file.path("DD8.RData"))
TP_MF2 <- TP_MF
TP_BP2 <- TP_BP
TP_CC2 <- TP_CC
Gene8 <- gene_interest

load(file.path("DD12.RData"))
TP_MF3 <- TP_MF
TP_BP3 <- TP_BP
TP_CC3 <- TP_CC
Gene12 <- gene_interest

load(file.path("DD16.RData"))
TP_MF4 <- TP_MF
TP_BP4 <- TP_BP
TP_CC4 <- TP_CC
Gene16 <- gene_interest

load(file.path("DD20.RData"))
TP_MF5 <- TP_MF
TP_BP5 <- TP_BP
TP_CC5 <- TP_CC
Gene20 <- gene_interest

load(file.path("DD24.RData"))
TP_MF6 <- TP_MF
TP_BP6 <- TP_BP
TP_CC6 <- TP_CC
Gene24 <- gene_interest

load(file.path("DDLL.RData"))
TP_MF7 <- TP_MF
TP_BP7 <- TP_BP
TP_CC7 <- TP_CC
GeneLL <- gene_interest

TP_M <- unlist(unique(c("MF", TP_MF1, TP_MF2, TP_MF3, TP_MF4, TP_MF5, TP_MF6, TP_MF7, "Unknown")))
TP_B <- unlist(unique(c("BP", TP_BP1, TP_BP2, TP_BP3, TP_BP4, TP_BP5, TP_BP6, TP_BP7, "Unknown")))
TP_C <- unlist(unique(c("CC", TP_CC1, TP_CC2, TP_CC3, TP_CC4, TP_CC5, TP_CC6, TP_CC7, "Unknown")))


# MAKE THIS A FUNCTION INSTEAD

allgenes_DD <- unique(rbind(Gene4, Gene8, Gene12, Gene16, Gene20, Gene24))
allgenes <- unique(rbind(Gene4, Gene8, Gene12, Gene16, Gene20, Gene24, GeneLL))
BP <- map2term20(allgenes, TP_B)
MF <- map2term20(allgenes, TP_M)
CC <- map2term20(allgenes, TP_C)

# make it just the primary interactions
Gene4 <- Gene4 %>%
  filter(Degree == "primary")

Gene8 <- Gene8 %>%
  filter(Degree == "primary")

Gene12 <- Gene12 %>%
  filter(Degree == "primary")

Gene16 <- Gene16 %>%
  filter(Degree == "primary")

Gene20 <- Gene20 %>%
  filter(Degree == "primary")

Gene24 <- Gene24 %>%
  filter(Degree == "primary")

# GeneInter <- GeneInter %>%
#   filter(Degree == "primary")
# 
# GeneTP <- GeneTP %>%
#   filter(Degree == "primary")
# 
# GeneTPnLL <- GeneTPnLL %>%
#   filter(Degree == "primary")

GeneLL <- GeneLL %>%
  filter(Degree == "primary")

# GeneLLnTP <- GeneLLnTP %>%
#   filter(Degree == "primary")

# NOTE: this changes the network to include only terms 
# with a chosen number genes or more, please use the network
# from the previous script (example: DD24.Rdata)
# for a complete network with all of the GO terms (second level only)
Gene4 <- remap2gene(Gene4, BP)
Gene4 <- remap2gene(Gene4, MF)
Gene4 <- remap2gene(Gene4, CC)
BP_4 <- map2term(Gene4, TP_B)
MF_4 <- map2term(Gene4, TP_M)
CC_4 <- map2term(Gene4, TP_C)

Gene8 <- remap2gene(Gene8, BP)
Gene8 <- remap2gene(Gene8, MF)
Gene8 <- remap2gene(Gene8, CC)
BP_8 <- map2term(Gene8, TP_B)
MF_8 <- map2term(Gene8, TP_M)
CC_8 <- map2term(Gene8, TP_C)


Gene12 <- remap2gene(Gene12, BP)
Gene12 <- remap2gene(Gene12, MF)
Gene12 <- remap2gene(Gene12, CC)
BP_12 <- map2term(Gene12, TP_B)
MF_12 <- map2term(Gene12, TP_M)
CC_12 <- map2term(Gene12, TP_C)


Gene16 <- remap2gene(Gene16, BP)
Gene16 <- remap2gene(Gene16, MF)
Gene16 <- remap2gene(Gene16, CC)
BP_16 <- map2term(Gene16, TP_B)
MF_16 <- map2term(Gene16, TP_M)
CC_16 <- map2term(Gene16, TP_C)


Gene20 <- remap2gene(Gene20, BP)
Gene20 <- remap2gene(Gene20, MF)
Gene20 <- remap2gene(Gene20, CC)
BP_20 <- map2term(Gene20, TP_B)
MF_20 <- map2term(Gene20, TP_M)
CC_20 <- map2term(Gene20, TP_C)


Gene24 <- remap2gene(Gene24, BP)
Gene24 <- remap2gene(Gene24, MF)
Gene24 <- remap2gene(Gene24, CC)
BP_24 <- map2term(Gene24, TP_B)
MF_24 <- map2term(Gene24, TP_M)
CC_24 <- map2term(Gene24, TP_C)


GeneLL <- remap2gene(GeneLL, BP)
GeneLL <- remap2gene(GeneLL, MF)
GeneLL <- remap2gene(GeneLL, CC)
BP_LL <- map2term(GeneLL, TP_B)
MF_LL <- map2term(GeneLL, TP_M)
CC_LL <- map2term(GeneLL, TP_C)


# GeneTP <- remap2gene(GeneTP, BP)
# GeneTP <- remap2gene(GeneTP, MF)
# GeneTP <- remap2gene(GeneTP, CC)
# BP_TP <- map2term(GeneTP, TP_B)
# MF_TP <- map2term(GeneTP, TP_M)
# CC_TP <- map2term(GeneTP, TP_C)
# 
# 
# GeneTPnLL <- remap2gene(GeneTPnLL, BP)
# GeneTPnLL <- remap2gene(GeneTPnLL, MF)
# GeneTPnLL <- remap2gene(GeneTPnLL, CC)
# BP_TPnLL <- map2term(GeneTPnLL, TP_B)
# MF_TPnLL <- map2term(GeneTPnLL, TP_M)
# CC_TPnLL <- map2term(GeneTPnLL, TP_C)
# 
# GeneLLnTP <- remap2gene(GeneLLnTP, BP)
# GeneLLnTP <- remap2gene(GeneLLnTP, MF)
# GeneLLnTP <- remap2gene(GeneLLnTP, CC)
# BP_LLnTP <- map2term(GeneLLnTP, TP_B)
# MF_LLnTP <- map2term(GeneLLnTP, TP_M)
# CC_LLnTP <- map2term(GeneLLnTP, TP_C)
# 
# 
# GeneInter <- remap2gene(GeneInter, BP)
# GeneInter <- remap2gene(GeneInter, MF)
# GeneInter <- remap2gene(GeneInter, CC)
# BP_Inter <- map2term(GeneInter, TP_B)
# MF_Inter <- map2term(GeneInter, TP_M)
# CC_Inter <- map2term(GeneInter, TP_C)


TP_B <- unique(as.list(BP$Parent))
TP_M <- unique(as.list(MF$Parent))
TP_C <- unique(as.list(CC$Parent))


# Save data -----
# save the legend data
# save(file = file.path("ALLTPLL_Legend.RData"), list = c("TP_B", "TP_M", "TP_C"))

# save each time point to be used for the network graphing
gene_interest <- Gene4
terms_go <- unique(rbind(BP_4, MF_4, CC_4))

# save(file = file.path("DD4_20.RData"), list = c("gene_interest", "terms_go"))
for (i in 1:length(terms_go$Parent)){
  terms_go$Genes[[i]] <- as.character(paste(unlist(terms_go$Genes[[i]]), collapse = ", "))
}
fwrite(terms_go, file = file.path("DD4_GOterms_primary.csv"), 
          row.names = FALSE)

gene_interest <- Gene8
terms_go <- unique(rbind(BP_8, MF_8, CC_8))
# save(file = file.path("DD8_20.RData"), list = c("gene_interest", "terms_go"))
for (i in 1:length(terms_go$Parent)){
  terms_go$Genes[[i]] <- as.character(paste(unlist(terms_go$Genes[[i]]), collapse = ", "))
}
fwrite(terms_go, file = file.path("DD8_GOterms_primary.csv"), 
          row.names = FALSE)

gene_interest <- Gene12
terms_go <- unique(rbind(BP_12, MF_12, CC_12))
# save(file = file.path("DD12_20.RData"), list = c("gene_interest", "terms_go"))
for (i in 1:length(terms_go$Parent)){
  terms_go$Genes[[i]] <- as.character(paste(unlist(terms_go$Genes[[i]]), collapse = ", "))
}
fwrite(terms_go, file = file.path("DD12_GOterms_primary.csv"), 
          row.names = FALSE)

gene_interest <- Gene16
terms_go <- unique(rbind(BP_16, MF_16, CC_16))
# save(file = file.path("DD16_20.RData"), list = c("gene_interest", "terms_go"))
for (i in 1:length(terms_go$Parent)){
  terms_go$Genes[[i]] <- as.character(paste(unlist(terms_go$Genes[[i]]), collapse = ", "))
}
fwrite(terms_go, file = file.path("DD16_GOterms_primary.csv"), 
          row.names = FALSE)

gene_interest <- Gene20
terms_go <- unique(rbind(BP_20, MF_20, CC_20))
# save(file = file.path("DD20_20.RData"), list = c("gene_interest", "terms_go"))
for (i in 1:length(terms_go$Parent)){
  terms_go$Genes[[i]] <- as.character(paste(unlist(terms_go$Genes[[i]]), collapse = ", "))
}
fwrite(terms_go, file = file.path("DD20_GOterms_primary.csv"), 
          row.names = FALSE)

gene_interest <- Gene24
terms_go <- unique(rbind(BP_24, MF_24, CC_24))
# save(file = file.path("DD24_20.RData"), list = c("gene_interest", "terms_go"))
for (i in 1:length(terms_go$Parent)){
  terms_go$Genes[[i]] <- as.character(paste(unlist(terms_go$Genes[[i]]), collapse = ", "))
}
fwrite(terms_go, file = file.path("DD24_GOterms_primary.csv"), 
          row.names = FALSE)

gene_interest <- GeneLL
terms_go <- unique(rbind(BP_LL, MF_LL, CC_LL))
# save(file = file.path("DDLL_20.RData"), list = c("gene_interest", "terms_go"))
for (i in 1:length(terms_go$Parent)){
  terms_go$Genes[[i]] <- as.character(paste(unlist(terms_go$Genes[[i]]), collapse = ", "))
}
fwrite(terms_go, file = file.path("LL_GOterms_primary.csv"), 
          row.names = FALSE)

# gene_interest <- GeneTP
# terms_go <- unique(rbind(BP_TP, MF_TP, CC_TP))
# # save(file = file.path("ALLTP_20.RData"), list = c("gene_interest", "terms_go"))
# for (i in 1:length(terms_go$Parent)){
#   terms_go$Genes[[i]] <- as.character(paste(unlist(terms_go$Genes[[i]]), collapse = ", "))
# }
# fwrite(terms_go, file = file.path("ALLTP_GOterms_primary.csv"), 
#           row.names = FALSE)


# Saving the intersections and differences between AllTP and LL

# terms_go <- unique(rbind(BP_TPnLL, MF_TPnLL, CC_TPnLL))
# for (i in 1:length(terms_go$Parent)){
#   terms_go$Genes[[i]] <- as.character(paste(unlist(terms_go$Genes[[i]]), collapse = ", "))
# }
# fwrite(terms_go, file = file.path("INTPnotLL.csv"),
#        row.names = FALSE)
# 
# terms_go <- unique(rbind(BP_LLnTP, MF_LLnTP, CC_LLnTP))
# for (i in 1:length(terms_go$Parent)){
#   terms_go$Genes[[i]] <- as.character(paste(unlist(terms_go$Genes[[i]]), collapse = ", "))
# }
# fwrite(terms_go, file = file.path("INLLnotTP.csv"),
#        row.names = FALSE)
# 
# terms_go <- unique(rbind(BP_Inter, MF_Inter, CC_Inter))
# for (i in 1:length(terms_go$Parent)){
#   terms_go$Genes[[i]] <- as.character(paste(unlist(terms_go$Genes[[i]]), collapse = ", "))
# }
# fwrite(terms_go, file = file.path("INBOTH.csv"),
#        row.names = FALSE)


# Proceed to the "Network_Graphing.R" file to create the final figures
