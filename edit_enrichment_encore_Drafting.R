# Multi-Omics Auto-GO Results
# By Hannah De los Santos
# Edited by Sam Altshuler
# Originated on: 1/20/20
# R 4.0
# load libraries and data ----

#https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script/35842176#35842176
# set working directory - only works in RStudio (with rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# If the packages aren't updated, the GO analysis will contain outdated terms
# such as Cell for CC (we already know it's in the cell, not a super useful term
# unless it is extracellular)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", version = "3.12")

# neurospora library
library(AnnotationHub)

# libraries for ontology enrichment
library(topGO)
library(stringr)
library(mygene)
library(data.table)

# ontology libraries;
library(AnnotationDbi)

# neurospora results (Loads data from other file)
#load(file.path("~/Hurley Lab/PracticeData.RData"))
add_background <- T # compare to the whole genome or not

# functions ----

prune_dag <- function(ex_graph, pvals, ont_pval_cut){
  sg <- ex_graph
  
  sg_levels <- buildLevels(sg)
  
  sg_relationships <- names(sg@edgeData@data)
  sg_relationships <- strsplit(sg_relationships,"|", fixed = T)
  
  # get children
  children <- as.character(lapply(sg_relationships, `[[`, 1))
  parents <- as.character(lapply(sg_relationships, `[[`, 2))
  
  # children parent dataframe
  cp.df <- data.frame(matrix(0,2,length(sg_relationships)))
  cp.df[1,] <- children
  cp.df[2,] <- parents
  rownames(cp.df) <- c("child","parent")
  # sort the dataframe by the children, most to least
  cp.df <- cp.df[,order(as.numeric(mget(as.character(cp.df[1,]), envir=sg_levels[["nodes2level"]])), decreasing = T)]
  
  # parallel array to the cp.df, keeps track of whether children are significant
  child_sig <- rep(F, ncol(cp.df))
  names(child_sig) <- cp.df[1,]
  keep <- rep(F, sg_levels$noOfNodes)
  names(keep) <- sg@nodes
  
  # get significance
  # pvals <- fc_results[["Damped"]][["BP"]][["go_test"]]@score # all pvalues
  pvals <- pvals[sg@nodes] # relevant pvalues
  keep <- pvals < ont_pval_cut # set only the significant pvalues
  child_sig <- pvals[as.character(cp.df[1,])] < .05
  
  # go through
  for (i in 1:ncol(cp.df)){
    if (child_sig[i]){
      keep[cp.df["parent",i]] <- T
      child_sig[cp.df["parent",i]] <- T
    }
  }
  
  prune_sg <- subGraph(names(keep)[keep],sg)
  
  return(prune_sg)
}

# ont_tax: string
# pval_type: string
# pval_cut: number
# note: has side affects, inherits from getting ontology file section
get_fc_results <- function(ont_tax, pval_type, pval_cut, ont_pval_type, ont_pval_cut, gene_focus, ont_group, map_sub.df, LOGICAL_VECT_OF_SIGNIFICANT_GENES){
  paste("Setting up ontology packages. Started on:",Sys.time())
  master_alias <- fread("Copy_5141_MasterEntrezAliasing.txt")
  gene_ids <- map_sub.df$entrezgene
  
  # get mapping for neurospora
  #(ont_tax == "5141")
  # neurospora crassa
  # neurospora library:
  ah <- AnnotationHub()
  nc_name <- which(AnnotationHub::query(ah,"OrgDb")$species == "Neurospora crassa_OR74A")
  
  # now we load in our species! nc will serve as our org.XX.eg.db package
  nc <- AnnotationHub::query(ah,"OrgDb")[[names(AnnotationHub::query(ah,"OrgDb"))[nc_name]]]
  
  kt <- "GID" # ncu numbers
  ont.df <- AnnotationDbi::select(nc,
                   keys = keys(nc,keytype = kt),
                   # columns = c("GID","GO_ID","GO_TERM_NAME"),
                   columns = c("GID", "ALIAS","GO","ONTOLOGY"),
                   keytype = kt)
  
  #Use aliasing doc to convert from Entrez to NCU
  for (i in 1:length(ont.df$GID)){
    j <- match(ont.df$GID[[i]], master_alias$ID)
    if (is.na(j)){
      next
    }
    ont.df$ALIAS[[i]] <- master_alias$NCU[[j]]
  }
  ont.df <- unique(ont.df)
  # ont.df <- ont.df %>%
  #  dplyr::filter(grepl('NCU', ont.df$ALIAS))
  # convert to a gene2go list
  gene2GO <- rep(list(c()),length(unique(ont.df$ALIAS)))
  names(gene2GO) <- unique(ont.df$ALIAS)
  for (i in 1:nrow(ont.df)){
    gene2GO[[ont.df$ALIAS[i]]] <- c(gene2GO[[ont.df$ALIAS[i]]], ont.df$GO[i]) 
  }
  
  
  # keep none significant AC category for each ont category
  no_sig <- list("BP"=c(),"CC"=c(),"MF"=c())
  
  # mod_types <- c("Linear","Exponential","ECHO","ECHO Linear","All.Osc","All.Non.Osc")
  # not 
  # fc.cats <- c("Damped", "Forced", "Harmonic", "All.Circ.wo.OE.RE")
  
  
  # no_sig[!c("BP","CC","MF") %in% ont_group] <- mod_types
  # now do analysis for each AC coefficient category
  # including all circadian without overexpressed/repressed, and all circadian
  fc_results <- list()
  gene_list <- factor(as.integer(gene_ids %in% 
                                   gene_ids[LOGICAL_VECT_OF_SIGNIFICANT_GENES])) # to pass to topgo
  # ex: log_vect <- c("A"=T,"B"=F,T,T)
  
  names(gene_list) <- gene_ids
  levels(gene_list) <- c(0,1)
  gene_list <<- gene_list
  # choices: Biological Process, Cellular Component, Molecular Function
  # now need to do it for each choice
  ontologies <- ont_group
  ontology_list <- list()
  m <- "dat"
  for (ont in ontologies){
    print(paste(m,ont))
    
    paste("Calculating enrichments for",m,ont,"ontology. Started on:",Sys.time())
    # make a topgo data object to find annotations
    go_data <- new("topGOdata",
                   description = paste(m,ont,"GO Data"),
                   ontology = ont, 
                   allGenes = gene_list,
                   nodeSize = 10,
                   annotationFun = annFUN.gene2GO,
                   gene2GO = gene2GO)
    
    if (length(go_data@graph@nodes) > 0){
      # run Fisher test for enrichment
      go_test <- runTest(go_data, algorithm = "classic", statistic = "fisher")
      
      go_levels <- buildLevels(go_data@graph) # build the dag graph so we can get levels
      
      # adjust for multiple hypothesis testing
      # before adjusting, PANTHER doesn't take into account terms that don't have at least 2 terms in sig
      go_results <- GenTable(go_data, classicFisher = go_test,
                             # orderBy = "classicFisher", 
                             # ranksOf = "classicFisher", 
                             topNodes = sum(score(go_test) <= 1))
      rownames(go_results) <- go_results$GO.ID
      go_results <- go_results[names(go_test@score),]
      # only keep the ones with at least 2 in significant
      has_2 <- go_results$Significant >= 2
      
      
      if (ont_pval_type == "BH_Adj_P-Value"){
        go_test@score[has_2] <- p.adjust(go_test@score[has_2], method = "BH")
      } # do not adjust if only pvalue is specified
      
      # prune graph based on pvalue
      prune_sg <- prune_dag(go_data@graph, go_test@score, ont_pval_cut)
      
      if (numNodes(prune_sg) > 0){
        # build new levels
        prune_sg_levels <- buildLevels(prune_sg)
        prune_levels <- as.numeric(as.list(prune_sg_levels$nodes2level))
        names(prune_levels) <- names(as.list(prune_sg_levels$nodes2level))
        # sort prune levels by name
        prune_levels <- prune_levels[order(names(prune_levels))]
        
        # generate a table of results for fisher test
        go_results <- GenTable(go_data, classicFisher = go_test,
                               # orderBy = "classicFisher", 
                               # ranksOf = "classicFisher", 
                               topNodes = sum(score(go_test) <= 1))
        # remove less thans
        go_results$classicFisher <- gsub("< ", "", go_results$classicFisher, fixed=TRUE)
        go_results$classicFisher <- as.numeric(go_results$classicFisher)
        
        # prune the graph
        go_data@graph <- prune_sg
        
        # remove results unrelated to the pruned nodes
        go_results <- go_results[go_results$GO.ID %in% prune_sg@nodes,]
        
        # sort by names
        go_results <- go_results[order(go_results$GO.ID),]
        
        # add level and sort by it
        go_results[,"Level"] <- 0
        go_results$Level <- prune_levels
        
        # sort table by level
        go_results <- go_results[order(go_results$Level),]
        
        # get fold enrichment
        go_results[["Fold.Enrichment"]] <- go_results$Significant/go_results$Expected
        
        # figure out whether everyone has a child
        # preallocate keep
        hasChild <- rep(F,nrow(go_results))
        names(hasChild) <- go_results$GO.ID
        
        # get graph
        sg <- go_data@graph
        
        # get parents and children -- parallel arrays
        sg_relationships <- names(sg@edgeData@data)
        sg_relationships <- strsplit(sg_relationships,"|", fixed = T)
        
        children <- as.character(lapply(sg_relationships, `[[`, 1))
        parents <- as.character(lapply(sg_relationships, `[[`, 2))
        
        hasChild <- names(hasChild) %in% parents
        names(hasChild) <- go_results$GO.ID
        
        go_results$hasChild <- hasChild
        
        # aggregate results
        go_list <- list(title = ont, go_data = go_data, 
                        go_test = go_test, go_results = go_results)
        # put in an overall list for each ontology
        ontology_list[[ont]] <- go_list
      } else {
        no_sig[[ont]] <- c(no_sig[[ont]],m)
      }
    } else {
      no_sig[[ont]] <- c(no_sig[[ont]],m)
    }
    
    fc_results[[m]] <- ontology_list
  }
  
  return(list(fc_results,no_sig))
}

# get data ----

# get subsets of certain range

is_all_range <- FALSE
low_range <- 20
high_range <- 24

# specifics based on organism ----
# Read in the background
all_genes <- read.csv(paste0("5141","_background.csv"),header = T,stringsAsFactors = F)
# Read in the csv from the Interactions_networking script
gene_interest <- read.csv("DDLL_Network.csv", header = T, stringsAsFactors = F)

# Only Nodes that we have a recorded NCU for and remove duplicates
# (the GO term enrichment needs NCU numbers)
final_df <- data.frame("Gene_Name" = gene_interest$NodeNCU, "Logical" = rep(T, nrow(gene_interest)), stringsAsFactors = F) %>%
  na.omit() %>%
  unique()


# add background genes not in list to total_results
if (add_background){
  id_type <- "GID"
  existing_bg <- all_genes[,id_type][all_genes[,id_type]!=""]
  
  # this will be the same for both
  background <- setdiff(existing_bg,final_df$Gene_Name)
  if (length(background) > 0){
    orig_len <- nrow(final_df)
    emp <- data.frame(matrix(NA,
                             length(background),
                             ncol(final_df)))
    
    colnames(emp) <- colnames(final_df)
    final_df <- rbind(final_df, emp)
    final_df$Gene_Name[-(1:orig_len)] <- background
    
    final_df$Logical[-(1:orig_len)] <- F
    
    # there are no duplicates
  }
}

# Node GO Analysis ----
# SOME MAPPING

paste("Mapping gene names. Started on:",Sys.time())

# neurospora has no mapping
map_sub.df <- data.frame(matrix(0,nrow(final_df),2))

colnames(map_sub.df) <- c("query","entrezgene")
# entrez gene is a misnomer for neurospora
map_sub.df$query <- map_sub.df$entrezgene  <- final_df$Gene_Name

map_sub.df <- cbind(map_sub.df)

gene_focus <- map_sub.df$query

# NO STRINGDB

# push map_sub.df to global to save
# get gene ontology enrichments for each fc category

pval_type <- "BH_Adj_P-Value"
pval_cut <- .05
ont_pval_type <- "BH_Adj_P-Value" # OR P-Value # adjustment type
ont_pval_cut <- .05 # pvalue cutoff
ont_group <- c("BP","MF","CC")

print(paste("Starting RNA", Sys.time()))
# The following line may need to be run once if the error:
# "error: package or namespace load failed for 'dplyr'"
# detach("package:dplyr", unload = T)
ont_list_n <- get_fc_results("5141",
                           pval_type,
                           pval_cut,
                           ont_pval_type,
                           ont_pval_cut,
                           gene_focus,
                           ont_group,
                           map_sub.df,
                           final_df$Logical)

fc_results_node <- ont_list_n[[1]]
no_sig_node <- ont_list_n[[2]]
#-------------------------------------------------------------------------------------
# Assigning GO Terms to the NCU's if possible (otherwise GO = Unknown)

library(topGO)
library(igraph)
# library(dplyr)
m <- "dat" 
# The different GO Categories
go_ont <- c("BP", "MF", "CC")
fam_pc <- data.frame("ChildID" = NA, "ParentID" = NA,
                     "ParentTerm" = NA, "GID" = NA, "Ont" = NA)
# Creating a vector to store second parent terms
twoparent <- vector(length = 0)
t <- 1 # Counter of the second parent terms
# Loop through each ontology category
for (ont in go_ont){
  dat <- fc_results_node[[m]][[ont]]$go_results$GO.ID[-(1)]
  # Temporary data frame for each ontology category
  temp.df <- data.frame("ChildID"= dat)
  # counter of the significant GO terms per category:
  i <- 1
  # Loop through each significantly enriched GO term to find the shortest path
  # to the second parent
  for (child in dat){
    ont_map <- c("BP"="GO:0008150",
                 "CC"="GO:0005575",
                 "MF"="GO:0003674")
    sg <- inducedGraph(fc_results_node[[m]][[ont]]$go_data@graph, child)
    # get the ID
    scndParent <- rev(names(shortest_paths(igraph.from.graphNEL(sg), child, to=ont_map[[ont]], weights = NA)$vpath[[1]]))[2]
    
    # get the name
    scndParent_Term <- fc_results_node[[m]][[ont]]$go_results[fc_results_node[[m]][[ont]]$go_results$GO.ID == scndParent,"Term"]
    # Save the term
    twoparent[[t]] <- scndParent_Term
    # save the parents
    temp.df$ParentID[[i]] <- scndParent
    temp.df$ParentTerm[[i]] <- scndParent_Term
    i <- i+1
    t <- t+1
  }
  # Find the genes that are associated with each GO Term
  for (j in 1:length(dat)){
    go_genes <- genesInTerm(fc_results_node[[m]][[ont]]$go_data, dat[[j]])
    # save them
    temp.df$GID[[j]] <- go_genes[[1]]
    temp.df$Ont[[j]] <- ont
    # ex: term1 -> A,b,c genes
  }
  # Add the temporary data frame to the data frame that combines all three ontology category
  fam_pc <- rbind(fam_pc, temp.df)
}
# Get rid of any NAs
fam_pc <- na.omit(fam_pc)

# Remove the duplicate second parent names
tp <- unique(twoparent)
MF_n <- fam_pc %>%
  dplyr::filter(Ont == "MF")
BP_n <- fam_pc %>%
  dplyr::filter(Ont == "BP")
CC_n <- fam_pc %>%
  dplyr::filter(Ont == "CC")

# Convert the data from a GO term with a list of associated genes to a gene with a list of 
# associated GO terms (from term -> a, b, c genes to A -> term1, term2, term3)
for (k in 1:length(gene_interest$NodeNCU)){
  i <- 1
  c <- 1
  
  #Assign Biological Processes terms to Genes
  gene_interest$GoID_Node_BP[[k]] <- vector("list",1)
  # Vector to store the ID's for this specific node
  GoID_Node <- vector(length = 0)
  # Loop through all the child IDs (except the first one which is the highest level GO Term)
  for (i in 1:length(BP_n$GID)){
    # If the Node is associated with that GO term, the term is added to the GOID_Node vector
    if (gene_interest$NodeNCU[[k]] %in% BP_n$GID[[i]]){
      GoID_Node[[c]] <- as.character(BP_n$ParentTerm[[i]])
      c <- c + 1
    }
  }
  # Remove duplicates
  GoID_Node <- as.list(unique(GoID_Node))
  # If the vector is not empty (has GO terms), the list of terms are assigned to that Node
  if (length(GoID_Node) > 0) {
    gene_interest$GoID_Node_BP[[k]] <- GoID_Node
  } else {
    # If no terms were found, the node is said to have unknown GO ID
    gene_interest$GoID_Node_BP[[k]] <- c("Unknown")
  }
  
  # Repeat for Molecular Function terms
  c <- 1
  gene_interest$GoID_Node_MF[[k]] <- vector("list",1)
  GoID_Node <- vector(length = 0)
  for (i in 1:length(MF_n$ChildID)){
    #print(i)
    if (gene_interest$NodeNCU[[k]] %in% MF_n$GID[[i]]){
      GoID_Node[[c]] <- as.character(MF_n$ParentTerm[[i]])
      c <- c + 1
    }
  }
  GoID_Node <- as.list(unique(GoID_Node))
  if (length(GoID_Node) > 0) {
    gene_interest$GoID_Node_MF[[k]] <- GoID_Node
  } else {
    gene_interest$GoID_Node_MF[[k]] <- c("Unknown")
  }
  
  #Repeat for Cellular Component terms
  c <- 1
  gene_interest$GoID_Node_CC[[k]] <- vector("list",1)
  GoID_Node <- vector(length = 0)
  for (i in 1:length(CC_n$ChildID)){
    #print(i)
    if (gene_interest$NodeNCU[[k]] %in% CC_n$GID[[i]]){
      GoID_Node[[c]] <- as.character(CC_n$ParentTerm[[i]])
      c <- c + 1
    }
  }
  GoID_Node <- as.list(unique(GoID_Node))
  if (length(GoID_Node) > 0) {
    gene_interest$GoID_Node_CC[[k]] <- GoID_Node
  } else {
    gene_interest$GoID_Node_CC[[k]] <- c("Unknown")
  }
}

# Source GO Analysis-----------------------------------------------------------------------------------------
#Same as above, but for the Sources instead of the nodes:
final_df_source <- data.frame("Gene_Name" = gene_interest$SourceNCU, "Logical" = rep(T, nrow(gene_interest)), stringsAsFactors = F) %>%
  na.omit() %>%
  unique()

# add background genes not in list to total_results
if (add_background){
  id_type <- "GID"
  existing_bg <- all_genes[,id_type][all_genes[,id_type]!=""]
  
  # this will be the same for both
  background <- setdiff(existing_bg,final_df_source$Gene_Name)
  if (length(background) > 0){
    orig_len <- nrow(final_df_source)
    emp <- data.frame(matrix(NA,
                             length(background),
                             ncol(final_df_source)))
    
    colnames(emp) <- colnames(final_df_source)
    final_df_source <- rbind(final_df_source, emp)
    final_df_source$Gene_Name[-(1:orig_len)] <- background
    
    final_df_source$Logical[-(1:orig_len)] <- F
    
    # there are no duplicates
  }
}

# SOME MAPPING 

paste("Mapping gene names. Started on:",Sys.time())

# neurospora has no mapping
map_sub_source.df <- data.frame(matrix(0,nrow(final_df_source),2))

colnames(map_sub_source.df) <- c("query","entrezgene")
# entrez gene is a misnomer for neurospora
map_sub_source.df$query <- map_sub_source.df$entrezgene  <- final_df_source$Gene_Name

map_sub_source.df <- cbind(map_sub_source.df)

gene_focus_s <- map_sub_source.df$query

# NO STRINGDB

# push map_sub.df to global to save
# get gene ontology enrichments for each fc category

pval_type <- "BH_Adj_P-Value"
pval_cut <- .05
ont_pval_type <- "BH_Adj_P-Value" # OR P-Value # adjustment type
ont_pval_cut <- .05 # pvalue cutoff
ont_group <- c("BP","MF","CC")

print(paste("Starting RNA", Sys.time()))
# detach("package:dplyr", unload = T)
ont_list_s <- get_fc_results("5141",
                             pval_type,
                             pval_cut,
                             ont_pval_type,
                             ont_pval_cut,
                             gene_focus,
                             ont_group,
                             map_sub_source.df,
                             final_df_source$Logical)

fc_results_s <- ont_list_s[[1]]
no_sig_s <- ont_list_s[[2]]
#-------------------------------------------------------------------------------------
# Assigning GO Terms to the NCU's if possible (otherwise GO = Unknown)

library(topGO)
library(igraph)
library(dplyr)
m <- "dat" 
go_ont <- c("BP", "MF", "CC")
fam_pc_s <- data.frame("ChildID" = NA, "ParentID" = NA,
                     "ParentTerm" = NA, "GID" = NA, "Ont" = NA)
twoparent_s <- vector(length = 0)
t <- 1
for (ont in go_ont){
  dat <- fc_results_s[[m]][[ont]]$go_results$GO.ID[-(1)]
  temp.df <- data.frame("ChildID"= dat)
  i <- 1
  for (child in dat){
    ont_map <- c("BP"="GO:0008150",
                 "CC"="GO:0005575",
                 "MF"="GO:0003674")
    sg <- inducedGraph(fc_results_s[[m]][[ont]]$go_data@graph, child)
    
    # get the ID
    scndParent <- rev(names(shortest_paths(igraph.from.graphNEL(sg), child, to=ont_map[[ont]], weights = NA)$vpath[[1]]))[2]
    # get the name
    scndParent_Term <- fc_results_s[[m]][[ont]]$go_results[fc_results_s[[m]][[ont]]$go_results$GO.ID == scndParent,"Term"]
    twoparent_s[[t]] <- scndParent_Term
    # save the parents
    temp.df$ParentID[[i]] <- scndParent
    temp.df$ParentTerm[[i]] <- scndParent_Term
    i <- i+1
    t <- t+1
  }
  j <- 1
  for (j in 1:length(dat)){
    go_genes <- genesInTerm(fc_results_s[[m]][[ont]]$go_data, dat[[j]])
    # save them
    temp.df$GID[[j]] <- go_genes[[1]]
    temp.df$Ont[[j]] <- ont
    # ex: term1 -> A,b,c genes
  }
  fam_pc_s <- rbind(fam_pc_s, temp.df)
}
Stp <- unique(twoparent_s)
MF_s <- fam_pc_s %>%
  dplyr::filter(Ont == "MF")
BP_s <- fam_pc_s %>%
  dplyr::filter(Ont == "BP")
CC_s <- fam_pc_s %>%
  dplyr::filter(Ont == "CC")

# map the genes to all the go terms they need
# example: A - term1,term2,term3
for (k in 1:length(gene_interest$SourceNCU)){
  i <- 1
  c <- 1
  gene_interest$GoID_Source_BP[[k]] <- vector("list",1)
  GoID_Source <- vector(length = 0)
  for (i in 1:length(BP_s$ChildID)){
    if (gene_interest$SourceNCU[[k]] %in% BP_s$GID[[i]]){
      GoID_Source[[c]] <- as.character(BP_s$ParentTerm[[i]])
      c <- c + 1
      
    }
  }
  GoID_Source <- as.list(unique(GoID_Source))
  if (length(GoID_Source) > 0) {
    gene_interest$GoID_Source_BP[[k]] <- GoID_Source
  } else {
    gene_interest$GoID_Source_BP[[k]] <- c("Unknown")
  }
  
  c <- 1
  gene_interest$GoID_Source_MF[[k]] <- vector("list",1)
  GoID_Source <- vector(length = 0)
  for (i in 1:length(MF_s$ChildID)){
    if (gene_interest$SourceNCU[[k]] %in% MF_s$GID[[i]]){
      GoID_Source[[c]] <- as.character(MF_s$ParentTerm[[i]])
      c <- c + 1
      
    }
  }
  GoID_Source <- as.list(unique(GoID_Source))
  if (length(GoID_Source) > 0) {
    gene_interest$GoID_Source_MF[[k]] <- GoID_Source
  } else {
    gene_interest$GoID_Source_MF[[k]] <- c("Unknown")
  }
  
  c <- 1
  gene_interest$GoID_Source_CC[[k]] <- vector("list",1)
  GoID_Source <- vector(length = 0)
  for (i in 1:length(CC_s$ChildID)){
    if (gene_interest$SourceNCU[[k]] %in% CC_s$GID[[i]]){
      GoID_Source[[c]] <- as.character(CC_s$ParentTerm[[i]])
      c <- c + 1
      
    }
  }
  GoID_Source <- as.list(unique(GoID_Source))
  if (length(GoID_Source) > 0) {
    gene_interest$GoID_Source_CC[[k]] <- GoID_Source
  } else {
    gene_interest$GoID_Source_CC[[k]] <- c("Unknown")
  }
}

#Overall list of all the GO terms in both the sources and nodes by ont category
TP_MF <- unique(append(MF_s$ParentTerm, MF_n$ParentTerm))
TP_BP <- unique(append(BP_s$ParentTerm, BP_n$ParentTerm))
TP_CC <- unique(append(CC_s$ParentTerm, CC_n$ParentTerm))
TP <- append(tp, Stp)


# Turning the list of GO terms into a string to be able to save the network as a CSV
GO_Network <- gene_interest
# Takes the list of terms for each protein and turns it into a string separated by commas
for (i in 1:length(GO_Network$Source)){
  GO_Network$GoID_Node_BP[[i]] <- as.character(paste(unlist(GO_Network$GoID_Node_BP[[i]]), collapse = ", "))
  GO_Network$GoID_Node_MF[[i]] <- as.character(paste(unlist(GO_Network$GoID_Node_MF[[i]]), collapse = ", "))
  GO_Network$GoID_Node_CC[[i]] <- as.character(paste(unlist(GO_Network$GoID_Node_CC[[i]]), collapse = ", "))
  GO_Network$GoID_Source_BP[[i]] <- as.character(paste(unlist(GO_Network$GoID_Source_BP[[i]]), collapse = ", "))
  GO_Network$GoID_Source_MF[[i]] <- as.character(paste(unlist(GO_Network$GoID_Source_MF[[i]]), collapse = ", "))
  GO_Network$GoID_Source_CC[[i]] <- as.character(paste(unlist(GO_Network$GoID_Source_CC[[i]]), collapse = ", "))
}

 # save results ----

# Save GO enrichment results
# Saved as .Rdata file due to use of lists in the GO_ID columns
save(file = file.path("enrichment_encore_DDLL_TEST_TEST.RData"), list = c("final_df","fc_results_s","fc_results_node",
                                                              "gene_interest", "TP", "TP_MF", "TP_BP", "TP_CC"))
fwrite(GO_Network, file = file.path("DDLL_Network_GO_TEST.csv"))
# Proceed to the "Legend_Creator.R" script to create the legends and 
# trim down the GO terms to only relevant terms (10 or more genes)
# NOTE: Utilize the network data table created here for all analysis beyond figure making!

