# Protein-Protein Interaction Networking
# By Sam Altshuler
# Originated on: 05/09/20
# R 3.6 to use STRINGDB v 10 for correct mapping
# Rtools35 
# load libraries and data ----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", version = "3.10")
##Install the string database
if (!requireNamespace("STRINGdb", quietly = TRUE))
  BiocManager::install("STRINGdb")

if (!("rvest" %in% installed.packages())) {
  install.packages("rvest")
}

# Libraries that need to be installed through BioConductor:
# - STRINGdb
# - topGO
# - mygene
# - AnnotationDbi
# Install these using this style:
# BiocManager::install("Package to be installed")

#Libraries for web scraping
library(rvest)

#Libraries for mapping
library(STRINGdb)
library(dplyr)
library(stringr)
library(data.table)


#Libraries for GO Mapping
library(topGO)
library(mygene)
library(igraph)

library(AnnotationDbi)

#https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script/35842176#35842176
# set working directory - only works in RStudio (with rstudioapi)
#setting your working directory to whatever folder the document is coming from
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in the data from a csv, the data should be in the form of two columns 
# The first column consisting of Uniprot Accession ID for Neurospora (example: Q1K502_NEUCR)
# The second column consisting of NCU#.
# Make sure that the number of NCU# matches the number of Neurospora Uniprot IDs
gene_n <- read.csv("DD4.csv", stringsAsFactors = F, na.strings = c("", NA))
colnames(gene_n) <- c("GID", "st_ID")
# Removing all "NA" NCU (no NCU found for the non-Neurospora proteins)
NCU <- na.omit(gene_n$st_ID)

gene_names <- gene_n
# Remove the non-Neurospora proteins from the list, 
# not needed if the list is prescreened 
# gene_names <- gene_n %>%
#   dplyr::filter(grepl('NEUCR', gene_n$GID))%>%
#   unique()

# Replace the NCU column with the list of NCU's with NA's omitted 
# (NA values correspond to non-Neurospora proteins). If there is a
# different number of NCU from Uniprot IDs at this point, an error
# will pop up indicating there are more of one or the other. Manually
# check that the two columns have the same number of rows.
gene_names$st_ID <- NCU

stringdb <- STRINGdb$new(version="10",
                         species=5141,
                         score_threshold=0,
                         input_directory="")

#-----------------------------------------------------------------------------------------
#Small example with just FRQ
frq <- "NCU02265"
mapdf <- data.frame("frq"="FRQ") # c("a","b")
# Creates the map with the string database (string database)
mapdf <- stringdb$map(mapdf, "frq", removeUnmappedRows = F)

# Create a data table using the neurospora crassa file from string
# Read in the NC file from string, this file contains the links between
# proteins (all links contained in the string database for NC)
links <- fread("5141.protein.links.full.v11.0.txt")

# Read in the protein aliases from String (linking the string IDs to 
# other forms of identification)
master_alias <- fread("Copy_5141_MasterEntrezAliasing.txt")


# Limit the aliasing to Uniprot IDs for conversion to gene names


# Create a map that has all proteins that interact with the protein in mapdf
# and shows how they  interact
temp <- links[links$protein1 %in% mapdf$STRING_id | links$protein2 %in% mapdf$STRING_id,]

# Can filter out interaction types
# temp <- temp[temp$homology > 0,]

#END OF EXAMPLE
#----------------------------------------------------------------------------------------

## use NCU to get what I can, then use Uniprot IDs to fill in the rest
map_all_3 <- stringdb$map(gene_names, "st_ID", removeUnmappedRows = F)

# Looks to see what proteins don't have a STRING ID from NCUs and uses
# the Uniprot accession ID to fill in the blanks
for (i in 1:length(map_all_3$GID)){
  if (is.na(map_all_3$STRING_id[[i]])){
    if (length(stringdb$mp(gene_names$GID[[i]])) > 0) {
      map_all_3$STRING_id[[i]] <- stringdb$mp(gene_names$GID[[i]])
    } else {
      map_all_3$STRING_id[[i]] <- NA
    }
  }
}
# If the above loop does not fill in all of the NAs, use the master
# aliasing document to fill in the remaining gaps (last resort since
# you want to use the STRINGdb directly rather than the doc which might
# be outdated)
for (i in 1:length(map_all_3$GID)){
  if (is.na(map_all_3$STRING_id[[i]])){
    j <- match(gene_names$st_ID[[i]], master_alias$NCU)
    if (is.na(j)){
      map_all_3$STRING_id[[i]] <- map_all_3$GID[[i]]
    } else{
      map_all_3$STRING_id[[i]] <- master_alias$StringID[[j]]
    }
  }
}

# If there are any gaps in the data (for instance a term having no 
# STRING id), that term is removed from the analysis. However, this 
# should never happen unless there was a mistake in the input data
map_all_3[complete.cases(map_all_3), ]

# read in the NC file from string
links <- fread("5141.protein.links.full.v11.0.txt")

# Create a map that has all proteins that interact with the proteins in map_all_3
# and shows how they  interact
temp1 <- links[links$protein1 %in% map_all_3$STRING_id | links$protein2 %in% map_all_3$STRING_id,]

# Filter out interaction types, make sure that you only have links that are 
# experimentally determined and score 990/1000 or higher in the combined score
temp1 <- temp1[temp1$experiments_transferred > 0,]
temp1 <- temp1[temp1$combined_score > 990,]


# Create a list of sources and targets (targets are called nodes in this script):
# Every source for primary interactions is FRQ (or the protein of interest)
source <- vector(length = 1)
node <- vector(length = 1)
connections <- vector(length = 1)
degree <- vector(length = 1)

# Assign the source and node values, if there is no STRING_id, use the uniprot ID
for (i in 1:(length(map_all_3$GID))){
  # Every source fo the primary interactions is FRQ
  source[[i]] <- map_all_3$STRING_id[[1]]
  degree[[i]] <- "primary"
  connections[[i]] <- "a"
  if (is.na(map_all_3$STRING_id[[i]])){
    node[[i]] <- map_all_3$GID[[i]]
  } else{
    node[[i]] <- map_all_3$STRING_id[[i]]
  }
}

# Add in the secondary interactors from string
# Since the interactions aren't directional, it doesn't matter
# if the protein is defined as source or node
source <- append(source, temp1$protein1)
node <- append(node, temp1$protein2)

# Label these interactions as secondary
for (i in (length(map_all_3$GID)+1):length(node)){
  connections[[i]] <- "a"
  degree[[i]] <- "secondary"
}

# Combine everything into a data frame with source and target nodes
interactions_pre <- data.table("Source" = source,
                           "Connection" = connections, 
                           "Node" = node, 
                           "Degree" = degree)

# Get rid of recipricol relationships 
# Example: x --> y and y --> x are the same because the interactions
# are not direction and this would double count the interaction
interactions_pre <- interactions_pre %>%
  mutate("Key" = paste0(pmin(source, node), pmax(source,node), sep= " "))

interactions <- interactions_pre[!duplicated(interactions_pre$Key),]

# Nullify the Key column as it is not needed
interactions$Key <- NULL

#------------------------------------------------------------------------------
# Naming the Source and Nodes from their string ID to protein names
# Get the protein names from the Aliases file

Source_name <- vector(length = 1)
Source_ncu <- vector(length = 1)

# Counter for how many string IDs do not have uniprot IDs in the file
p <- 0

# Loop through all source IDs 
for (i in 1:length(interactions$Source)){
  k <- 0
  # Loop through the aliases file until it finds a matching string ID
  for (j in 1:length(master_alias$StringID)){
    # If a matching string ID is found, the Uniprot ID is added to 
    # a vector storing the names in the order of the string IDs
    if (master_alias$StringID[[j]] == interactions$Source[[i]]) {
      Source_name[[i]]<- master_alias$Name[[j]]
      Source_ncu[[i]] <- master_alias$NCU[[j]]
      k <- 1
      # Don't continue this inner loop once a match is found
      break
    }
  }
  # If not match is found, the name is not available
  # These IDs would have to manually be added
  if (k == 0){
    p <- p+1
    Source_name[[i]] <- NA
  }
}
cat("Total Source StringID's unable to be converted to Uniprot: ", p)
# Attach the string IDs to the Uniprot IDs
Source_df <- data.frame("StringID" = interactions$Source, "NCU" = Source_ncu,
                        "Name" = Source_name)

# Get rid of repeats and omit the NA values (No NCU or no name)
Source_df_n <- unique(Source_df)
Source_df_uniq <- na.omit(Source_df_n)

# Repeat the process from above, but with the target nodes instead
Node_name <- vector(length = 1)
Node_ncu <- vector(length = 1)
l <- 0
for (x in 1:length(interactions$Node)){
  k <- 0
  for (y in 1:length(master_alias$StringID)){
    if (master_alias$StringID[[y]] == interactions$Node[[x]]) {
      Node_name[[x]]<- master_alias$Name[[y]]
      Node_ncu[[x]]<- master_alias$NCU[[y]]
      k <- 1
      break
    }
  }
  if (k == 0){
    l <- l+1
    Node_name[[x]] <- NA
  }
}
Node_df <- data.frame("StringID" = interactions$Node, "NCU"= Node_ncu,
                      "Name" = Node_name)
cat("Total Node StringID's unable to be converted to Uniprot: ", l)

Node_df_n <- unique(Node_df)
Node_df_uniq <- na.omit(Node_df_n)


#-------------------------------------------------------------------------------
# Assign the gene names and NCUs to the interactions


# Assign the gene name and NCU to each source in the network
for (q in 1:length(interactions$Source)){
  t <- 1
  # Loop through the sources data frame to find the match (STRING ID to STRING ID)
  for (w in 1:length(Source_df_uniq$StringID)){
    # Once a match is found, give it the gene name and NCU
    if (interactions$Source[[q]] == Source_df_uniq$StringID[[w]]){
      interactions$SourceName[[q]] <- as.character(Source_df_uniq$Name[[w]])
      interactions$SourceNCU[[q]] <- as.character(Source_df_uniq$NCU[[w]])
      t <- 0
      break
    }
  }
  # If there is no match, the name and NCU are set to NA 
  # This also can be manually edited in the final CSV
  if (t == 1){
    interactions$SourceName[[q]] <- NA
    interactions$SourceNCU[[q]] <- NA
  }
}

#Same as above, but for the nodes
for (q in 1:length(interactions$Node)){
  t <- 1
  w <- 1
  for (w in 1:length(Node_df_uniq$StringID)){
    if (interactions$Node[[q]] == Node_df_uniq$StringID[[w]]){
      interactions$NodeNCU[[q]] <- as.character(Node_df_uniq$NCU[[w]])
      interactions$NodeName[[q]] <- as.character(Node_df_uniq$Name[[w]])
      t <- 0
      break
    }
  }
  if (t == 1){
    interactions$NodeName[[q]] <- NA
    interactions$NodeNCU[[q]] <- NA
  }
}

#Get rid of any duplicate interactions 
interactions <- unique(interactions)

# adds in known values 
for(i in 1:length(interactions$Source)){
  if (interactions$Source[[i]] %in% map_all_3$STRING_id){
    for (j in 1:length(map_all_3$STRING_id)){
      if(is.na(map_all_3$STRING_id[[j]])){
        map_all_3$STRING_id[[j]] <- map_all_3$GID[[j]]
      }
      if (interactions$Source[[i]] == map_all_3$STRING_id[[j]]){
        interactions$SourceNCU[[i]] <- map_all_3$st_ID[[j]]
        break
      }
    }
  }
}
for(i in 1:length(interactions$Node)){
  if (interactions$Node[[i]] %in% map_all_3$STRING_id){
    for (j in 1:length(map_all_3$STRING_id)){
      if(is.na(map_all_3$STRING_id[[j]])){
        map_all_3$STRING_id[[j]] <- map_all_3$GID[[j]]
      }
      if (interactions$Node[[i]] == map_all_3$STRING_id[[j]]){
        interactions$NodeNCU[[i]] <- map_all_3$st_ID[[j]]
        if (is.na(interactions$NodeName[[i]])){
          interactions$NodeName[[i]] <- map_all_3$st_ID[[j]]
        }
      }
    }
  }
}

#---------------------------------------------------------------------------------------
# save results, this CSV will show up in the file designated here
write.csv(interactions, file = file.path("DD4_Network_Update.csv"), 
          row.names = FALSE)
print("Protein Protein Interaction Network saved")

# Proceed to the "edit_enrichment_encore.R" script to add in GO term analysis
