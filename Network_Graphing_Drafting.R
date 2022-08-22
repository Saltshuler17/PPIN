# Protein-Protein Interaction Networking Graphing
# By Sam Altshuler
# Originated on: 05/09/20

# Replaces the need for Cytoscape
library(igraph)
library(qgraph)
library(RColorBrewer)
library(pdftools)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#Load in the data from the GO enrichment file
load(file.path("DD16.RData"))
load(file.path("ALLTP_Old_Legend.RData"))
# Interact_map is a copy of the data frame with all of the information about the interactions
interact_map <- gene_interest

# If you want just the primary interactions in your figure, uncomment the filter function
# interact_map <- interact_map %>%
#   dplyr::filter(Degree == "primary")

# Data frame with only the source and target nodes
sn.df <- data.frame("Source" = interact_map$Source, "Node" = interact_map$Node, stringsAsFactors = F)
# give each source it's name
sname <- unique(data.frame("id" = sn.df$Source, "label" = interact_map$SourceName, stringsAsFactors = F))
i <- 1
for (i in 1:length(sname$id)){
  # If there is no name for the source in the interact_map data frame, it is given
  # the string ID as a name
  if (is.na(sname$id[[i]])){
    sname$label[[i]] <- sname$id[[i]]
  }
}
# Repeat above with the nodes
nname <- unique(data.frame("id" = sn.df$Node, "label" = interact_map$NodeName, stringsAsFactors = F))
for (i in 1:length(nname)){
  if (is.na(nname$label[[i]])){
    nname$label[[i]] <- nname$id[[i]]
  }
}

# Data frame with all of the vertices on the graph and their names
# (A vertex is a single point on the graph)
na_me <- unique(rbind(sname,nname))
#na_me <- na_me[-157,]
# Initialize the network graph based on the interactions from the 
# sn.df and the vertex names from na_me
mynetwork <- graph.data.frame(sn.df, directed = FALSE,
                              vertices = na_me) #, 
                      
# Add the GO terms to the source and node name data frames
sname_go <- data.frame("id" = sn.df$Source, 
                       "label" = interact_map$SourceName)
sname_go$GO_MF <- interact_map$GoID_Source_MF
sname_go$GO_BP <- interact_map$GoID_Source_BP
sname_go$GO_CC <- interact_map$GoID_Source_CC
sname_go <- unique(sname_go)
nname_go <- data.frame("id" = sn.df$Node, 
                       "label" = interact_map$NodeName)#,
                       #"GO" = interact_map$GoID_Node))
nname_go$GO_MF <- interact_map$GoID_Node_MF
nname_go$GO_BP <- interact_map$GoID_Node_BP
nname_go$GO_CC <- interact_map$GoID_Node_CC
nname_go <- unique(nname_go)
na_me_go <- unique(rbind(sname_go, nname_go))
# Add 'Unknown' to the list of GO terms and remove duplicates
TP_MF <- append(TP_M, "Unknown")
TP_MF <- unique(TP_MF)
TP_BP <- append(TP_B, "Unknown")
TP_BP <- unique(TP_BP)
TP_CC <- append(TP_C, "Unknown")
TP_CC <- unique(TP_CC)
# Interaction matrix with the GO terms as rows and the vertices as columns
val_MF <- matrix(ncol = length(na_me_go$id), nrow = length(TP_MF), dimnames = list(TP_MF, na_me_go$id))
val_BP <- matrix(ncol = length(na_me_go$id), nrow = length(TP_BP), dimnames = list(TP_BP, na_me_go$id))
val_CC <- matrix(ncol = length(na_me_go$id), nrow = length(TP_CC), dimnames = list(TP_CC, na_me_go$id))
# Loop through every vertex 
for (i in 1:length(na_me_go$id)){
  # Loop through every GO term
  for (j in 1:length(TP_MF)){
    # If the go term is associated with the vertex, that spot on the matrix 
    # is assigned a 1 (otherwise it will be a zero)
    if (TP_MF[[j]] %in% na_me_go$GO_MF[[i]]){
      val_MF[j,i] <- 1
    } else {
      val_MF[j,i] <- 0 
    }
  }
}
for (i in 1:length(na_me_go$id)){
  # Loop through every GO term
  for (j in 1:length(TP_BP)){
    # If the go term is associated with the vertex, that spot on the matrix 
    # is assigned a 1 (otherwise it will be a zero)
    if (TP_BP[[j]] %in% na_me_go$GO_BP[[i]]){
      val_BP[j,i] <- 1
    } else {
      val_BP[j,i] <- 0 
    }
  }
}
for (i in 1:length(na_me_go$id)){
  # Loop through every GO term
  for (j in 1:length(TP_CC)){
    # If the go term is associated with the vertex, that spot on the matrix 
    # is assigned a 1 (otherwise it will be a zero)
    if (TP_CC[[j]] %in% na_me_go$GO_CC[[i]]){
      val_CC[j,i] <- 1
    } else {
      val_CC[j,i] <- 0 
    }
  }
}
# Turn the interaction matrix into a data frame and then a list
values_MF <- as.list(as.data.frame(val_MF))
values_BP <- as.list(as.data.frame(val_BP))
values_CC <- as.list(as.data.frame(val_CC))

# Initialize the PDF that will be getting created, choose the title
# Error will pop up if this file is already open, make sure it is closed
# before trying to make any edits to the graph

V(mynetwork)$label <- na_me$label
pdf("DDLL_Old_labels20_BP.pdf")

# Create the layout format:
# Can edit the area and the repulsion radius to fit each specific graph
# based on the number of nodes (may require some trial and error). However
# editing the vertex size in the plot() function will offer best results
e <- get.edgelist(mynetwork, names = F)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(mynetwork), area = 200* vcount(mynetwork),
                                       repulse.rad = (vcount(mynetwork)^2.6), layout.control = "spring")
# Create the colors for each GO Term
nbcol <- length(TP_BP)-1
pal <- colorRampPalette(brewer.pal(8,"RdYlBu"))(nbcol)
pal <- append(pal, "#D3D3D3")
# Create the actual plot:
# Each vertex will be a pie chart showing what GO terms it is associated with
# Edit vertex.size and vertex.label.cex to change the vertex size and label size
# defaults for vertex size is 30 and label size is 1
plot(mynetwork,layout=l , rescale= T, vertex.shape = "pie", vertex.pie = values_BP,
     vertex.pie.color = list(pal), #edge.width = 0.15,
     vertex.size= 3.0, vertex.frame.color = NA, 
     vertex.label.font = 2, 
     vertex.label.color = "black", vertex.label.cex = 0.11)

# Add a legend to associate the colors with the GO terms
legend("topleft",
       legend = (TP_BP),
       fill = pal,
       cex = 0.28)

# Run dev.off() to create the file! Automatically saves the file
dev.off()

# Initialize the PDF that will be getting created, choose the title
# Error will pop up if this file is already open, make sure it is closed
# before trying to make any edits to the graph
pdf("DDLL_Old_labels20_MF.pdf") #, width = 100, height = 100)

# Create the layout format:
# Can edit the area and the repulsion radius to fit each specific graph
# based on the number of nodes (may require some trial and error). However
# editing the vertex size in the plot() function will offer best results
e <- get.edgelist(mynetwork, names = F)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(mynetwork), area = 200* vcount(mynetwork),
                                       repulse.rad = (vcount(mynetwork)^2.6), layout.control = "spring")
# Create the colors for each GO Term
nbcol <- length(TP_MF)-1
pal <- colorRampPalette(brewer.pal(8,"RdYlBu"))(nbcol)
pal <- append(pal, "#D3D3D3")
# Create the actual plot:
# Each vertex will be a pie chart showing what GO terms it is associated with
# Edit vertex.size and vertex.label.cex to change the vertex size and label size
# defaults for vertex size is 30 and label size is 1
plot(mynetwork,layout=l , rescale= T, vertex.shape = "pie", vertex.pie = values_MF,
     vertex.pie.color = list(pal), #edge.width = 0.15,
     vertex.size= 3.0, vertex.frame.color = NA, 
     vertex.label.font = 2,
     vertex.label.color = "black", vertex.label.cex = 0.11)

# Add a legend to associate the colors with the GO terms
legend("topleft",
       legend = (TP_MF),
       fill = pal,
       cex = 0.28)

# Run dev.off() to create the file! Automatically saves the file
dev.off()

# Initialize the PDF that will be getting created, choose the title
# Error will pop up if this file is already open, make sure it is closed
# before trying to make any edits to the graph
pdf("DDLL_Old_labels20_CC.pdf")

# Create the layout format:
# Can edit the area and the repulsion radius to fit each specific graph
# based on the number of nodes (may require some trial and error). However
# editing the vertex size in the plot() function will offer best results
e <- get.edgelist(mynetwork, names = F)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(mynetwork), area = 200* vcount(mynetwork),
                                       repulse.rad = (vcount(mynetwork)^2.6), layout.control = "spring")
# Create the colors for each GO Term
nbcol <- length(TP_CC)-1
pal <- colorRampPalette(brewer.pal(8,"RdYlBu"))(nbcol)
pal <- append(pal, "#D3D3D3")
# Create the actual plot:
# Each vertex will be a pie chart showing what GO terms it is associated with
# Edit vertex.size and vertex.label.cex to change the vertex size and label size
# defaults for vertex size is 30 and label size is 1
plot(mynetwork,layout= l, rescale= T, vertex.shape = "pie", vertex.pie = values_CC, 
     vertex.pie.color = list(pal), #edge.width = 0.15,
     vertex.size= 3.0, vertex.frame.color = NA, 
     vertex.label.font = 2,
     vertex.label.color = "black", vertex.label.cex = 0.11)

# Add a legend to associate the colors with the GO terms
legend("topleft",
       legend = (TP_CC),
       fill = pal,
       cex = 0.28)

# Run dev.off() to create the file! Automatically saves the file
dev.off()

# pdf("LLTP_Legend_MF.pdf")
# nbcol <- length(TP_MF)
# pal <- colorRampPalette(brewer.pal(8,"RdYlBu"))(nbcol+1)
# plot(NULL, xaxt = "n", yaxt = "n", byt = "n", ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
# legend("topleft",
#        legend = (TP_MF),
#        fill = pal, bty = "n",
#        cex = 0.28)
# dev.off()

#Convert from PDF to .tiff

pdf_convert("DDLL_Old_labels20_BP.pdf", format = "tiff", dpi = 1450)
pdf_convert("DDLL_Old_labels20_MF.pdf", format = "tiff", dpi = 1450)
pdf_convert("DDLL_Old_labels20_CC.pdf", format = "tiff", dpi = 1450)

# You're Done! Enjoy your figures and have a swell day