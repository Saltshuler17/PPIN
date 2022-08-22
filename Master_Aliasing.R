## Master Aliasing Creator

library(dplyr)
library(stringr)
library(data.table)
library(tibble)
library(naniar)


aliases <- fread("5141.protein.aliases.v11.0.txt")
colnames(aliases) <- c("StringID", "ID", "ID_type")

# Limit the aliasing to Uniprot IDs for conversion to gene names
aliases <- aliases %>%
  dplyr::group_by(StringID) %>%
  filter(ID_type == "Ensembl_UniProt" | ID_type == "Ensembl_EntrezGene")

aliases_uni <- aliases %>% filter(ID_type == "Ensembl_UniProt") %>%
  dplyr::rename(Uniprot = ID, Unnecessary = ID_type)

aliases_entrez <- aliases %>% filter(ID_type == "Ensembl_EntrezGene")
# Data frame with the NCU for every entrez ID that has one
NC_entrez2ncu <- read.csv("NC_ENTREZ2NCU.csv", stringsAsFactors = F)

entrez2ncu <- data.frame("ID" = NC_entrez2ncu$Input.ID, "NCU" = NC_entrez2ncu$Gene.ID, 
                         "Name" = NC_entrez2ncu$Gene.Name.or.Symbol, stringsAsFactors = F)
entrez2ncu <- unique(entrez2ncu)

for (i in 1:length(entrez2ncu$Name)){
  if (entrez2ncu$Name[[i]] == 'N/A'){
    entrez2ncu$Name[[i]] = entrez2ncu$NCU[[i]]
  }
}

master_alias <- merge(aliases_entrez, entrez2ncu, by = "ID", all = T)
master_alias <- merge(master_alias, aliases_uni, by = "StringID", all = T)
master_alias$Unnecessary <- NULL



write.csv(master_alias, file = file.path("5141_MasterEntrezAliasing.csv"),
          row.names = FALSE)

write.table(master_alias, file = "5141_MasterEntrezAliasing.txt", sep = "\t",
            row.names = FALSE)
