## assign cell to clonotype

library(Seurat)
library(tidyverse)
library(broom)

# Update lib name and processed seurat object every time
lib_name = "D13_D33_D45_D47_CD4_Proliferated_L2"
seurat_obj <- readRDS("robjects/cellbender_CCregressed_D13_D33_D45_D47_CD4_Proliferated_L2_nF_500_4000_nC_150_post_demultiplexed_full_processed.rds")

# data from sequencing facility
tcr <- read.csv(paste0("../Data/", lib_name, "/filtered_contig_annotations.csv"))
clono <- read.csv(paste0("../Data/", lib_name, "/clonotypes.csv"))

# keep TCR and clonotype info for each barcode
tcr <- tcr %>% 
  distinct(barcode, .keep_all = T) %>% 
  dplyr::select(barcode, raw_clonotype_id) %>% 
  dplyr::rename(clonotype_id = raw_clonotype_id) %>% 
  merge(., clono) %>%
  tibble::column_to_rownames("barcode")

# keep only the barcodes that have tcr info
seurat_obj <- seurat_obj[,rownames(tcr)]

# add tcr info to object meta data
clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)

# extract meta data
meta <- clono_seurat@meta.data %>% 
  tibble::rownames_to_column("Cell") %>% 
  dplyr::select(-HTO_maxID, -HTO_secondID, -HTO_margin, -HTO_classification, 
                -doublet_logLikRatio, -donor_id, -best_singlet, -best_doublet,
                -prob_doublet, -prob_max)

write.table(meta, paste0("Expt2_", lib_name, "_metadata.txt"), row.names = F, quote = F)