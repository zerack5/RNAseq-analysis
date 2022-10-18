library(Seurat)
sce1 <- CreateSeuratObject(Read10X('C:/Users/a/Desktop/as分选课题/公共单细胞/singel cell data/AS'), "as")
sce2 <- CreateSeuratObject(Read10X('C:/Users/a/Desktop/as分选课题/公共单细胞/singel cell data/HC'), "hc")


require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)

# peak-bc matrix
mex_dir_path <- "C:/Users/a/Desktop/as分选课题/公共单细胞/singel cell data/AS"

mtx_path <- paste(mex_dir_path, "GSM4771370_AS_PBMC_matrix.mtx", sep = '/')
feature_path <- paste(mex_dir_path, "GSM4771370_AS_PBMC_peaks.bed", sep = '/')
barcode_path <- paste(mex_dir_path, "GSM4771370_AS_PBMC_barcodes.tsv", sep = '/')

features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)

mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)

# tf-bc matrix
mex_dir_path <- "C:/Users/a/Desktop/as分选课题/公共单细胞/singel cell data/HC"

mtx_path <- paste(mex_dir_path, "GSM4771371_NC_PBMC_matrix.mtx", sep = '/')
feature_path <- paste(mex_dir_path, "GSM4771371_NC_PBMC_peaks.bed", sep = '/')
barcode_path <- paste(mex_dir_path, "GSM4771371_NC_PBMC_barcodes.tsv", sep = '/')

features <- readr::read_tsv(feature_path, col_names = c('feature', 'common_name'))
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)

mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)








