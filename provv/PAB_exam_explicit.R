# Load libraries
library(Seurat)
library(Signac)
library(Matrix)
library(readr)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(dplyr)
#library(rtracklayer) per i ::


# Step 1: Load sparse matrix and convert to full matrix
data_path <- "Data/matrix/"

matrix <- readMM(file = paste0(data_path, "matrix.mtx.gz"))
features <- read_tsv(file = paste0(data_path, "features.tsv.gz"), col_names = c("id", "name", "type", "chr", "start", "end"), show_col_types = F)
barcodes <- read_tsv(file = paste0(data_path, "barcodes.tsv.gz"), col_names = F, show_col_types = F) 

colnames(matrix) <- barcodes$X1
rownames(matrix) <- features$id
matrix <- as.matrix(matrix)

dt <- as.data.table(matrix, row.names = rownames(matrix))
dt[, id := rownames(matrix)]
setcolorder(dt, c("id", setdiff(names(dt), "id")))


# Step 2: Split gene expression and ATAC-seq data
dt_genes <- dt[grepl("ENSG", dt$id)]
dt_atac <- rbind(dt[grepl("chr", dt$id)],dt[!grepl("ENSG|chr", dt$id)])


# Step 3: Summarize Data
genes_summary <- dt_genes[, rowSums(.SD), .SDcols = -"id"]
names(genes_summary) <- dt_genes$id

atac_summary <- dt_atac[, rowSums(.SD), .SDcols = -"id"]
names(atac_summary) <- dt_atac$id


# Step 4: Create GenomicRanges
features_dt <- as.data.table(features)
features_genes <- features_dt[grepl("Gene", features_dt$type)]
features_atac <- features_dt[grepl("Peaks", features_dt$type)]

features_genes$chr[is.na(features_genes$chr)] <- "unknown"

gr_genes <- GRanges(
    seqnames = features_genes$chr,
    ranges = IRanges(
        start = features_genes$start,
        end = features_genes$end
    ),
    gene_id = dt_genes$id,
    name = features_genes$name,
    sum = genes_summary
)

gr_atac <- GRanges(
    seqnames = features_atac$chr,
    ranges = IRanges(
        start = features_atac$start,
        end = features_atac$end
    ),
    peack_id = features_atac$id,
    sum = atac_summary
)


# Step 5: Gene Annotation for ATACseq data
h38 <- rtracklayer::import("Data/Homo_sapiens.GRCh38.114.gtf.gz")
seqlevels(h38) <- ifelse(grepl("^([1-9]|1[0-9]|2[0-2]|X|Y)$", seqlevels(h38)), paste0("chr", seqlevels(h38)), seqlevels(h38)) # [1-22] would not work but just check from 1 to 2, with the second 2 ignored
h38_protein_coding <- h38[h38$type == "gene" & h38$gene_biotype == "protein_coding"]

overlap_atac <- findOverlaps(gr_atac, h38_protein_coding)

gr_atac$nearest_gene <- NA
gr_atac$nearest_gene[queryHits(overlap_atac)] <- h38_protein_coding$gene_id[subjectHits(overlap_atac)]


# Step 6: Finalize expression data
gr_genes$gene_symbol <- h38_protein_coding$gene_name[match(gr_genes$gene_id, h38_protein_coding$gene_id)]
gr_genes <- gr_genes[!is.na(gr_genes$gene_symbol)]


# Step 7: Data Normalization and Integration
dt_genes_cmp <- NormalizeData(dt_genes[,-"id"], normalization.method = "LogNormalize", scale.factor = 1e6) %>% as.data.table
dt_genes_cmp[, id := dt_genes$id]
setcolorder(dt_genes_cmp, c("id", setdiff(names(dt_genes_cmp), "id")))

dt_atac_cmp <- NormalizeData(dt_atac[,-"id"], normalization.method = "LogNormalize", scale.factor = 1e6) %>% as.data.table
dt_atac_cmp[, id := dt_atac$id]
setcolorder(dt_atac_cmp, c("id", setdiff(names(dt_atac_cmp), "id")))




# Step 8: Visualization

