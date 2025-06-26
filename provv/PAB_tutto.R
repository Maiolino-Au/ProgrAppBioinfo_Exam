# Load libraries
library(Seurat)
library(Signac)
library(Matrix)
library(readr)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(rtracklayer)

# Step 1 ___________________________________________________________________________
s1.make.dt <- function(
    data_path = "/sharedFolder/Data/matrix/"
) {
    # Load the sparse matrix
    matrix <- readMM(file = paste0(data_path, "matrix.mtx.gz"))
    features <- read_tsv(file = paste0(data_path, "features.tsv.gz"), 
                         col_names = c("id", "name", "type", "chr", "start", "end"), 
                         show_col_types = F)
    barcodes <- read_tsv(file = paste0(data_path, "barcodes.tsv.gz"), 
                         col_names = F, 
                         show_col_types = F) 

    # Assign names to colums and rows and transform the object "matrix" in a matrix
    colnames(matrix) <- barcodes$X1
    rownames(matrix) <- features$id
    matrix <- as.matrix(matrix)

    # Genearate a data.table, a id columns containg the ids of the various genes/peaks is added 
    dt <- as.data.table(matrix, row.names = rownames(matrix))
    dt[, id := rownames(matrix)]
    setcolorder(dt, c("id", setdiff(names(dt), "id"))) # id is put as the first column


    return(dt)
}

dt <- s1.make.dt(data_path = "/sharedFolder/Data/matrix/")


# Step 2 ___________________________________________________________________________
s2.subset.dt <- function(
    string,
    data
) {
    # Use grepl to create a logic that indicates which elemnent of the data.table to take
    sub <- data[grepl(string, data$id)]
    return(sub)
}

dt_genes <- s2.subset.dt(string = "ENSG", data = dt)
dt_atac <- s2.subset.dt(string = "chr", data = dt)


# Step 3 ___________________________________________________________________________
s3.summarize.data <- function(
    data
) {
    summary <- data[, rowSums(.SD), .SDcols = -"id"]
    names(summary) <- data$id
    return(summary)
}

genes_summary <- s3.summarize.data(data = dt_genes)
atac_summary <- s3.summarize.data(data = dt_atac)


# Step 4 ___________________________________________________________________________
s4.create.GenomicRanges <- function(
    type,
    path_f = "/sharedFolder/Data/matrix/features.tsv.gz"
) {
    # Load Features
    features <- as.data.table(
        read_tsv(
            file = path_f, 
            col_names = c("id", "name", "type", "chr", "start", "end"), 
            show_col_types = F
        )
    )

    # Subset features for genes or peaks
    if (type %in% c("Genes", "genes")) {
        if (type != "Genes") {type = "Genes"}
        
        features <- features[grepl("Gene", features$type)]
        features$chr[is.na(features$chr)] <- "unknown"

        summary = genes_summary
    } else if (type %in% c("Peaks", "peaks", "ATAC", "Atac", "atac")) {
        if (type != "Peaks") {type = "Peaks"}
        
        features <- features[grepl("Peaks", features$type) & grepl("chr", features$chr)]
        features$chr[is.na(features$chr)] <- "unknown"

        summary = atac_summary
    } else {
        print("Error, invalid value for type: must be 'Gene' or 'Peaks'")
        return ("")
    }

    # Create GRanges object
    gr <- GRanges(
        seqnames = features$chr,
        ranges = IRanges(
            start = features$start,
            end = features$end
        )
    )

    # Add metadata 
    if (type == "Genes") {
        gr$gene_sum = summary
        gr$gene_id <- features$id
        gr$gene_name <- features$name
    } else if (type == "Peaks") {
        gr$peak_sum = summary
        gr$peak_id <- features$id
    }

    return(gr)
}

gr_genes <- s4.create.GenomicRanges("genes")
gr_atac <- s4.create.GenomicRanges("atac")


# Step 5 ___________________________________________________________________________
s5.GR.protein.coding <- function(
    path = "/sharedFolder/Data/Homo_sapiens.GRCh38.114.gtf.gz"
) {
    # Import the GTF file
    data <- rtracklayer::import("Data/Homo_sapiens.GRCh38.114.gtf.gz")

    # Standardize the names (in the files the chromosomes are called by number/letter)
    seqlevels(data) <- ifelse(
        grepl(
            "^([1-9]|1[0-9]|2[0-2]|X|Y)$", # [1-22] would not work but just check from 1 to 2, with the second 2 ignored
            seqlevels(data)
        ), 
        paste0(
            "chr", 
            seqlevels(data)
        ), 
        seqlevels(data)
    )

    # Filter for protein coding genes
    data_coding <- data[data$type == "gene" & data$gene_biotype == "protein_coding"]

    return(data_coding)
}

h38_coding <- s5.GR.protein.coding()

s5.remap.atac <- function(
    atac_data = gr_atac
) {
    # Find the atac peaks that overlap with coding genes
    overlaps <- findOverlaps(atac_data, h38_coding)
    
    # Create a list of the genes for which there is an overlap
    gene_list <- split(h38_coding$gene_id[subjectHits(overlaps)], queryHits(overlaps))
    
    # Assign the genes to the corresponding peaks
    # In the case a peak is overlapping multiple genes they are added as a comma separated list
    atac_data$overlapping_genes <- NA
    atac_data$overlapping_genes[as.integer(names(gene_list))] <- sapply(gene_list, paste, collapse = ",") 

    return(atac_data)
}

gr_atac <- s5.remap.atac()


# Step 6 ___________________________________________________________________________
s6.finalize.exp.data <- function(
    exp_data = gr_genes
) {
    # Add a column to with gene symbols by matching the gene_id to the annotation from the GTF file
    exp_data$gene_symbol <- h38_coding$gene_name[match(exp_data$gene_id, h38_coding$gene_id)]

    # # Remove rows where gene_symbol is NA
    exp_data <- exp_data[!is.na(exp_data$gene_symbol)]

    return(exp_data)
}

gr_genes <- s6.finalize.exp.data()


# Step 7 ___________________________________________________________________________
# Step 7 - Normalize _______________________________________________________________
s7.normalize.cmp <- function(data) {
    log2((data / sum(data) * 1e6) + 1)
}

genes_cpm <- s7.normalize.cmp(data = genes_summary)
atac_cpm <- s7.normalize.cmp(data = atac_summary)

gr_genes$gene_cpm <- genes_cpm[match(gr_genes$gene_id, names(genes_cpm))]
gr_atac$peak_cpm  <- atac_cpm[match(gr_atac$peak_id, names(atac_cpm))]

# Step 7 - Merge ___________________________________________________________________
s7.merge <- function(
    exp_data = gr_genes,
    atac_data = gr_atac
) {
    merged <- merge(
      as.data.table(exp_data)[, .(gene_id, gene_symbol, chr = seqnames, gene_cpm)],
      as.data.table(atac_data)[, .(overlapping_genes, peak_cpm)],
      by.x = "gene_id", by.y = "overlapping_genes", all.x = TRUE
    )

    return(merged)
}

merged_data <- s7.merge()

# Step 7 - Summary _________________________________________________________________
s7.summary.table.peaks <- function(
    atac_data = gr_atac
) {
    unmerged <- as.data.table(atac_data)[is.na(overlapping_genes)]
    unmerged<- unmerged[, .(peak_id, seqnames, start, end, peak_cpm)]
    colnames(unmerged)[colnames(unmerged) == "seqnames"] <- "chr"
        
    return(unmerged)
}

unmerged_peaks <- s7.summary.table.peaks()

s7.summary.table.genes <- function(
    exp_data = gr_genes,
    merged = merged_data
) {
    genes_no_atac <- merged[is.na(peak_cpm) | peak_cpm == 0]
    genes_no_atac <- genes_no_atac[, .(gene_id, gene_symbol, gene_cpm)]
    
    return(genes_no_atac)
}

genes_no_atac <- s7.summary.table.genes()

s7.summary.plot.peaks <- function(
    data_merged = merged_data,
    data_unmerged = unmerged_peaks
) {
    # Provide "a plot of peak intensity distribution chromosome by chromosome"
    # For peaks that could be merged or for the ones that could not?
    

    merged <- ggplot(data_merged[!is.na(peak_cpm)], aes(x=chr, y=peak_cpm)) + 
        geom_boxplot() + theme(axis.text.x = element_text(angle=90)) +
        labs(title = "Peak intensity distribution for merged peaks")

    unmerged <- ggplot(data_unmerged[!is.na(peak_cpm)], aes(x=chr, y=peak_cpm)) + 
        geom_boxplot() + theme(axis.text.x = element_text(angle=90)) +
        labs(title = "Peak intensity distribution for unmerged peaks")
    
    #Save
    dir_results <- "/sharedFolder/Results"
    if (!dir.exists(dir_results)) {dir.create(dir_results)}
    
    ggsave("/sharedFolder/Results/s7_peaks_intesnisty_merged.png", plot = merged)
    ggsave("/sharedFolder/Results/s7_peaks_intesnisty_unmerged.png", plot = unmerged)
    
    return(list(merged, unmerged))
}

plot_s7_pea <- s7.summary.plot.peaks()

s7.summary.plot.genes <- function(
    data = merged_data
) {
    plot <- list()
    
    # atac
    data_atac <- data[peak_cpm > 0 & grepl("chr", data$chr)]
    merged <- ggplot(data_atac, aes(x=chr, y=gene_cpm)) + 
        geom_boxplot() + theme(axis.text.x = element_text(angle=90)) +
        labs(title = "Expression of Genes with ATAC Peaks")
    
    # no atac
    data_no_atac <- data[(is.na(peak_cpm)| peak_cpm == 0) & grepl("chr", data$chr)]
    unmerged <- ggplot(data_no_atac, aes(x=chr, y=gene_cpm)) + 
        geom_boxplot() + theme(axis.text.x = element_text(angle=90)) +
        labs(title = "Expression of Genes without ATAC Peaks")

    #Save
    dir_results <- "/sharedFolder/Results"
    if (!dir.exists(dir_results)) {dir.create(dir_results)}
    
    ggsave("/sharedFolder/Results/s7_gene_expression_merged.png", plot = merged)
    ggsave("/sharedFolder/Results/s7_gene_expression_unmerged.png", plot = unmerged)
    
    return(list(merged, unmerged))
}

plot_s7_genes <- s7.summary.plot.genes(merged_data)

s7.summary.chr.wise <- function(

) {
    peak_chr_summary <- unmerged_peaks[, .N, by = chr]
    colnames(peak_chr_summary) <- c("chr", "N_unmerged_peaks")
    
    gr_genes_dt <- as.data.table(gr_genes)
    genes_no_atac <- merge(genes_no_atac, gr_genes_dt[, .(gene_id, seqnames)], by = "gene_id", all.x = TRUE)
    gene_chr_summary <- genes_no_atac[, .N, by = seqnames] 
    colnames(gene_chr_summary) <- c("chr", "N_unmerged_genes")
    
    chr_wise_summary <- merge(gene_chr_summary,peak_chr_summary, all = T)

    return(chr_wise_summary)
}

chr_wise_summary <- s7.summary.chr.wise()


# Step 8 ___________________________________________________________________________
s8.scatter.plots <- function(chr_num) {
    plot <- ggplot(
        merged_data[!is.na(peak_cpm) & merged_data$chr == paste0("chr", chr_num)], 
        aes(x=gene_cpm, y=peak_cpm)
    ) +
      geom_point(alpha=0.5) +
      labs(
          x = "Gene Expression CPM (log2)", 
          y = "ATAC Peaks CPM (log2)", 
          title = paste("Chromosme", chr_num)
      ) +
      #facet_wrap(~ chr) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return(plot)
}

s8.save.plots <- function() {
    plot_list <- lapply(c(as.character(1:22), "X", "Y"), s8.scatter.plots)

    dir_results <- "/sharedFolder/Results"
    if (!dir.exists(dir_results)) {dir.create(dir_results)}

    for (i in 1:(length(plot_list)/12)) {
        n <- 12*(i-1)+1
        combined <- cowplot::plot_grid(plotlist = plot_list[(n):(n+11)], ncol = 3) 
        ggsave(paste0(
            "/sharedFolder/Results/s8_plot_",
            i,
            ".png"
        ), plot = combined, width = 1920*2, height = 1080*3.6, units = "px")
    }
}

s8.save.plots()