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
    features <- as.data.table(
        read_tsv(
            file = path_f, 
            col_names = c("id", "name", "type", "chr", "start", "end"), 
            show_col_types = F
        )
    )
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

    gr <- GRanges(
        seqnames = features$chr,
        ranges = IRanges(
            start = features$start,
            end = features$end
        ),
        sum = summary
    )

    if (type == "Genes") {
        gr$gene_id <- features_genes$id
        gr$gene_name <- features_genes$name
    } else if (type == "Peaks") {
        gr$peak_id <- features_atac$id
    }

    return(gr)
}


gr_genes <- s4.create.GenomicRanges("genes")
gr_atac <- s4.create.GenomicRanges("atac")


# Step 5 ___________________________________________________________________________
s5.GR.protein.coding <- function(
    path = "/sharedFolder/Data/Homo_sapiens.GRCh38.114.gtf.gz"
) {
    data <- rtracklayer::import("Data/Homo_sapiens.GRCh38.114.gtf.gz")

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

    data_coding <- data[data$type == "gene" & data$gene_biotype == "protein_coding"]

    return(data_coding)
}

h38_coding <- s5.GR.protein.coding()

