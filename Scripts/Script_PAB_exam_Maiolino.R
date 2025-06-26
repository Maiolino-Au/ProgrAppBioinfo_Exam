# Load libraries
library(PABMaiolinoPackage)

# Step 1: Load sparse matrix and convert to full matrix
dt <- s1.make.dt(data_path = "/sharedFolder/Data/matrix/")

# Step 2: Split gene expression and ATAC-seq data
dt_genes <- s2.subset.dt(string = "ENSG", data = dt)
dt_atac <- s2.subset.dt(string = "chr", data = dt)

# Step 3: Summarize Data
genes_summary <- s3.summarize.data(data = dt_genes)
atac_summary <- s3.summarize.data(data = dt_atac)

# Step 4: Create GenomicRanges
gr_genes <- s4.create.GenomicRanges("genes")
gr_atac <- s4.create.GenomicRanges("atac")

# Step 5: Gene Annotation for ATACseq data
# Step 5 - Create GFT
h38_coding <- s5.GR.protein.coding()
# Step 5: Remap ATAC
gr_atac <- s5.remap.atac()

# Step 6 Finalize expression data
gr_genes <- s6.finalize.exp.data()

# Step 7 Data Normalization and Integration
# Step 7 - Normalize 
genes_cpm <- s7.normalize.cmp(data = genes_summary)
atac_cpm <- s7.normalize.cmp(data = atac_summary)

gr_genes$gene_cpm <- genes_cpm[match(gr_genes$gene_id, names(genes_cpm))]
gr_atac$peak_cpm  <- atac_cpm[match(gr_atac$peak_id, names(atac_cpm))]

# Step 7 - Merge 
merged_data <- s7.merge()

# Step 7 - Summary 
unmerged_peaks <- s7.summary.table.peaks()

genes_no_atac <- s7.summary.table.genes()

plot_s7_peak <- s7.summary.plot.peaks()

plot_s7_genes <- s7.summary.plot.genes()

chr_wise_summary <- s7.summary.chr.wise()

# Step 8: Visualization
s8.save.plots()
