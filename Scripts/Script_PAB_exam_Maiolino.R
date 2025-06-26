# Load libraries
library(Maiolino_PAB_exam)

# Step 1: Load sparse matrix and convert to full matrix
dt <- s1.make.dt(data_path = "/sharedFolder/Data/matrix/")


# Step 2: Split gene expression and ATAC-seq data
dt_genes <- s2.subset.dt(string = "ENSG", data = dt)
dt_atac <- s2.subset.dt(string = "chr", data = dt)


# Step 3: Summarize Data___
genes_summary <- s3.summarize.data(data = dt_genes)
atac_summary <- s3.summarize.data(data = dt_atac)


# Step 4: Create GenomicRanges
gr_genes <- s4.create.GenomicRanges("genes")
gr_atac <- s4.create.GenomicRanges("atac")


# Step 5 ___________________________________________________________________________
h38_coding <- s5.GR.protein.coding()

gr_atac <- s5.remap.atac()


# Step 6 ___________________________________________________________________________
gr_genes <- s6.finalize.exp.data()


# Step 7 ___________________________________________________________________________
# Step 7 - Normalize _______________________________________________________________
genes_cpm <- s7.normalize.cmp(data = genes_summary)
atac_cpm <- s7.normalize.cmp(data = atac_summary)

gr_genes$gene_cpm <- genes_cpm[match(gr_genes$gene_id, names(genes_cpm))]
gr_atac$peak_cpm  <- atac_cpm[match(gr_atac$peak_id, names(atac_cpm))]

# Step 7 - Merge ___________________________________________________________________
merged_data <- s7.merge()

# Step 7 - Summary _________________________________________________________________
unmerged_peaks <- s7.summary.table.peaks()

genes_no_atac <- s7.summary.table.genes()

s7.summary.plot.peaks <- function(
plot_s7_pea <- s7.summary.plot.peaks()

plot_s7_genes <- s7.summary.plot.genes(merged_data)

chr_wise_summary <- s7.summary.chr.wise()


# Step 8 ___________________________________________________________________________
s8.save.plots()