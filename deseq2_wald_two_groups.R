library (openxlsx)
library (DESeq2)
library (ggplot2)


## See Github RNA-Seq_mouse/gene_annotation.R
system ("cp /projects/ncrrbt_share_la/dev_pipe/gencode.vM32.annotation.txt .")

anno <- read.delim ("gencode.vM32.annotation.txt")

anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)


## normal STAR results from the RNA-Seq IIT pipeline
a <- read.xlsx ("star_gene_raw_counts.xlsx")
a <- a[ ,grep ("Gene|bam", colnames (a))]

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 

a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("Aligned.out.bam", "", colnames (a))
colnames (a) <- gsub ("star.", "", colnames (a))


counts <- annot <- a

annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]

row.names (counts) <- counts$Geneid
counts <- counts[ ,grep ("S", colnames (counts))]

samples <- data.frame (matrix (nrow=dim (counts)[2], ncol=3))
colnames (samples) <- c("sample", "condition", "sex")
samples$sample <- colnames (counts)
samples$sex <- as.factor (substr (gsub ("\\..*" , "", colnames (counts)), 1, 1))
samples$condition <- as.factor (substr (gsub ("\\..*" , "", colnames (counts)), 2, 2))

samples <- samples[!samples$sample %in% c("mF.01_S10", "fF.02_S11"), ]
samples

counts <- counts[ ,colnames (counts) %in% samples$sample]

idx <- match (samples$sample, colnames (counts))
samples <- samples[idx, ]
stopifnot (samples$sample == colnames (counts))



## TEspeX counts from the outfile.txt

tesp <- read.delim ("outfile.txt", row.names=1)
colnames (tesp) <- gsub ("_R.*", "", colnames (tesp))
tesp <- tesp[ ,colnames (tesp) %in% colnames (counts)]
stopifnot (colnames (tesp) == colnames (counts))

counts <- rbind (counts, tesp)



## DESeq2 

dds <- DESeqDataSetFromMatrix(countData = round (counts), colData = samples, design = ~ sex +condition)
                                 
# keep <- rowSums(counts(dds)) >= 10
keep <- rowSums(counts(dds) >= 10) >= dim (counts)[2]/2
dds <- dds[keep,]
dds

# first contrast of interest (F vs C)
ddsLRT <- DESeq(dds, test="LRT", full=~sex+condition, reduced=~sex)
resultsNames(ddsLRT)

res <- results(ddsLRT, contrast=list("condition_F_vs_C"), test="Wald")

res <- merge (data.frame (res), counts (dds), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid", all.x=TRUE)

res$gene_name [is.na (res$gene_name)] <- res$Row.names [is.na (res$gene_name)]
res$external_gene_name [is.na (res$external_gene_name)] <- res$Row.names [is.na (res$external_gene_name)]
res$gene_type[is.na (res$gene_type)] <- paste ("transposon", gsub (".*#", "", res$Row.names [is.na (res$gene_type)]), sep=":")

colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

write.xlsx (res, "hippocampus_deseq2_FentanylvsControl_with_transposons_differential_expression.xlsx", rowNames=F)


# second contrast of interest (W vs C)
ddsLRT <- DESeq(dds, test="LRT", full=~sex+condition, reduced=~sex)
resultsNames(ddsLRT)

res <- results(ddsLRT, contrast=list("condition_W_vs_C"), test="Wald")

res <- merge (data.frame (res), counts (dds), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid", all.x=TRUE)

res$gene_name [is.na (res$gene_name)] <- res$Row.names [is.na (res$gene_name)]
res$external_gene_name [is.na (res$external_gene_name)] <- res$Row.names [is.na (res$external_gene_name)]
res$gene_type[is.na (res$gene_type)] <- paste ("transposon", gsub (".*#", "", res$Row.names [is.na (res$gene_type)]), sep=":")

colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

write.xlsx (res, "hippocampus_deseq2_WithdrawalvsControl_with_transposons_differential_expression.xlsx", rowNames=F)



# third contrast of the two previous interaction terms, using the list() style of contrasts !!!
ddsLRT <- DESeq(dds, test="LRT", full=~sex+condition, reduced=~sex)
resultsNames(ddsLRT)

res <- results(ddsLRT, contrast=list("condition_W_vs_C", "condition_F_vs_C"), test="Wald")  ## This is equivalent to W/F


res <- merge (data.frame (res), counts (dds), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

write.xlsx (res, "hippocampus_deseq2_WithdrawalvsFentanyl_with_transposons_differential_expression.xlsx", rowNames=F)












