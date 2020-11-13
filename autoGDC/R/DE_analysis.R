list.of.packages <- c("ggplot2", "DESeq2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(DESeq2)

deseq2_r <- function(counts_df, design_matrix, design_formula) {
    rownames(design_matrix) <- NULL

    # Remove small reads
    counts_df[is.na(counts_df)] <- 0
    counts_df = counts_df[rowSums(counts_df==0)<2, ]
    
    dds <- DESeqDataSetFromMatrix(countData = counts_df,
                                  colData = design_matrix,
                                  design = design_formula)
    dds <- DESeq(dds)

    resLFC <- results(dds, alpha=0.05)
    resOrdered <- resLFC[order(resLFC$padj),]
    resOrdered <- resOrdered[resOrdered$log2FoldChange > 2,]
  
    # To maintain rows when switching to pandas
    resOrdered <- as.data.frame(resOrdered)
    resOrdered["gene"] <- rownames(resOrdered)
    return(resOrdered)
}
