suppressPackageStartupMessages(library("DESeq2"))               # Load the DESeq2 Library.
suppressPackageStartupMessages(library("BiocParallel"))         # Enable Parallelization.
register(MulticoreParam(48))                                    # Register the amount of CPU cores to use.

params = commandArgs(trailingOnly = TRUE)
countPath = params[1]
conditions = c("co_culture", "co_culture", "co_culture", "monoculture", "monoculture", "monoculture", "monoculture", "co_culture")
targetPath = params[2]


string_splitter <- function(string, pattern, indices) {
    elements <- strsplit(string, pattern)[[1]]
    result_elements <- elements[indices]
    result <- paste(result_elements, collapse = pattern)
    return(result)
}


countFile = read.csv(countPath, skip = 1, row.names = 1, sep="\t")

names <- c()
for (filename in colnames(countFile)[6:length(colnames(countFile))])
{
  names <- cbind(names, string_splitter(string_splitter(filename, "reads.", 2), ".sam"))
}

countMatrix <- data.frame(countFile[,6:length(colnames(countFile))])
colnames(countMatrix) = names

sampleInformation = matrix(colnames(countMatrix))                                                          
sampleInformation = data.frame(cbind(sampleInformation, conditions, seq(1,length(sampleInformation))))   
colnames(sampleInformation) = c("sample", "condition", "sampleID")                                

dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = sampleInformation, design = ~ condition)   # Initialize the deseq2 object using the countMatrix object and the sampleInformation object.
dds$condition <- factor(dds$condition, levels = c("monoculture","co_culture"))                                                 # Correctly set the levels of the contrasted conditions of the deseq2 object.

dds <- DESeq(dds, quiet = TRUE, parallel=TRUE, BPPARAM=MulticoreParam(48))                                  # Start calculations.
res <- results(dds, parallel=TRUE, BPPARAM=MulticoreParam(48))                                              # Retrive results.
resOrdered <- res[order(res$padj),]                                                                         # Order results, sorted ascending by the adjusted p-Value column of the result object.

write.csv(as.data.frame(resOrdered), file= paste(targetPath, "deseq2_result.csv", sep=""))         # Write the result to the disk as a csv file.
write.csv(as.data.frame(dds$sizeFactor), file= paste(targetPath, "sample_size_factors.csv", sep=""))         # Write the sample size normalization factors to the disk as a csv file.

###2D PCA Plots
suppressPackageStartupMessages(library("ggplot2"))
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
symbol_numbers = c(15:25,0:14,33:127)

generate_PCA_plots <- function(ntop = 10)
{
    data <- plotPCA(vsd,intgroup=c("condition","sample"), ntop = ntop, returnData=TRUE)
    percentVar <- round(100 * attr(data, "percentVar"))
    ggplot(data,aes(PC1, PC2,color=condition,shape=sample),margin=c(13, 13)) +
        scale_shape_manual(values=symbol_numbers[1:8]) +
        geom_point(size=4) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed() +
        theme(legend.box = "horizontal", plot.title = element_text(size = 20, face = "bold", hjust=0.5)) +
        ggtitle(paste("Principal Compontant Analysis", "\nbased on the", ntop, "most variant features."))
    ggsave(paste(targetPath,"PCA_top_",ntop,".pdf",sep=""), width = 15, height = 15)
    ggsave(paste(targetPath,"PCA_top_",ntop,".png",sep=""), width = 15, height = 15)
} 

generate_PCA_plots(10)
generate_PCA_plots(50)
generate_PCA_plots(100)
generate_PCA_plots(500)

###MA Plot
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[[2]], type="normal")

pdf(paste(targetPath,"MA_plot.pdf",sep=""), width = 8, height = 8)
plotMA(resLFC, ylim = c(-3, 3))
dev.off()

###Volcano Plot
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
pdf(paste(targetPath,"volcano_plot.pdf",sep=""), width = 8, height = 8)
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
dev.off()
