##Instalar pacotes. Utilizar vers??o atualizada do R, 3.3.0
source("https://bioconductor.org/biocLite.R")
biocLite("IRanges")

source("https://bioconductor.org/biocLite.R")
biocLite("airway")

source("https://bioconductor.org/biocLite.R")
biocLite("Rsamtools")

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicAlignments")

source("https://bioconductor.org/biocLite.R")
biocLite("BiocParallel")

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

install.packages("RColorBrewer")

install.packages("pheatmap")

install.packages("ggplot2")

install.packages("PoiClaClu")

source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")

source("https://bioconductor.org/biocLite.R")
biocLite("AnnotationDbi")

source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")

source("https://bioconductor.org/biocLite.R")
biocLite("sva")

source("https://bioconductor.org/biocLite.R")
biocLite("fission")

source("https://bioconductor.org/biocLite.R")
biocLite("BiocStyle")

install.packages("knitr")

install.packages("rmarkdown")

##Importar arquivos
library("BiocStyle")
library("knitr")
library("rmarkdown")
library("gplots")
#options(width=100)

#opts_chunk$set(message = FALSE, error = FALSE, warning = FALSE, fig.width=5, fig.height=5)

library("airway")

#dir <- system.file("extdata", package="airway", mustWork=TRUE)
#dir
#list.files(dir)
#csvfile <- file.path(dir,"sample_table.csv")
#(sampleTable <- read.csv(csvfile,row.names=1))

## ----eval=FALSE----------------------------------------------------------
## sampleTable <- read.csv(csvfile,row.names=1)
## sampleTable

#(filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam")))

#file.exists(filenames)

## ------------------------------------------------------------------------
library("Rsamtools")
#bamfiles <- BamFileList(filenames, yieldSize=2000000)
#bamfiles
#seqinfo(bamfiles[1:8])

## ------------------------------------------------------------------------
library("GenomicFeatures")

#gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
#txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character() )
#txdb
#ebg <- exonsBy(txdb, by="gene")
#ebg

## ------------------------------------------------------------------------
library("GenomicAlignments")
library("BiocParallel")

## ------------------------------------------------------------------------
#register(SerialParam())

## ------------------------------------------------------------------------
#se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE)
#se


## ------------------------------------------------------------------------
data("airway")
se <- airway
se
dim(se)
assayNames(se)
?assay
head(assay(se), 3)
colSums(assay(se))
rowRanges(se)
str(metadata(rowRanges(se)))
colData(se)
colData(se) <- DataFrame(sampleTable)
se
se$cell
se$dex
se$dex <- relevel(se$dex, "untrt")
se$dex

## Inicio expessao diferencial
library("DESeq2")

#dds <- DESeqDataSet(se , design = ~ cell + dex)
countdata <- assay(se)
head(countdata, 3)
coldata <- colData(se)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ cell + dex)
nrow(ddsMat)
dds <- ddsMat[ rowSums(counts(ddsMat)) > 1, ]
nrow(dds)

## Transformacao rlog. Os valores de expressao muito baixos sao reduzidas em relacao as medias dos genes para todas as amostras
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)

## ----rldplot, fig.width=7, fig.height=3.75-------------------------------
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)

## ------------------------------------------------------------------------
##???((x1 ??? x2)?? + (y1 ??? y2)??).
sampleDists <- dist( t( assay(rld) ) )
sampleDists

## ------------------------------------------------------------------------
library("pheatmap")
library("RColorBrewer")

## ----distheatmap, fig.width=6, fig.height=4.5----------------------------
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## ----airwayDE------------------------------------------------------------
dds <- DESeq(dds)
(res <- results(dds))
mcols(res, use.names=TRUE)
summary(res)

res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

results(dds, contrast=c("cell", "N061011", "N61311"))
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
sum(res$padj < 0.1, na.rm=TRUE)
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])

## ----plotcounts----------------------------------------------------------
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("dex"))

## ----ggplotcountsjitter, fig.width=6, fig.height=4.5---------------------
data <- plotCounts(dds, gene=topGene, intgroup=c("dex","cell"), returnData=TRUE)
ggplot(data, aes(x=dex, y=count, color=cell)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=3)

## ----plotma--------------------------------------------------------------
plotMA(res, ylim=c(-5,5))

## ----plotmalabel---------------------------------------------------------
plotMA(resLFC1, ylim=c(-5,5))
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

## ----histpvalue2---------------------------------------------------------
hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")

## ------------------------------------------------------------------------
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)
topVarGenes

## ----genescluster--------------------------------------------------------
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("cell","dex")])
pheatmap(mat, annotation_col=df)

