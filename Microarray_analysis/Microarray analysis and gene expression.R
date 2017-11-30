#http://www.ncbi.nlm.nih.gov/geo/browse/
#GEO2R 
#library(limma)
#library(GEOquery)

myGSE = "GSE45643" # qin, demo 

# call data

GSE45643 <- getGEO("GSE45643", GSEMatrix =TRUE)
GPL10558 <- getGEO("GPL10558", GSEMatrix =TRUE)
if (length(GSE45643) > 1) idx <- grep("GPL10558", attr(GSE45643, "Symbol")) else idx <- 1
gSET <- GSE45643[[idx]]
ex <- exprs(gSET) #This is the expression matrix

# Find out probes 
dictionary = gSET@featureData@data[, c('ID', "Symbol")]  #This is a lookup table for probe ID and gene symbol
rownames(ex)<-dictionary$Symbol

######## choosing a phenotype and compare ir
GSE45643 <- getGEO('GSE45643',GSEMatrix = TRUE)
colnames(pData(phenoData(GSE45643[[1]])))
pData(phenoData(GSE35452[[1]]))[,c(1:32)]
my_dataset<-ex[,c(1:6)]

################ welcome limma

# average of each probe
eset2 <-avereps(my_dataset,ID=rownames(my_dataset))
key<-c(rep("control",3),rep("E2",3))
key <- factor(key, levels = unique(key))
design1 <- model.matrix(~0 + key)
colnames(design1) <- levels(key)


# limma and ebayes

fit <- lmFit(eset2, design1,method = "robust")
contrast.matrix <- makeContrasts("control-E2", "E2-control", levels=design1)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2,trend = TRUE)


outputfit2NR2 <- topTable(fit2, adjust="BH", coef=1, genelist=rownames(eset2), number=nrow(fit2$coefficients))
outputfit2NR<-outputfit2NR2[abs(outputfit2NR2$logFC) > 1.5,]



