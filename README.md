# DE_Genes_and_ERVs_together
R code to perform differential expression of genes and ERVs together, but normalised according to genes only.  If you are working with genes only, the first part of the script can be run independently. This follows after the https://github.com/rkabiljo/RNASeq_Genes_ERVs pipeline.

## Dependencies

<br>1.phenotype.txt: the phenotype table with all the variables needed for the analysis.  In our case: Status, Sex, Age, PMD, RIN, Site_Specimen_Collected
<br>2.merged_cellular.txt: matrix with raw merged cellular counts,in our case coming from HTseq
<br>3.merged_erv.txt: matrix with raw merged ERV counts,in our case coming from ERVMap 

## Load Phenotype Data, categorise variables
<br> if it is the first time, derive the categories needed.  If the table with categories, sizeFactor column and SV1 exists, skip to that part below <br>

```
#in this table we need columns Age, PMD, Sex,Status, RIN...  Adjust if the names are different

phe <- read.table("phenotype.txt",row.names=1,header=T) 
#categorise these as DEseq needs that
phe$AgeCat <- cut(phe$Age, 5, labels = c('1','2','3','4','5'),include.lowest = TRUE, right = TRUE, dig.lab = 3,ordered_result = FALSE)
phe$PMDCat <- cut(phe$PMD, 5, labels = c('1','2','3','4','5'),include.lowest = TRUE, right = TRUE, dig.lab = 3,ordered_result = FALSE)
phe$RINCat <- cut(phe$RIN, 5, labels = c('1','2','3','4','5'),include.lowest = TRUE, right = TRUE, dig.lab = 3,ordered_result = FALSE)

#these need to be factors for DEseq
phe$Sex<- as.factor(phe$Sex)
phe$Status<- as.factor(phe$Status)
phe$Sex<- as.factor(phe$Sex)
phe$AgeCat<- as.factor(phe$AgeCat)
phe$PMDCat<- as.factor(phe$PMDCat)
phe$RINCat <- as.factor(phe$RINCat)

#for TargetALS only, there is also Site_Specimen_Collected
#phe$Site_Specimen_Collected <- as.factor(phe$Site_Specimen_Collected)

#set status 1 as the baseline, in my case these are controls
library("magrittr")
phe$Status %<>% relevel("1")
```

## Normalise, Filter, Derive SV1 and SizeFactors from gene counts

### Read gene counts, normalise
```
library(DESeq2)
gene_cts <- read.table("geneMatrixEdited.txt",row.names=1,header=T)

phe_g<-phe
phe_g<-phe_g[colnames(gene_cts),]
dds_genes <- DESeqDataSetFromMatrix(countData = gene_cts, colData = phe_g , design = ~ Sex + AgeCat + PMDCat + RINCat + Status)

#In TargetALS there is also site where the secimen was collected
#dds_genes <- DESeqDataSetFromMatrix(countData = gene_cts, colData = phe_g , design = ~ Sex + AgeCat + PMDCat + RINCat + Site_Specimen_Collected + Status)

#this estimates size factors and performs normalisation according to these
dds_genes <- estimateSizeFactors(dds_genes)
```

### filter rows with low counts
```
idx <- rowSums( counts(dds_genes, normalized=T) >= 5 ) >= 10
dds_genes <- dds_genes[idx,]
```

### derive surrogate variables

```
library("sva")
dat  <- counts(dds_genes, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ Sex + AgeCat + PMDCat + RINCat + Status, colData(dds_genes))
#for TargetALS add site
# mod  <- model.matrix(~ Sex + AgeCat + PMDCat + RINCat + + Site_Specimen_Collected  + Status, colData(dds_genes))

mod0 <- model.matrix(~   1, colData(dds_genes))
nsv <- num.sv(dat,mod)
nsv
svseq <- svaseq(dat, mod, mod0, n.sv = nsv)
dds_genes$SV1 <- svseq$sv
#if nsv is more than 1, adjust that and use as many as there are

```

### update the design to use these surrogate variables

```
design(dds_genes) <- ~ Sex + AgeCat + PMDCat + RINCat + SV1 + Status
#for TargetALS add site
design(dds_genes) <- ~ Sex + AgeCat + PMDCat + RINCat + Site_Specimen_Collected + SV1 + Status
```

### write the table with everything we derived so far, so that the next time we can start from that table
```
write.table(colData(dds_genes), file = "samples.design.updated.txt",quote=FALSE,sep="\t")

```

# If the table with SV1 and sizeFactor exist, start from HERE

## Read the table saved - as above

```
phe <- read.table("samples.design.updated.txt",row.names=1,header=T,check.names=FALSE)
```

### Prepare the data: factorise, relevel

```
phe$Sex<- as.factor(phe$Sex)
phe$Status<- as.factor(phe$Status)
phe$Sex<- as.factor(phe$Sex)
phe$AgeCat<- as.factor(phe$AgeCat)
phe$PMDCat<- as.factor(phe$PMDCat)
phe$RINCat <- as.factor(phe$RINCat)
#for TargetALS only, there is also Site_Specimen_Collected
#phe$Site_Specimen_Collected <- as.factor(phe$Site_Specimen_Collected)

#set status 1 as the baseline, in my case these are controls
library("magrittr")
phe$Status %<>% relevel("1")
```
### Read in the gene matrix, define design and filter as before <br>

```
gene_cts <- read.table("merged_cellular.txt",row.names=1,header=T)

phe_g<-phe
phe_g<-phe_g[colnames(gene_cts),]
dds_genes <- DESeqDataSetFromMatrix(countData = gene_cts, colData = phe_g , design = ~ Sex + AgeCat + PMDCat + RINCat + SV1 + Status)

#In TargetALS there is also site where the secimen was collected
#dds_genes <- DESeqDataSetFromMatrix(countData = gene_cts, colData = phe_g , design = ~ Sex + AgeCat + PMDCat + RINCat + Site_Specimen_Collected + SV1 + Status)

#filter our rows with low counts
idx <- rowSums( counts(dds_genes, normalized=T) >= 5 ) >= 10
dds_genes <- dds_genes[idx,]

```

## Run DE - genes only
```
library("IHW")
dds_genes <- DESeq(dds_genes)
resIHW <- results(dds_genes, filterFun=ihw)
resIHWOrdered <- resIHW[order(resIHW$pvalue),]
#explore
sum(resIHW$padj < 0.05, na.rm=TRUE)
head(resIHWOrdered,10)

write.table(resIHWOrdered, "genesDEres.txt",sep="\t")
```
<b>If you are working only with genes, stop here :)

# Run DE - merge genes and ERVs

<br>This assumes that phenotype data (phe) is already processed: factorised, and SV1 and SizeFactors are in that phe table. See Section 'If the table with SV1 and sizeFactor exist, start from HERE' <br>

```
gene_cts <- read.table("merged_cellular.txt",row.names=1,header=T,check.names=FALSE)
gene_cts<-gene_cts[,rownames(phe)]
#this is raw erv counts, NOT the ones normalised in ERVmap script
cts <- read.table("merged_erv.txt",row.names=1,header=T,check.names=FALSE)
cts<-cts[,rownames(phe)]
gene_erv_cts<-rbind(gene_cts,cts)

#sanity check
dim(gene_erv_cts)
dim(gene_cts)
dim(cts)

dds_genes <- DESeqDataSetFromMatrix(countData = gene_erv_cts, colData = phe , design = ~ Sex + AgeCat + PMDCat + RINCat + SV1 + Status)
#target
#dds_genes <- DESeqDataSetFromMatrix(countData = gene_erv_cts, colData = phe , design = ~ Sex + AgeCat + PMDCat + RINCat + Site_Specimen_Collected + SV1 + Status)
```

### filter low counts
```
idx <- rowSums( counts(dds_genes, normalized=T) >= 5 ) >= 10
dds_genes <- dds_genes[idx,]
dds_genes
```

### run DE on the merged genes and ervs
```
dds_genes <- DESeq(dds_genes)
resIHW <- results(dds_genes, filterFun=ihw)
resIHWOrdered <- resIHW[order(resIHW$pvalue),]
```
### explore and write results
```
sum(resIHW$padj < 0.05, na.rm=TRUE)
sum(resIHW$padj < 0.05 & abs(resIHW$log2FoldChange)>0.2, na.rm=TRUE)
write.table(resIHWOrdered, "genes_ervs_DE_res.txt",sep="\t")

```
If needed, separate gene and ERV results.  In my example, every line which does not have 'ENSG' will be an ERV line and I am saving these separately
```

library(tidyverse)
erv_res <- as.data.frame(resIHWOrdered) %>%
  filter(!grepl("ENSG", rownames(resIHWOrdered)))
write.table(erv_res, "ervs_DE_res.txt",sep="\t")
```
