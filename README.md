# DE_Genes_and_ERVs_together
this performs differential expression of genes and ERVs together, but normalised according to genes only

<br>
# Load Phenotype Data - step1, categorise variables
<br> if it is the first time, derive the categories needed.  If the table with categories, sizeFactor column and SV1 exists, skip to that part below <br>

```
#in this table we need columns Age, PMD, Sex,Status, RIN...  Adjust if the names are different

phe <- read.table("phenotype.txt",row.names=1,header=T) 
#categorise these as DEseq needs that
phe$AgeCat <- cut(phe$Age, 5, labels = c('1','2','3','4','5'),include.lowest = TRUE, right = TRUE, dig.lab = 3,ordered_result = FALSE)
phe$PMDCat <- cut(phe$PMD, 5, labels = c('1','2','3','4','5'),include.lowest = TRUE, right = TRUE, dig.lab = 3,ordered_result = FALSE)

#these need to be factors for DEseq
phe$Sex<- as.factor(phe$Sex)
phe$Status<- as.factor(phe$Status)
phe$Sex<- as.factor(phe$Sex)
phe$AgeCat<- as.factor(phe$AgeCat)
phe$PMDCat<- as.factor(phe$PMDCat)

#set status 1 as the baseline, in my case these are controls
library("magrittr")
phe$Status %<>% relevel("1")
```
<br>
# Derive SV1 and SizeFactors

<br>


<br>
if the table exist, load it
```
code
```

# Run DE
```
code
```

#If needed, separate gene and ERV results
```
code
```
