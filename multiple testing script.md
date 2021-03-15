##Analysis of the effect of dexamethasone on gene set expression

#1 look at the data structure
```{r}
BiocManager::install("GEOquery")
library(GEOquery)
g <- getGEO("GSE34313")
e <- g[[1]] 
dim(e)
```

```{r}
e$condition = e$characteristics_ch1.2
e$condition=factor(e$condition)
levels(e$condition)=c("dex24","dex4","cont") #re-label 
table(e$condition)
```

```{r}
library(rafalib)
mypar(1,2)
boxplot(exprs(e),range=0) # the data is normalized 

```


#2 Test the expression of gene set of immune response between control and 4 hours after dexamethasone use

```{r}
lvls=c("cont","dex4") # indicator
es=e[,e$condition%in%lvls] 
```

```{r}
library(limma)
design <- model.matrix(~ es$condition) 
fit <- lmFit(es, design=design) # linear model for statistical testing
fit <- eBayes(fit) 
tt <- topTable(fit, coef=2, genelist=fData(es)$GENE_SYMBOL)
tt # shows that CSF2, LIF, CCL2..., have negative fold change 4 hours after dexamethasone use

```


#3 test single gene set expression diffrence between 2 groups
```{r}
idx <- grep("GO:0006955", fData(es)$GO_ID) #GO:0006955:immune response
r1 <- roast(es, idx, design)
r1 # down & mixed up/down are significantly different between 2 groups

```


#4 test multiple gene sets
```{r}
library(org.Hs.eg.db)
go2eg <- as.list(org.Hs.egGO2EG)# human genes by GO as a list 
govector=unlist(go2eg)# all the human GO as characters

```

```{r}
golengths=sapply(go2eg,length)
idxvector <- match(govector, fData(es)$GENE) # which one is in es data
table(is.na(idxvector)) # see the matched one
idx <- split(idxvector,rep(names(go2eg), golengths)) # match GO:00XXXX with genes

```

```{r}
# filter the gene sets with less than 10 genes to prevent significant results by small number of genes 

idxclean <- lapply(idx, function(x) x[!is.na(x)])# which id is matched
idxlengths <- sapply(idxclean, length) # length of the matched gene set
idxsub <- idxclean[idxlengths > 10] 
```

```{r}
r2 <- mroast(es, idxsub, design)
r2 <- r2[order(r2$PValue.Mixed),] 
r2
```
