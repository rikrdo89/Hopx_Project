# 1. Importing Data --------------
library(GEOquery)

#2010 Dev Bio Trivedi et al
geo = "GSE23700" # chnage the dataset to be downloaded
mydir = paste("./", geo, collapse = "") 
if(!dir.exists(mydir)) dir.create(mydir) # creates a dir with the neame of the dataset if it doesnt exist

GSE23700 <- getGEO(geo, destdir = mydir)

GSE23700
class(GSE23700)
length(GSE23700)

# correctExp <- function(X) {
#   if(class(X) == "list" & length(X) ==1 ) X <- X[[1]]
# }
# correctExp(GSE23700)

class(GSE23700)
length(GSE23700)

GSE23700 <- GSE23700[[1]]
class(GSE23700)
length(GSE23700)

# 2. Assigning phenoData, featureData, and expressionDaata-------------
pd <- pData(GSE23700)
fd <- fData(GSE23700)
eset <- exprs(GSE23700)


# 3. Filtering GSMs based on scientific question
#toMatch <- c("wild-type","hopx-null")
ls <- grep("wild-type|hopx-null", pd$title, ignore.case = TRUE)
View(ls)

groupGenotype <- rep(c("wt", "null"), each = 3 )
?rep
View(groupGenotype)
colnames(GSE23700)
NewGSE <- GSE23700[,ls]
colnames(NewGSE)
class(NewGSE)
NewGSE

# 4. Assigning phenoData, featureData, and expressionDaata based on filtered DATA -------------
pd <- pData(NewGSE)
fd <- fData(NewGSE)
eset <- exprs(NewGSE)

# 5. Data distribution, is it normailized? -----------

boxplot(eset, outline=FALSE)
# log2 convertion
esetlog2 <- log2(exprs(NewGSE))
boxplot(esetlog2, outline=FALSE)
View(fd)



ls2 <- grep("hopx", fd$`Gene Symbol`, ignore.case = TRUE)
boxplot(esetlog2[ls2[1],]~groupGenotype, main = "Hopx Expresion by Array")

ls2
fd[ls2,c(1,2,3,10,11)]
colnames(fd)
plot(esetlog2[ls2[2],], main = "Hopx_1")

# 6. FIlter Top 75% -------------
exprs(NewGSE) <- log2(exprs(NewGSE))
library(genefilter)
FilteredGSE <- varFilter(NewGSE, var.cutoff = 0.75) # FIlter only top 50% most variable genes across samples.
class(FilteredGSE)
nrow(FilteredGSE)
nrow(FilteredGSE)/nrow(NewGSE)
# manking new fetures and phenotype dbs
esetF <- exprs(FilteredGSE)
pd <- pData(FilteredGSE)
fd <- fData(FilteredGSE)

View(esetF[c(3462,9051),])
nrow(fd)


# 8. Cluster --------------

# 9. PCA analysis ------------
GroupColors <- rep(c("red", "blue"), each = 3 )

pca <- prcomp(t(exprs(FilteredGSE)))
class(pca)
head(pca$x)
plot(pca)
plot(pca$x[,1], pca$x[,2], xlab = "PCA1", ylab = "PCA2", pch=16, col = as.character(GroupColors), cex =2)
legend("topleft", fill = c("red", "blue"), legend =c("wt", "Hopx-null"))
?plot

?legend
groupGenotype
# 10. DE analysis using model.matrix ------------
library(limma)
design <- model.matrix(groupGenotype)
View(design)
class(design)
colnames(design) <- c("Hopx_null", "WT")

fit <- lmFit(esetF, design) # lmFit from limma package!!!
names(fit)

head(fit$coefficients)

# defining contrast matrix
contrast.matrix <- makeContrasts(Hopx_null-WT, levels = design)
?makeContrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
names(fit2)
fit2$coefficients[c(3462,9051),]
fit2.bayes <- eBayes(fit2)

names(fit2.bayes)
fit2.bayes$coefficients
topTable(fit2.bayes)
?topTable

DA_results <- topTable(fit2.bayes, number = Inf)
class(DA_results)
View(DA_results)

# 11. Adding annotation to DE table -------
colnames(fd)
anno <- fd[, c("Gene Title", "Gene Symbol", "ENTREZ_GENE_ID")]
View(anno)
fit2.bayes$genes <- anno
names(fit2.bayes)
class(fit2.bayes)
View(fit2.bayes)


decideTests(fit2.bayes) #, adjust.method = "none")

?decideTests
table(decideTests(fit2.bayes))
sum(abs(decideTests(fit2.bayes))==1)

# 12. Volcanoplot -------------
volcanoplot(fit2.bayes, highlight = 15, names = fit2.bayes$genes$`Gene Symbol`, main = "Hopx_Null - WT")
?volcanoplot
names(fit2.bayes)
head(fit2.bayes$coefficients)

write.fit(fit2.bayes, file = "de-results.txt",adjust="BH")

DE_results2 <- read.delim("de-results.txt")
View(DE_results2)


?write
# Annotation --------
BiocManager::install("org.Mm.eg.db", version = "3.8")
library(org.Mm.eg.db)
columns(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
genedb <- select(org.Mm.eg.db, keys = c("Hopx", "Tp53", "Brca1", "Pten"), keytype = "SYMBOL", columns = c("ENTREZID", "UNIGENE", "GENENAME", "MGI", "UNIPROT")) #case sensitive!
class(genedb)
       
