### Prepare the RNA-Seq table ###
### This is just a standard processing of file
de <- read.csv("~/Documents/GitHub/RBPInper/inst/Joe/sfpq-kd_pvalues.csv")
de <- de[1:7]
de <- de[c(1, 7, 2:6)]
names(de) <- c("gene_id", "symbol", "GSE149370", "GSE157622", "GSE157622", "ENCSR782MXN", "ENCSR535YPK")
de[is.na(de)] <- 1
### Save to file
write.csv(de, file = "rna-seq.csv", row.names = F)

#### Sort the lack of p-values in encode narrow peak files ####
bed <- read.delim("~/Desktop/data_clean/ENCFF923EYZ.bed", header=FALSE)
premac2bed(ENCFF073RTN, nam="ENCFF923EYZ.bed", soure="encode")

### Process the bed files ###
path <- "/Users/chineduanene/Desktop/data_clean/"
### List files with full path
listpath <- as.vector(list.files(path, pattern = "\\.bed$", full.names = TRUE))
### Set gtf path
gtf <- "/Users/chineduanene/Desktop/Ens.GRCh38.99.gtf"
processed <- prebed(bed=listpath, gtf=gtf)
names(processed)[2] <- "symbol"
### Load the RNA-Seq knockdown data
rna.seq <- read.csv("~/Desktop/data_clean/rna-seq.csv")
# Set the id and symbol
idsym <- rbind(processed[1:2], rna.seq[1:2])
idsym <- idsym[!duplicated(idsym$gene_id),]

final <- merge(processed[-2], rna.seq[-2], by="gene_id",
               all.x=T, all.y=T)
final[is.na(final)] <- 1
final <- merge(idsym, final, by="gene_id", all.x=F,
               all.y=F)
write.csv(final, file="pre_integration.csv", row.names = F)


### Run the final analysis
predata <- read.csv("~/Documents/GitHub/RBPInper/inst/data/pre_integration.csv")
infodata <- read.csv("~/Documents/GitHub/RBPInper/inst/data/info.csv")

infoRNA <- infodata[infodata$method %in% c("eCLIP", "RIP-Seq", "RNA-Seq"), ]
infoDNA <- infodata[infodata$method %in% c("ChIP", "RNA-Seq"), ]

resultRNA <- RBPInper::rbpinper.run(evi=predata, info = infoRNA)
resultDNA <- RBPInper::rbpinper.run(evi=predata, info = infoDNA)


## Get and count the interactors
mergdf <- cbind(resultRNA@L1.result, resultRNA@L2.result)
this <- names(mergdf)[1:6]
for(i in this){
  mergdf[i] <- ifelse(mergdf[[i]] <= 0.05, "Hit", "")
}
mergdf <- mergdf[c(1:6, 10)]
names(mergdf)[7] <- "Global"
mergdf <- colSums(mergdf == "Hit")
barplot(mergdf, names.arg = colnames(df),
        col = c(rep("skyblue", 6), "darkred"),
        main = "", cex.names = 0.7, las = 2,
        ylim = c(0, 10000),
        xlab = "Interactome", ylab = "Hit Count")
#
## Get target column
resultRNA <- resultRNA@L2.result
resultDNA <- resultDNA@L2.result

#
resultRNA <- resultRNA[resultRNA$call == "Hit", ]
resultDNA <- resultDNA[resultDNA$call == "Hit", ]

resultRNA <- merge(predata[1:2], resultRNA, by="gene_id")
resultDNA <- merge(predata[1:2], resultDNA, by="gene_id")

check <- merge(resultRNA, resultDNA, by = "gene_id", all.y = T)
check[is.na(check)] <- "dna"
table(check$call.x)
check <- check[check$call.x == "dna", ]
check <- check[c(1, 6:9)]
names(check) <- c("gene_id", "symbol", "global", "adjP", "call")

## Save as supplementary
write.csv(resultRNA, file = "SFPQ_RNA.csv", row.names = F)
write.csv(resultDNA, file = "SFPQ_DNA.csv", row.names = F)
write.csv(check, file = "Exclusive_SFPQ_DNA.csv", row.names = F)

# Gene ontology on the RNA interaction
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(enrichplot)
library(DOSE)

# Ontology single
ego <- enrichGO(gene          = resultRNA$symbol,
                universe      = predata$symbol,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01,
                readable      = T)
##
d <- godata('org.Hs.eg.db', ont="BP")
ego2 <- pairwise_termsim(ego, method = "Wang", semData = d)
#
pdf(width = 13, height = 8)
treeplot(ego2, showCategory = 30)
dev.off()
#  heatplot(ego,  showCategory = 2, symbol="dot", label_format=1)
ontology <- ego@result
write.csv(ontology, file = "biological processes_RNA.csv", row.names = F)

### Get interesting target to validate
validate <- res2@compareClusterResult
validate <- ontology[ontology$Description %in% c("Wnt signaling pathway",
                                                 "heterochromatin formation"), ][["geneID"]]

validate <- c(strsplit(validate[1], split = "/"),
           strsplit(validate[2], split = "/"))
#
valiframe1 <- data.frame("gene_name" = validate[[1]])
valiframe1$onto <- "Wnt signaling pathway"
#
valiframe2 <- data.frame("gene_name" = validate[[2]])
valiframe2$onto <- "heterochromatin formation"

# Bind full
validate <- rbind(valiframe1, valiframe2)
write.csv(validate, file = "Validation_target_SFPQ_RNA.csv", row.names = F)

## Load the PAR-CLIP independent validation
U2OS <- read.csv("~/Documents/GitHub/RBPInper/inst/data/independent validation/GSM3104133_PARCLIP_U2OS.csv")
HeLa <- read.csv("~/Documents/GitHub/RBPInper/inst/data/independent validation/GSM3104135_PARCLIP_HeLa.csv")

movbase <- resultRNA[c(1, 2)]
movbase$RBPInper <- "RBPInper"

U2OS <- U2OS[c(17,16)]
U2OS <- U2OS[!duplicated(U2OS$gene_id),]
U2OS$Parclip <- "Parclip"
U2OS <- merge(U2OS, movbase, by="gene_id", all.x = T,
              all.y = T)
# Fill the NA symbols
U2OS$symbol <- ifelse(is.na(U2OS$symbol.x),
                      U2OS$symbol.y, U2OS$symbol.x)



HeLa <- HeLa[c(17,16)]
HeLa <- HeLa[!duplicated(HeLa$gene_id),]
HeLa$Parclip <- "Parclip"
HeLa <- merge(HeLa, movbase, by="gene_id", all.x = T,
              all.y = T)
# Fill the NA symbols
HeLa$symbol <- ifelse(is.na(HeLa$symbol.x),
                      HeLa$symbol.y, HeLa$symbol.x)

## Subset the columns
U2OS <- U2OS[c(1, 6, 5, 3)]
HeLa <- HeLa[c(1, 6, 5, 3)]
## Add the cell group
U2OS$cell <- "U2OS"
HeLa$cell <- "HeLa"
## Create a group for plot
parvalid <- rbind(U2OS, HeLa)
parvalid$group <- paste(parvalid$RBPInper, parvalid$Parclip,
                        sep = "_")
parvalid$group <- gsub("_NA", "", parvalid$group)
parvalid$group <- gsub("NA_", "", parvalid$group)

## Get the tables for each
parv1 <- parvalid[!parvalid$group %in% c("RBPInper"),]
parv1$test <- ifelse(parv1$group == "Parclip", "out", "in")
parU <- parv1[parv1$cell == "U2OS", ]
parH <- parv1[parv1$cell == "HeLa", ]

tabbU <- table(parU$cell, parU$test)
tabbH <- table(parH$cell, parH$test)
prop.table(tabbU)
prop.table(tabbH)
tabb <- table(parv1$cell, parv1$test)

Chiq1 <- chisq.test(tabbU)
Chiq2 <- chisq.test(tabbH)
#
propT <- prop.table(tabb, 1)
#
mosaicplot(propT, shade = F,
           las = 1, off = 6,
           main = "",
           sub = "p 1.576e-64 | 5.302e-74",
           ylab = "RBPInper",
           xlab = "PAR-CLIP",
           type = "pearson",
           color = c("darkred", "gray"))


########### Start of large supplementary table binding ###########################
# Classify eclip directly because of method
predata2 <- resultRNA@meta

predata2$ENCFF417EZT <- ifelse(predata2$ENCFF417EZT <= 0.05, "Hit", "")
predata2$ENCFF960OTE <- ifelse(predata2$ENCFF960OTE <= 0.05, "Hit", "")

# Do the RNA-Seq with padjust
this2 <- c("Omera", "Omera1", "GSE149370", "GSE157622", "GSE157622.1", "ENCSR782MXN", "ENCSR535YPK")

for(i in this2){
  predata2[i] <- p.adjust(predata2[[i]])
  predata2[i] <- ifelse(predata2[[i]] <= 0.05, "Hit", "")
}

predata2 <- predata2[c("gene_id", "symbol", "ENCFF417EZT", "ENCFF960OTE", "Omera",
                       "Omera1", "GSE149370", "GSE157622",
                       "GSE157622.1", "ENCSR782MXN", "ENCSR535YPK")]
## add the global
predata2$global <- resultRNA@L2.result$call


## Next process the independent data
U2OS <- U2OS[c(17,16)]
U2OS <- U2OS[!duplicated(U2OS$gene_id),]
U2OS$U2OS <- "Hit"

HeLa <- HeLa[c(17,16)]
HeLa <- HeLa[!duplicated(HeLa$gene_id),]
HeLa$HeLa <- "Hit"

extern <- merge(U2OS, HeLa,  by="gene_id", all.x = T,
                all.y = T)

# Fill the NA symbols
extern$symbol <- ifelse(is.na(extern$symbol.x),
                        extern$symbol.y, extern$symbol.x)
extern <- extern[c("gene_id","symbol", "U2OS", "HeLa")]
extern[is.na(extern)] <- ""

## Merge with the integrated
predata3 <- merge(predata2, extern, by="gene_id",
                  all.x = T, all.y = T)
# Remove one line with missing Esemble id odd
predata3 <- predata3[!predata3$gene_id == "", ]

# Fill in missing symbol again
predata3$symbol <- ifelse(is.na(predata3$symbol.x),
                          predata3$symbol.y, predata3$symbol.x)
predata3 <- predata3[c("gene_id", "symbol", "ENCFF417EZT", "ENCFF960OTE",
                       "Omera",   "Omera1",   "GSE149370", "GSE157622",
                       "GSE157622.1", "ENCSR782MXN", "ENCSR535YPK", "global", "U2OS",
                       "HeLa" )]

# Next load and add the literature evidence
lit <- read.csv("~/Documents/GitHub/RBPInper/inst/data/known_interactors_lit.csv")
lit <- lit[ "ensemblID"]
names(lit) <- "gene_id"
lit$Lit <- "Hit"
#### Merge it in
predata3 <- merge(predata3, lit, by="gene_id",
                  all.x = T, all.y = T)
# Change NA to space
predata3[is.na(predata3)] <- ""

# save to file
write.csv(predata3, file = "meta_results_indviduals.csv", row.names = F)
########### End of large supplementary table binding ###########################


## Now compute hit overlap
# U2OS
per_hit <- sapply(predata3[3:12], function(col) sum(col == "Hit" &
                                                predata3[["U2OS"]] == "Hit") /
                    sum(predata3[["U2OS"]] == "Hit") * 100)

# Create a new data frame with the results
U2OS_ovlp <- data.frame(Column = names(per_hit), overlap = per_hit)
U2OS_ovlp$group <- "U2OS"

# HeLa
per_hit2 <- sapply(predata3[3:12], function(col) sum(col == "Hit" &
                                                      predata3[["HeLa"]] == "Hit") /
                    sum(predata3[["HeLa"]] == "Hit") * 100)

# Create a new data frame with the results
HeLa_ovlp <- data.frame(Column = names(per_hit2), overlap = per_hit2)
HeLa_ovlp$group <- "HeLa"

## Benchmark
bench <- rbind(U2OS_ovlp, HeLa_ovlp)
bench$Column <- gsub("global", "RBPInper", bench$Column)
bench <- bench[order(bench$overlap, bench$group), ]

# Convert X_Labels to a factor with custom order
bench$Column <- factor(bench$Column, levels = c("RBPInper", "Omera", "Omera1",
                                                "ENCFF417EZT", "ENCFF960OTE",
                                                "GSE157622", "GSE157622.1",
                                                "ENCSR535YPK", "ENCSR782MXN",
                                                "GSE149370"))

# Plot the overlap
# Load the ggplot2 package
library(ggplot2)

# Create a grouped bar plot
cus_col <- c("HeLa" = "darkblue", "U2OS" = "darkred")

ggplot(bench, aes(x = Column, y = overlap, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = cus_col) +
  labs(title = "",
       y = "Percentage Overlap",
       x = "Interactome",
       fill = "group") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## Get aggreagte for reporting
aggregate(overlap ~ Column, data = bench, FUN = mean)



########
# Automatic literature curation for SFPQ
res_all_genes$Gene_Words <- apply(res_all_genes, 1, function(column) {
  extract_genes(column["Title"], pt=c("/", "-", ":", "·", ",", '"'),
                remove = c("RNA", "DNA", "SFPQ", "PSF", "3D"))
})

#creating a table with frequency of filtered words
gene_table <- table(unlist(strsplit(res_all_genes$Gene_Words, split = " ")))
unique_gene_table <- data.frame(symbol = names(gene_table), Freq = as.vector(gene_table))

#_genes to ensemblID___====
unique_gene_table$ensemblID <- mapIds(org.Hs.eg.db,
                                      keys = unique_gene_table$symbol,
                                      keytype = "SYMBOL",
                                      column = "ENSEMBL")

lit_genes <- unique_gene_table[!is.na(unique_gene_table$ensemblID),]

write.csv(lit_genes, file="lit_extract_genes.csv", row.names = F)
