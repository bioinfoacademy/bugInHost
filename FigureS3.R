library(readxl)
library(ggpubr)
library(matrixStats)
library(GGally)
library(ggbiplot)
library(cowplot)
library(janitor)
library(ggrepel)
library(biomaRt)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)

library(tidyverse)
#library(org.Mm.eg.db)
library(refGenome)


my_fn <- function(data, mapping, method="p", use="pairwise", ...){
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}


files_counts <- sort(list.files(path="~/PATH/", pattern="BF.*.counts"))
files_coverage <- sort(list.files(path="~/PATH/", pattern="BF.*.coverage"))

metadata <- tibble(sampleID = colnames(readCountTable_all), mouseID = rep(c(1:3), each=6), site = rep(c("lumen","mucus","tissue"), each=2, times=3), enrichment = rep(c("nonHS","HS"),times=9))


# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

# read GTF file into ensemblGenome object
read.gtf(ens, "./data/reference/Mus_musculus.GRCm38.82.chr.gtf")
tableSeqids(ens)

# create table of genes
mm10_gene <- getGenePositions(ens)
dim(mm10_gene)
head(mm10_gene)
mm10_gene[grep("0610005C13Rik", mm10_gene$gene_name),"end"] - mm10_gene[grep("0610005C13Rik", mm10_gene$gene_name),"start"]
mm10_gene$geneLength <- mm10_gene$end - mm10_gene$start

readCoverageTable_all <- tibble()
i=1
readCoverageTable_one <- read_tsv(paste0("~/PATH/",files_coverage[i]), col_names = F)
readCoverageTable_all <- readCoverageTable_one[,8]
for(i in 2:length(files_coverage)){
  readCoverageTable_one <- read_tsv(paste0("~/PATH/",files_coverage[i]), col_names = F)
  readCoverageTable_all <- bind_cols(readCoverageTable_all, readCoverageTable_one[,8])
}
colnames(readCoverageTable_all) <- sub("\\.host.*$","",files_coverage)
readCoverageTable_name_length <- readCoverageTable_one[,c(4,7)]


readCountTable_all <- tibble()
i=1
readCountTable_one <- read_tsv(paste0("~/PATH/",files_counts[i]), col_names = F)
readCountTable_all <- readCountTable_one[,2]
for(i in 2:length(files_counts)){
  readCountTable_one <- read_tsv(paste0("~/PATH/",files_counts[i]), col_names = F)
  readCountTable_all <- bind_cols(readCountTable_all, readCountTable_one[,2])
}
colnames(readCountTable_all) <- gsub("\\.host.*$","",files_counts)
readCountTable_summary <- readCountTable_all[grepl("^__",readCountTable_one$X1),]
readCountTable_summary_rowName <- readCountTable_one[grepl("^__",readCountTable_one$X1),1]

readCountTable_all <- readCountTable_all[!grepl("^__",readCountTable_one$X1),]
readCountTable_all_rowName <- readCountTable_one[!grepl("^__",readCountTable_one$X1),1]

#plot(density(rowMeans(readCountTable_all)), xlim=c(0,3000))


#plot total read count of each sample
#data.frame(colSum=colSums(readCountTable_all)) %>% tibble::rownames_to_column("sampleID") %>% ggplot(aes(x=sampleID, y=colSum)) + geom_label_repel(aes(label=sampleID)) + geom_point()

# ggpairs(log10(readCountTable_all+1) %>% sample_n(100, replace=F), 
#         upper = list(continuous = my_fn),
#         lower = list(continuous = "smooth"))


readCountTable_name_length <- readCoverageTable_name_length[match(readCountTable_all_rowName$X1,readCoverageTable_name_length$X4),]
colnames(readCountTable_name_length) <- c("name","length")

readCountTable_all_name <- bind_cols(readCountTable_name_length %>% dplyr::select(-length), readCountTable_all)
readCountTable_all_name_melt <- gather(readCountTable_all_name, sampleID, readCount, -name)

readCountTable_all_name_melt_metadata <- left_join(readCountTable_all_name_melt, metadata, by=c("sampleID"="sampleID"))

readCountTable_all_name_melt_metadata_averageReadCount <- readCountTable_all_name_melt_metadata %>% dplyr::select(-mouseID, -sampleID) %>% group_by(name, site, enrichment) %>% dplyr::summarise(averageReadCount = mean(readCount)) %>% ungroup() %>% mutate(siteEnrichment=paste0(site,"_",enrichment)) %>% dplyr::select(-site, -enrichment) %>% spread(siteEnrichment, averageReadCount)

readCountTable_all_name_melt_metadata_averageReadCount
readCountTable_all_name_melt_metadata_averageReadCount2 <- left_join(readCountTable_all_name_melt_metadata_averageReadCount, mm10_gene %>% as_tibble(), by=c("name"="gene_name")) 
#remove MT genes
readCountTable_all_name_melt_metadata_averageReadCount2 <- readCountTable_all_name_melt_metadata_averageReadCount2 %>% filter(!grepl("MT",seqid))


# ggpairs((log10(readCountTable_all_name_melt_metadata_averageReadCount2 %>% dplyr::select(contains("HS"))+1) %>% sample_n(1000, replace=F)) , 
#         #upper = list(continuous = my_fn),
#         upper=list(continuous = wrap("cor", size = 6)),
#         lower = list(continuous = "smooth")) #

## with all 46403 genes
# ggpairs((log10(readCountTable_all_name_melt_metadata_averageReadCount[,-1]+1)) , 
#         #upper = list(continuous = my_fn),
#         upper=list(continuous = wrap("cor", size = 6)),
#         lower = list(continuous = "smooth")) #

## Lumen
dataForPlot <-  readCountTable_all_name_melt_metadata_averageReadCount2 %>% dplyr::select(lumen_HS,lumen_nonHS,name,seqid,geneLength,gene_biotype)
colnames(dataForPlot)[1:2] <- gsub("^.*_","",colnames(dataForPlot)[1:2])
dataForPlot_SD <- dataForPlot[which(abs(log10(dataForPlot$HS+1)-log10(dataForPlot$nonHS+1)) > 4*sd(log10(dataForPlot$HS+1)-log10(dataForPlot$nonHS+1))),]
nrow(dataForPlot_SD)/nrow(dataForPlot)*100

ggplot(data=dataForPlot, aes(x=log10(nonHS+1),y=log10(HS+1))) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(readCount+1)\nnonHS", y="HS\nlog10(readCount+1)")+geom_abline(intercept =0 , slope = 1, color = "gray")+geom_point(data=dataForPlot_SD, aes(x=log10(nonHS+1),y=log10(HS+1)),alpha=1, color="blue") +geom_point(data=dataForPlot %>% filter(nonHS==0 | HS==0), aes(x=log10(nonHS+1),y=log10(HS+1)),alpha=1, color="red") 

geneOnAxisX <- dataForPlot %>% dplyr::filter(nonHS!=0 & HS==0)
nrow(geneOnAxisX)/nrow(dataForPlot)*100
sort(table(geneOnAxisX$gene_biotype))
nrow(geneOnAxisX %>% filter(gene_biotype=="protein_coding"))/nrow(dataForPlot)*100
nrow(geneOnAxisX %>% filter(gene_biotype!="protein_coding"))/nrow(dataForPlot)*100
geneOnAxisY <- dataForPlot %>% dplyr::filter(nonHS==0 & HS!=0)
nrow(geneOnAxisY)/nrow(dataForPlot)*100
nrow(geneOnAxisY %>% filter(gene_biotype=="protein_coding"))/nrow(dataForPlot)*100
nrow(geneOnAxisY %>% filter(gene_biotype!="protein_coding"))/nrow(dataForPlot)*100
geneOnAxisXY <- dataForPlot %>% dplyr::filter(nonHS==0 & HS==0)
nrow(geneOnAxisXY)/nrow(dataForPlot)*100
nrow(geneOnAxisXY %>% filter(gene_biotype=="protein_coding"))/nrow(dataForPlot)*100
nrow(geneOnAxisXY %>% filter(gene_biotype!="protein_coding"))/nrow(dataForPlot)*100

ggplot(data=dataForPlot_SD, aes(x=log10(nonHS+1),y=log10(HS+1))) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(readCount+1)\nnonHS", y="HS\nlog10(readCount+1)")+geom_abline(intercept =0 , slope = 1, color = "gray") + geom_text_repel(data=dataForPlot_SD , aes(label=gene_biotype))
sort(table(dataForPlot_SD$seqid))
sort(table(dataForPlot_SD$gene_biotype))

#correlation
cor(log10(dataForPlot$HS+1), log10(dataForPlot$nonHS+1), method="pearson")
cor(log10(dataForPlot %>% filter(HS!=0 & nonHS!=0) %>% select(1:2) +1), method="pearson")

## Mucus
dataForPlot <-  readCountTable_all_name_melt_metadata_averageReadCount2 %>% dplyr::select(mucus_HS,mucus_nonHS,name,seqid,geneLength,gene_biotype)
colnames(dataForPlot)[1:2] <- gsub("^.*_","",colnames(dataForPlot)[1:2])
dataForPlot_SD <- dataForPlot[which(abs(log10(dataForPlot$HS+1)-log10(dataForPlot$nonHS+1)) > 4*sd(log10(dataForPlot$HS+1)-log10(dataForPlot$nonHS+1))),]
nrow(dataForPlot_SD)/nrow(dataForPlot)*100

ggplot(data=dataForPlot, aes(x=log10(nonHS+1),y=log10(HS+1))) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(readCount+1)\nnonHS", y="HS\nlog10(readCount+1)")+geom_abline(intercept =0 , slope = 1, color = "gray")+geom_point(data=dataForPlot_SD, aes(x=log10(nonHS+1),y=log10(HS+1)),alpha=1, color="blue") +geom_point(data=dataForPlot %>% filter(nonHS==0 | HS==0), aes(x=log10(nonHS+1),y=log10(HS+1)),alpha=1, color="red") 

geneOnAxisX <- dataForPlot %>% dplyr::filter(nonHS!=0 & HS==0)
nrow(geneOnAxisX)/nrow(dataForPlot)*100
sort(table(geneOnAxisX$gene_biotype))
nrow(geneOnAxisX %>% filter(gene_biotype=="protein_coding"))/nrow(dataForPlot)*100
nrow(geneOnAxisX %>% filter(gene_biotype!="protein_coding"))/nrow(dataForPlot)*100
geneOnAxisY <- dataForPlot %>% dplyr::filter(nonHS==0 & HS!=0)
nrow(geneOnAxisY)/nrow(dataForPlot)*100
nrow(geneOnAxisY %>% filter(gene_biotype=="protein_coding"))/nrow(dataForPlot)*100
nrow(geneOnAxisY %>% filter(gene_biotype!="protein_coding"))/nrow(dataForPlot)*100
geneOnAxisXY <- dataForPlot %>% dplyr::filter(nonHS==0 & HS==0)
nrow(geneOnAxisXY)/nrow(dataForPlot)*100
nrow(geneOnAxisXY %>% filter(gene_biotype=="protein_coding"))/nrow(dataForPlot)*100
nrow(geneOnAxisXY %>% filter(gene_biotype!="protein_coding"))/nrow(dataForPlot)*100


ggplot(data=dataForPlot_SD, aes(x=log10(nonHS+1),y=log10(HS+1))) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(readCount+1)\nnonHS", y="HS\nlog10(readCount+1)")+geom_abline(intercept =0 , slope = 1, color = "gray") + geom_text_repel(data=dataForPlot_SD, aes(label=seqid))
#geneLength, name, 
table(dataForPlot_SD$seqid)
table(dataForPlot_SD$gene_biotype)
#correlation
cor(log10(dataForPlot$HS+1), log10(dataForPlot$nonHS+1), method="pearson")
cor(log10(dataForPlot %>% filter(HS!=0 & nonHS!=0) %>% select(1:2) +1), method="pearson")

## Tissue
dataForPlot <-  readCountTable_all_name_melt_metadata_averageReadCount2 %>% dplyr::select(tissue_HS, tissue_nonHS,name,seqid,geneLength,gene_biotype)
colnames(dataForPlot)[1:2] <- gsub("^.*_","",colnames(dataForPlot)[1:2])
dataForPlot_SD <- dataForPlot[which(abs(log10(dataForPlot$HS+1)-log10(dataForPlot$nonHS+1)) > 4*sd(log10(dataForPlot$HS+1)-log10(dataForPlot$nonHS+1))),]
nrow(dataForPlot_SD)/nrow(dataForPlot)*100

ggplot(data=dataForPlot, aes(x=log10(nonHS+1),y=log10(HS+1))) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(readCount+1)\nnonHS", y="HS\nlog10(readCount+1)")+geom_abline(intercept =0 , slope = 1, color = "gray")+geom_point(data=dataForPlot_SD, aes(x=log10(nonHS+1),y=log10(HS+1)),alpha=1, color="blue") +geom_point(data=dataForPlot %>% filter(nonHS==0 | HS==0), aes(x=log10(nonHS+1),y=log10(HS+1)),alpha=1, color="red") + geom_label_repel(data=dataForPlot_SD %>% filter(nonHS>10^2.5), aes(x=log10(nonHS+1),y=log10(HS+1), label=geneLength),alpha=1, color="blue")

ggplot(data=dataForPlot, aes(x=log10(nonHS+1),y=log10(HS+1))) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(TPM+1)\nnonHS", y="HS\nlog10(TPM+1)")+geom_abline(intercept =0 , slope = 1, color = "gray") + geom_label_repel(data=dataForPlot %>% filter(grepl("Gm23935",name)), aes(x=log10(nonHS+1),y=log10(HS+1), label=geneLength),alpha=1, color="blue")+geom_point(data=dataForPlot %>% filter(grepl("Gm23935",name)), aes(x=log10(nonHS+1),y=log10(HS+1), label=geneLength),alpha=1, color="blue")

geneOnAxisX <- dataForPlot %>% dplyr::filter(nonHS!=0 & HS==0)
nrow(geneOnAxisX)/nrow(dataForPlot)*100
sort(table(geneOnAxisX$gene_biotype))
nrow(geneOnAxisX %>% filter(gene_biotype=="protein_coding"))/nrow(dataForPlot)*100
nrow(geneOnAxisX %>% filter(gene_biotype!="protein_coding"))/nrow(dataForPlot)*100
geneOnAxisY <- dataForPlot %>% dplyr::filter(nonHS==0 & HS!=0)
nrow(geneOnAxisY)/nrow(dataForPlot)*100
nrow(geneOnAxisY %>% filter(gene_biotype=="protein_coding"))/nrow(dataForPlot)*100
nrow(geneOnAxisY %>% filter(gene_biotype!="protein_coding"))/nrow(dataForPlot)*100
geneOnAxisXY <- dataForPlot %>% dplyr::filter(nonHS==0 & HS==0)
nrow(geneOnAxisXY)/nrow(dataForPlot)*100
nrow(geneOnAxisXY %>% filter(gene_biotype=="protein_coding"))/nrow(dataForPlot)*100
nrow(geneOnAxisXY %>% filter(gene_biotype!="protein_coding"))/nrow(dataForPlot)*100


ggplot(data=dataForPlot_SD, aes(x=log10(nonHS+1),y=log10(HS+1))) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(readCount+1)\nnonHS", y="HS\nlog10(readCount+1)")+geom_abline(intercept =0 , slope = 1, color = "gray") + geom_text_repel(data=dataForPlot_SD, aes(label=seqid))
#geneLength, name, 
table(dataForPlot_SD$seqid)
table(dataForPlot_SD$gene_biotype)
#correlation
cor(log10(dataForPlot$HS+1), log10(dataForPlot$nonHS+1), method="pearson")
cor(log10(dataForPlot %>% filter(HS!=0 & nonHS!=0) %>% select(1:2) +1), method="pearson")

cs <- colSums(readCountTable_all) #Total mapped reads per sample
rpkm <- t(t(readCountTable_all/readCountTable_name_length$length*(10^9))/cs)
rpkm1 <- rpkm %>% as_tibble() %>% filter(complete.cases(rpkm))
colSums(rpkm1)
data.frame(colSum=colSums(rpkm1)) %>% add_rownames("sampleID") %>% ggplot(aes(x=sampleID, y=colSum)) + geom_label_repel(aes(label=sampleID))

RPKM_Table_name <- bind_cols(readCountTable_name_length %>% filter(complete.cases(rpkm)), rpkm1)
ggpairs((log10(RPKM_Table_name[,-c(1,2)]+1)) %>% sample_n(500, replace=F), 
        upper = list(continuous = my_fn),
        lower = list(continuous = "smooth"))


RPKM_Table_name_melt <- gather(RPKM_Table_name, "sampleID", "RPKM", -name, -length)

RPKM_Table_name_melt_metadata <- left_join(RPKM_Table_name_melt, metadata, by=c("sampleID"="sampleID"))

RPKM_Table_name_melt_metadata_averageTPM <- RPKM_Table_name_melt_metadata %>% dplyr::select(-length,-mouseID, -sampleID) %>% group_by(name, site, enrichment) %>% dplyr::summarise(averageRPKM = mean(RPKM)) %>% ungroup() %>% mutate(siteEnrichment=paste0(site,"_",enrichment)) %>% dplyr::select(-site, -enrichment) %>% spread(siteEnrichment, averageRPKM)

# ggpairs((log10(RPKM_Table_name_melt_metadata_averageTPM[,-1]+1) %>% sample_n(1000, replace=F)) , 
#         upper = list(continuous = my_fn),
#         lower = list(continuous = "smooth")) #

x <- readCountTable_all/readCountTable_name_length$length
naIndex <- which(!complete.cases(x))
x1 <- x[complete.cases(x),]

y <- t(t(x1)*1e6/colSums(x1))
colSums(y)

TMP_Table_name <- bind_cols(readCountTable_name_length %>% filter(complete.cases(x)), y %>% as_tibble())
ggpairs((log10(TMP_Table_name[,-c(1,2)]+1)) %>% sample_n(1000, replace=F), 
        upper = list(continuous = my_fn),
        lower = list(continuous = "smooth"))

TMP_Table_name_melt <- gather(TMP_Table_name, "sampleID", "TPM", -name, -length)
TMP_Table_name_melt_metadata <- left_join(TMP_Table_name_melt, metadata, by=c("sampleID"="sampleID"))
TMP_Table_name_melt_metadata_averageTPM <- TMP_Table_name_melt_metadata %>% dplyr::select(-length,-mouseID, -sampleID) %>% group_by(name, site, enrichment) %>% dplyr::summarise(averageTPM = mean(TPM)) %>% ungroup() %>% mutate(siteEnrichment=paste0(site,"_",enrichment)) %>% dplyr::select(-site, -enrichment) %>% spread(siteEnrichment, averageTPM)


my_diagonalLine <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "black",alpha=0.5) +
    geom_abline(intercept =0 , slope = 1, color = "blue", ...)
  p
  # # grab data
  # x <- eval_data_col(data, mapping$x)
  # y <- eval_data_col(data, mapping$y)
  # ggplot(aes(x=x,y=y))+geom_point()+geom_abline(intercept =0 , slope = 1)
}
lowerFn <- function(data, mapping, method = "loess", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "blue", alpha=0.5) +
    geom_smooth(method = method, color = "red", ...)
  p
}


# ggpairs((log10((TMP_Table_name_melt_metadata_averageTPM %>% filter(!grepl("^mt-", name)) %>% dplyr::select(-name))+1) ) %>% sample_n(5000, replace=F), 
#         upper = list(continuous = wrap("cor", size = 6)),
#         lower = list(continuous = my_diagonalLine )) #%>% dplyr::select(contains("lumen")


## Lumen
ggplot(data=log10((TMP_Table_name_melt_metadata_averageTPM %>% filter(!grepl("^mt-", name)) %>% dplyr::select(-name))+1) %>% dplyr::select(contains("lumen")), aes(x=lumen_nonHS,y=lumen_HS)) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(TPM+1)\nnonHS", y="HS\nlog10(TPM+1)")+geom_abline(intercept =0 , slope = 1, color = "gray")

cor(log10((TMP_Table_name_melt_metadata_averageTPM %>% filter(!grepl("^mt-", name)) %>% dplyr::select(-name))+1) %>% dplyr::select(contains("lumen")))

## Mucus
ggplot(data=log10((TMP_Table_name_melt_metadata_averageTPM %>% filter(!grepl("^mt-", name)) %>% dplyr::select(-name))+1) %>% dplyr::select(contains("mucus")), aes(x=mucus_nonHS,y=mucus_HS)) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(TPM+1)\nnonHS", y="HS\nlog10(TPM+1)")+geom_abline(intercept =0 , slope = 1, color = "gray")
cor(log10((TMP_Table_name_melt_metadata_averageTPM %>% filter(!grepl("^mt-", name)) %>% dplyr::select(-name))+1) %>% dplyr::select(contains("mucus")))

## Tissue
ggplot(data=log10((TMP_Table_name_melt_metadata_averageTPM %>% filter(!grepl("^mt-", name)) %>% dplyr::select(-name))+1) %>% dplyr::select(contains("tissue")), aes(x=tissue_nonHS,y=tissue_HS)) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(TPM+1)\nnonHS", y="HS\nlog10(TPM+1)")+geom_abline(intercept =0 , slope = 1, color = "gray")
cor(log10((TMP_Table_name_melt_metadata_averageTPM %>% filter(!grepl("^mt-", name)) %>% dplyr::select(-name))+1) %>% dplyr::select(contains("tissue")))


### TPM
tpm = function (counts, effective_lengths) {
  rate = log(counts) - log(effective_lengths)
  exp(rate - log(sum(exp(rate))) + log(1E6))
}
readCountTable_name_length_noMT <- readCountTable_name_length[!grepl("^mt-",readCountTable_name_length$name),]
readCountTable_all_noMT <- readCountTable_all[!grepl("^mt-",readCountTable_name_length$name),]

allTPM <- list()
for(i in 1:ncol(readCountTable_all_noMT)){
  allTPM[[i]] <- tpm(unlist(readCountTable_all_noMT[,i]), unlist(readCountTable_name_length_noMT$length))
}
allTPM <- as.data.frame(allTPM) %>% as_tibble()
colnames(allTPM) <- colnames(readCountTable_all_noMT)

ggpairs((log10(allTPM %>% select(ends_with("_HS"))+1) %>% sample_n(500, replace=F)) , 
        #upper = list(continuous = my_fn),
        upper=list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = "smooth")) #
ggpairs((log10(allTPM %>% select(-ends_with("_HS"))+1) %>% sample_n(500, replace=F)) , 
        #upper = list(continuous = my_fn),
        upper=list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = "smooth")) #



colnames(readCountTable_all)
i=11
colnames(readCountTable_all)[c(i,i+1)]
TPM_nonHS <- tpm(unlist(readCountTable_all[,i]), unlist(readCountTable_name_length$length)[])
readCount_nonHS <- unlist(readCountTable_all[,i])

j=i+1
j=17
TPM_HS <- tpm(unlist(readCountTable_all[,j]), unlist(readCountTable_name_length$length)[])
readCount_HS <- unlist(readCountTable_all[,j])

dataForPlot<-data.frame(nonHS=TPM_nonHS, HS=TPM_HS)
ggplot(data=dataForPlot, aes(x=log10(nonHS+1), y=log10(HS+1))) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(TPM+1)\nnonHS", y="HS\nlog10(TPM+1)")+geom_abline(intercept =0 , slope = 1, color = "gray")
cor(log10(dataForPlot$nonHS+1), log10(dataForPlot$HS+1), method="spearman")
cor(log10(dataForPlot$nonHS+1), log10(dataForPlot$HS+1), method="pearson")

cor(TPM_nonHS, TPM_HS, method="pearson")

geneOnAxisX <- dataForPlot %>% dplyr::filter(nonHS!=0 & HS==0)
nrow(geneOnAxisX)/nrow(dataForPlot)*100
sort(table(geneOnAxisX$gene_biotype))
geneOnAxisY <- dataForPlot %>% dplyr::filter(nonHS==0 & HS!=0)
nrow(geneOnAxisY)/nrow(dataForPlot)*100

ggplot(data=data.frame(nonHS=readCount_nonHS, HS=readCount_HS), aes(x=log10(nonHS+1), y=log10(HS+1))) + geom_point(alpha=0.5) + theme_bw() + labs(x="log10(readCount+1)\nnonHS", y="HS\nlog10(readCount+1)")+geom_abline(intercept =0 , slope = 1, color = "gray")
cor(readCount_nonHS, readCount_HS, method="spearman")
cor(readCount_nonHS, readCount_HS, method="pearson")

# average TPM of HS vs nonHS excluding genes with no more than 10 reads across all samples
averageTPM_removal10_ls <- list()
TPM_removal10_ls <- list()
readCount_df <- as.data.frame(readCountTable_all_name)
row.names(readCount_df) <- readCount_df$name

j=1
for(site in c("lumen","mucus","tissue")){
  readCount_df_oneSite <- readCount_df[,grep(site,colnames(readCount_df))]
  readCount_df_oneSite_removal10 <- readCount_df_oneSite[((apply(readCount_df_oneSite[,1:3],1, min) > 10)*1 + (apply(readCount_df_oneSite[,4:6],1, min) > 10)*1)>0,]; dim(readCount_df_oneSite_removal10)
  readCount_df_oneSite_removal10_geneLength <- (merge(readCount_df_oneSite_removal10,mm10_gene[,c("gene_name","geneLength","gene_biotype")], by.x="row.names",by.y="gene_name"))
  
  TPM_oneSite_removal10 <- list()
  for(i in 2:7){
    readsPerKilobase_oneSample <- readCount_df_oneSite_removal10_geneLength[,i]/(readCount_df_oneSite_removal10_geneLength[,"geneLength"]/1000)
    perMillionScalingFactor_oneSample <- sum(readsPerKilobase_oneSample)/1e6
    TPM_oneSample <- readsPerKilobase_oneSample/perMillionScalingFactor_oneSample
    TPM_oneSite_removal10[[i-1]] <- TPM_oneSample
  }
  TPM_oneSite_removal10 <- do.call("cbind",TPM_oneSite_removal10)
  colnames(TPM_oneSite_removal10) <- paste0(rep(c("HS","nonHS"),each=3), c(1:3))
  TPM_removal10_ls[[j]] <- data.frame(TPM_oneSite_removal10,locusTag=readCount_df_oneSite_removal10_geneLength$Row.names, biotype=readCount_df_oneSite_removal10_geneLength$gene_biotype)
  averageTPM_oneSite_removal10 <- data.frame(HS=rowMeans(TPM_oneSite_removal10[,1:3]), nonHS=rowMeans(TPM_oneSite_removal10[,4:6]), geneLength=readCount_df_oneSite_removal10_geneLength$geneLength, locusTag=readCount_df_oneSite_removal10_geneLength$Row.names, biotype=readCount_df_oneSite_removal10_geneLength$gene_biotype)
  averageTPM_removal10_ls[[j]] <- averageTPM_oneSite_removal10
  j=j+1
}
names(TPM_removal10_ls) <- names(averageTPM_removal10_ls) <- c("lumen","mucus","tissue")



# FigS3 
ggplots_ls <- list()
for(i in 1:length(averageTPM_removal10_ls)){
  averageTPM_oneSite_removal10 <- averageTPM_removal10_ls[[i]]
  site <- names(averageTPM_removal10_ls)[i]
  #p<-ggplot(data=averageTPM_oneSite_removal10, aes(y=log10(HS+1),x=log10(nonHS+1), label=locusTag))+theme_classic()+coord_fixed()+labs(title=paste0(site,": ",nrow(averageTPM_oneSite_removal10)," genes remaining"))+geom_text(angle = 45, size =3, alpha=0.5, color="blue")+geom_point(aes(size=geneLength),alpha=1,shape=1)+xlim(0,6)+ylim(0,6)
  p<-ggplot(data=averageTPM_oneSite_removal10, aes(y=log10(HS+1),x=log10(nonHS+1), label=locusTag))+theme_classic()+coord_fixed()+labs(title=paste0(site,"; cor=",round(cor(log10(averageTPM_oneSite_removal10[,1:2]+1))[1,2],2)), x="log10(TPM+1)\nnonHS", y="HS\nlog10(TPM+1)\n")+geom_point(alpha=0.4,size=1,shape=1)+scale_x_continuous(breaks=seq(0,6,1), limits=c(0,6)) + scale_y_continuous(breaks=seq(0,6,1), limits=c(0,6))+geom_abline(intercept = 0, slope = 1, linetype="dashed",color="grey")
  print(p)
  ggplots_ls[[i]] <- p
}
lapply(averageTPM_removal10_ls,function(x) nrow(x))
figure <- ggarrange(ggplots_ls[[1]],ggplots_ls[[2]],ggplots_ls[[3]],
                    ncol =3, nrow =1)
print(figure)



# FigS_HS_vs_nonHS_individual_mice for host RNA
ggplots_ls <- list()
k=1
for(i in 1:length(TPM_removal10_ls)){
  TPM_oneSite_removal10 <- TPM_removal10_ls[[i]]
  site <- names(TPM_removal10_ls)[i]
  for(j in 1:3){
    TPM_oneSite_removal10_oneMouse <- TPM_oneSite_removal10[,c(j,j+3,7)]
    colnames(TPM_oneSite_removal10_oneMouse)[1:2] <- c("HS","nonHS")
    p<-ggplot(data=TPM_oneSite_removal10_oneMouse, aes(y=log10(HS+1),x=log10(nonHS+1), label=locusTag))+theme_classic()+coord_fixed()+labs(title=paste0(site,"-mouse",j,"; cor=",round(cor(log10(TPM_oneSite_removal10_oneMouse[,1:2]+1))[1,2],2)), x="log10(TPM+1)\nnonHS", y="HS\nlog10(TPM+1)\n")+geom_point(alpha=0.4,size=1,shape=1)+scale_x_continuous(breaks=seq(0,6,1), limits=c(0,6)) + scale_y_continuous(breaks=seq(0,6,1), limits=c(0,6))+geom_abline(intercept = 0, slope = 1, linetype="dashed",color="grey")
    print(p)
    ggplots_ls[[k]] <- p
    k <- k+1
  }
}
figure <- do.call("ggarrange", c(ggplots_ls, nrow=3,ncol=3))
print(figure)
