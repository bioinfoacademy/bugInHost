library(tidyverse)

# get the raw read count table
readCount_df <- read.csv("~/PATH/FullList_RawReadCounts_EachMouseSite.csv", row.names = 1)
readCount_df
head(readCount_df)

# get the gene length of each locustag
latestLocusTag<-read.csv("~/PATH/GCF_000025985.1_ASM2598v1_genomic.gff.locusTag", header=F, sep="\t")
head(latestLocusTag)
getLengthLocusTag <- data.frame(do.call("rbind",strsplit(as.character(latestLocusTag$V1),"-")))
getLengthLocusTag$locusTag <- latestLocusTag$V2
getLengthLocusTag$geneLength <- as.numeric(as.character(getLengthLocusTag$X3)) - as.numeric(as.character(getLengthLocusTag$X2)) + 1
head(getLengthLocusTag)

# average TPM of HS vs nonHS excluding genes with no more than 10 reads across all samples
averageTPM_ls <- list()
TPM_ls <- list()
j=1
for(site in c("lumen","mucus","tissue")){
  readCount_df_oneSite <- readCount_df[,grep(site,colnames(readCount_df))]
  readCount_df_oneSite_geneLength <- (merge(readCount_df_oneSite, getLengthLocusTag[,c("locusTag","geneLength")], by.x="row.names",by.y="locusTag"))
  
  TPM_oneSite <- list()
  for(i in 2:7){
    readsPerKilobase_oneSample <- readCount_df_oneSite_geneLength[,i]/(readCount_df_oneSite_geneLength[,"geneLength"]/1000)
    perMillionScalingFactor_oneSample <- sum(readsPerKilobase_oneSample)/1e6
    TPM_oneSample <- readsPerKilobase_oneSample/perMillionScalingFactor_oneSample
    TPM_oneSite[[i-1]] <- TPM_oneSample
  }
  TPM_oneSite <- do.call("cbind",TPM_oneSite)
  colnames(TPM_oneSite) <- paste0(rep(c("HS","nonHS"),each=3), c(1:3))
  TPM_ls[[j]] <- data.frame(TPM_oneSite,locusTag=readCount_df_oneSite_geneLength$Row.names)
  averageTPM_oneSite <- data.frame(HS=rowMeans(TPM_oneSite[,1:3]), nonHS=rowMeans(TPM_oneSite[,4:6]), geneLength=readCount_df_oneSite_geneLength$geneLength, locusTag=readCount_df_oneSite_geneLength$Row.names)
  averageTPM_ls[[j]] <- averageTPM_oneSite
  j=j+1
}
names(TPM_ls) <- names(averageTPM_ls) <- c("lumen","mucus","tissue")

# Plot figure S5 nonHS lumen vs mucus ---------------------------------------------------------------
averageTPM_twoSites <- merge(averageTPM_ls[["lumen"]] %>% select(nonHS, locusTag), averageTPM_ls[["mucus"]] %>% select(nonHS, locusTag), by="locusTag")
colnames(averageTPM_twoSites)[2:3] <- c("lumen","mucus")
ggplot(data=averageTPM_twoSites, aes(x=log10(lumen+1),y=log10(mucus+1), label=locusTag))+theme_classic()+coord_fixed()+labs(title=paste0("cor=",round(cor(averageTPM_twoSites[,2:3], method="spearman")[1,2],2)), x="log10(TPM+1)\nnonHS lumen", y="nonHS mucus\nlog10(TPM+1)\n")+geom_point(alpha=0.4,size=1,shape=1)+scale_x_continuous(breaks=seq(0,6,1), limits=c(0,6)) + scale_y_continuous(breaks=seq(0,6,1), limits=c(0,6))+geom_abline(intercept = 0, slope = 1, linetype="dashed",color="grey")

# Plot figure S5 nonHS lumen vs tissue ---------------------------------------------------------------
averageTPM_twoSites <- merge(averageTPM_ls[["lumen"]] %>% select(nonHS, locusTag), averageTPM_ls[["tissue"]] %>% select(nonHS, locusTag), by="locusTag")
colnames(averageTPM_twoSites)[2:3] <- c("lumen","tissue")
ggplot(data=averageTPM_twoSites, aes(x=log10(lumen+1),y=log10(tissue+1), label=locusTag))+theme_classic()+coord_fixed()+labs(title=paste0("cor=",round(cor(averageTPM_twoSites[,2:3], method="spearman")[1,2],2)), x="log10(TPM+1)\nnonHS lumen", y="nonHS tissue\nlog10(TPM+1)\n")+geom_point(alpha=0.4,size=1,shape=1)+scale_x_continuous(breaks=seq(0,6,1), limits=c(0,6)) + scale_y_continuous(breaks=seq(0,6,1), limits=c(0,6))+geom_abline(intercept = 0, slope = 1, linetype="dashed",color="grey")

# Plot figure S5 HS lumen vs mucus ---------------------------------------------------------------
averageTPM_twoSites <- merge(averageTPM_ls[["lumen"]] %>% select(HS, locusTag), averageTPM_ls[["mucus"]] %>% select(HS, locusTag), by="locusTag")
colnames(averageTPM_twoSites)[2:3] <- c("lumen","mucus")
ggplot(data=averageTPM_twoSites, aes(x=log10(lumen+1),y=log10(mucus+1), label=locusTag))+theme_classic()+coord_fixed()+labs(title=paste0("cor=",round(cor(averageTPM_twoSites[,2:3], method="spearman")[1,2],2)), x="log10(TPM+1)\nHS lumen", y="HS mucus\nlog10(TPM+1)\n")+geom_point(alpha=0.4,size=1,shape=1)+scale_x_continuous(breaks=seq(0,6,1), limits=c(0,6)) + scale_y_continuous(breaks=seq(0,6,1), limits=c(0,6))+geom_abline(intercept = 0, slope = 1, linetype="dashed",color="grey")

# Plot figure S5 HS lumen vs tissue ---------------------------------------------------------------
averageTPM_twoSites <- merge(averageTPM_ls[["lumen"]] %>% select(HS, locusTag), averageTPM_ls[["tissue"]] %>% select(HS, locusTag), by="locusTag")
colnames(averageTPM_twoSites)[2:3] <- c("lumen","tissue")
ggplot(data=averageTPM_twoSites, aes(x=log10(lumen+1),y=log10(tissue+1), label=locusTag))+theme_classic()+coord_fixed()+labs(title=paste0("cor=",round(cor(averageTPM_twoSites[,2:3], method="spearman")[1,2],2)), x="log10(TPM+1)\nHS lumen", y="HS tissue\nlog10(TPM+1)\n")+geom_point(alpha=0.4,size=1,shape=1)+scale_x_continuous(breaks=seq(0,6,1), limits=c(0,6)) + scale_y_continuous(breaks=seq(0,6,1), limits=c(0,6))+geom_abline(intercept = 0, slope = 1, linetype="dashed",color="grey")
