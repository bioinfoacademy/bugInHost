library(tidyverse)

# get the raw read count table
readCount_df <- read.csv("~/Google Drive/bugInHost.Mazmanian/docs/manuscript/revisions_bug_in_host/files_from_wenchi/FullList_RawReadCounts_EachMouseSite.csv", row.names = 1)
readCount_df
head(readCount_df)

# get the gene length of each locustag
latestLocusTag<-read.csv("~/Google Drive/bugInHost.Mazmanian/data/reference/GCF_000025985.1_ASM2598v1_genomic.gff.locusTag", header=F, sep="\t")
head(latestLocusTag)
getLengthLocusTag <- data.frame(do.call("rbind",strsplit(as.character(latestLocusTag$V1),"-")))
getLengthLocusTag$locusTag <- latestLocusTag$V2
getLengthLocusTag$geneLength <- as.numeric(as.character(getLengthLocusTag$X3)) - as.numeric(as.character(getLengthLocusTag$X2)) + 1
head(getLengthLocusTag)

# average TPM of HS vs nonHS excluding genes with no more than 10 reads across all samples
averageTPM_removal10_ls <- list()
TPM_removal10_ls <- list()
j=1
for(site in c("lumen","mucus","tissue")){
  readCount_df_oneSite <- readCount_df[,grep(site,colnames(readCount_df))]
  #readCount_df_oneSite_removal10 <- readCount_df_oneSite[(apply(readCount_df_oneSite,1, min) > 10),]
  readCount_df_oneSite_removal10 <- readCount_df_oneSite[((apply(readCount_df_oneSite[,1:3],1, min) > 10)*1 + (apply(readCount_df_oneSite[,4:6],1, min) > 10)*1)>0,]; dim(readCount_df_oneSite_removal10)
  readCount_df_oneSite_removal10_geneLength <- (merge(readCount_df_oneSite_removal10,getLengthLocusTag[,c("locusTag","geneLength")], by.x="row.names",by.y="locusTag"))
  
  TPM_oneSite_removal10 <- list()
  for(i in 2:7){
    readsPerKilobase_oneSample <- readCount_df_oneSite_removal10_geneLength[,i]/(readCount_df_oneSite_removal10_geneLength[,"geneLength"]/1000)
    perMillionScalingFactor_oneSample <- sum(readsPerKilobase_oneSample)/1e6
    TPM_oneSample <- readsPerKilobase_oneSample/perMillionScalingFactor_oneSample
    #   x2_new <- x_newlengthNormalized/perMillionScalingFactor
    TPM_oneSite_removal10[[i-1]] <- TPM_oneSample
  }
  TPM_oneSite_removal10 <- do.call("cbind",TPM_oneSite_removal10)
  colnames(TPM_oneSite_removal10) <- paste0(rep(c("HS","nonHS"),each=3), c(1:3))
  TPM_removal10_ls[[j]] <- data.frame(TPM_oneSite_removal10,locusTag=readCount_df_oneSite_removal10_geneLength$Row.names)
  averageTPM_oneSite_removal10 <- data.frame(HS=rowMeans(TPM_oneSite_removal10[,1:3]), nonHS=rowMeans(TPM_oneSite_removal10[,4:6]), geneLength=readCount_df_oneSite_removal10_geneLength$geneLength, locusTag=readCount_df_oneSite_removal10_geneLength$Row.names)
  averageTPM_removal10_ls[[j]] <- averageTPM_oneSite_removal10
  j=j+1
}
names(TPM_removal10_ls) <- names(averageTPM_removal10_ls) <- c("lumen","mucus","tissue")

ggplots_ls <- list()
k=1
for(i in 1:length(TPM_removal10_ls)){
  TPM_oneSite_removal10 <- TPM_removal10_ls[[i]]
  site <- names(TPM_removal10_ls)[i]
  for(j in 1:3){
    TPM_oneSite_removal10_individual <- TPM_oneSite_removal10[,c(j,j+3)]
    colnames(TPM_oneSite_removal10_individual) <- c("HS", "nonHS")
    p<-ggplot(data=TPM_oneSite_removal10_individual, aes(y=log10(HS+1),x=log10(nonHS+1)))+theme_classic()+coord_fixed()+labs(title=paste0(site,j,"; cor=",round(cor(TPM_oneSite_removal10_individual[,1:2])[1,2],2)), x="log10(TPM)\nnonHS", y="HS\nlog10(TPM)\n")+geom_point(alpha=0.4,size=1,shape=1)+scale_x_continuous(breaks=seq(0,6,1), limits=c(0,6)) + scale_y_continuous(breaks=seq(0,6,1), limits=c(0,6))+geom_abline(intercept = 0, slope = 1, linetype="dashed",color="grey")
    #print(p)
    ggplots_ls[[k]] <- p
    k <- k+1
  }
}
length(ggplots_ls)
# Plot figure S2  ---------------------------------------------------------------
print(ggplots_ls)
