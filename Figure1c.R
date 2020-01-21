library(tidyverse)

coverageFiles <- grep("BF.*psorted.coverage$", list.files("~/PATH/", full.names=T),value=T)
#length(coverageFiles)
coverageListALL=NULL
for(i in 1:length(coverageFiles)){
  oneCoverage <- read_table2(coverageFiles[i], col_names = F)
  coverageListALL[[i]] <- oneCoverage$X7
}
length(coverageListALL)
names(coverageListALL) <- sub("\\.bug.*$","",sub("^.*BF","BF",coverageFiles))
coverageListALL_metadata <- tibble(group = names(coverageListALL))
coverageListALL_metadata$enrichment <- ifelse(grepl("HS", coverageListALL_metadata$group), "HS", "nonHS")
coverageListALL_metadata$site <- sub("BF[1-3]","",sub("_.*$","",coverageListALL_metadata$group))
coverageListALL_metadata$SN <- c(1:nrow(coverageListALL_metadata))

# lumen nonHS
selectedGroups <- coverageListALL_metadata %>% filter(site=="lumen" & enrichment=="nonHS") %>% select(group) %>% unlist()
coverage_lumen_nonHS <- coverageListALL[selectedGroups] %>% do.call(cbind, .) %>% rowMeans()
# lumen HS
selectedGroups <- coverageListALL_metadata %>% filter(site=="lumen" & enrichment=="HS") %>% select(group) %>% unlist()
coverage_lumen_HS <- coverageListALL[selectedGroups] %>% do.call(cbind, .) %>% rowMeans()

# mucus nonHS
selectedGroups <- coverageListALL_metadata %>% filter(site=="mucus" & enrichment=="nonHS") %>% select(group) %>% unlist()
coverage_mucus_nonHS <- coverageListALL[selectedGroups] %>% do.call(cbind, .) %>% rowMeans()
# mucus HS
selectedGroups <- coverageListALL_metadata %>% filter(site=="mucus" & enrichment=="HS") %>% select(group) %>% unlist()
coverage_mucus_HS <- coverageListALL[selectedGroups] %>% do.call(cbind, .) %>% rowMeans()

# tissue nonHS
selectedGroups <- coverageListALL_metadata %>% filter(site=="tissue" & enrichment=="nonHS") %>% select(group) %>% unlist()
coverage_tissue_nonHS <- coverageListALL[selectedGroups] %>% do.call(cbind, .) %>% rowMeans()
# tissue HS
selectedGroups <- coverageListALL_metadata %>% filter(site=="tissue" & enrichment=="HS") %>% select(group) %>% unlist()
coverage_tissue_HS <- coverageListALL[selectedGroups] %>% do.call(cbind, .) %>% rowMeans()

mean_coverage_all <- tibble(lumen_nonHS=coverage_lumen_nonHS, lumen_HS=coverage_lumen_HS, mucus_nonHS=coverage_mucus_nonHS, mucus_HS=coverage_mucus_HS, tissue_nonHS=coverage_tissue_nonHS, tissue_HS=coverage_tissue_HS) %>% gather("group","coverage")

mean_coverage_all$enrichment <- ifelse(grepl("nonHS", mean_coverage_all$group), "nonHS", "HS")
mean_coverage_all$enrichment <- factor(mean_coverage_all$enrichment, levels=c("nonHS","HS"))
mean_coverage_all$site <- sub("BF[1-3]","",sub("_.*$","",mean_coverage_all$group))

# Plot figure 1c ---------------------------------------------------------------
ggplot(data=mean_coverage_all, aes(x=site,y=coverage, fill=enrichment)) + geom_boxplot(outlier.shape = NA)+scale_fill_brewer(palette="Paired") + theme_bw()+labs(fill="", y="Average Read Coverage")
