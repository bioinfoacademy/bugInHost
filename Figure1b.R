library(tidyverse)

mappingResult <- read_csv("~/PATH/overallMappingBugHostUnaligned.csv", col_names=F)
head(mappingResult)
colnames(mappingResult)=c("sample","alignedToBug","alignedToHost","unaligned")

mappingResult_percentage <- (mappingResult %>% select(-sample))/rowSums(mappingResult %>% select(-sample))
mappingResult_percentage <- bind_cols(mappingResult %>% select(sample), mappingResult_percentage) %>% dplyr::filter(!grepl("^ccf", sample))
mappingResult_percentage$group <- gsub("^BF[1-3]","",mappingResult_percentage$sample)

mappingResult_percentage_summary <- mappingResult_percentage %>% select(alignedToBug, group) %>% group_by(group) %>% summarise(mean = mean(alignedToBug), sd = sd(alignedToBug))
mappingResult_percentage_summary$enrichment <- ifelse(grepl("HS", mappingResult_percentage_summary$group), "HS", "nonHS")
mappingResult_percentage_summary$enrichment <- factor(mappingResult_percentage_summary$enrichment, levels=c("nonHS","HS"))
mappingResult_percentage_summary$site <- gsub("_.*$","",mappingResult_percentage_summary$group)

# Plot figure 1b ---------------------------------------------------------------
ggplot(mappingResult_percentage_summary, aes(x=site, y=mean, fill=enrichment)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + theme(legend.position="top")+
  scale_fill_brewer(palette="Paired")+ylab("% bacterial RNA") + theme_bw()+labs(fill="")
