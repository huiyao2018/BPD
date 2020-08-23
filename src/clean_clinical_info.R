library(tidyverse)
library(stringr)
## read data
colname_df <- openxlsx::read.xlsx("/cluster/home/chenhy/project/BPD/data/BPD_临床信息汇总20200723_tochy.xlsx",sheet=4,colNames=F)
bpd_info <- openxlsx::read.xlsx("/cluster/home/chenhy/project/BPD/data/BPD_临床信息汇总20200723_tochy.xlsx",sheet=1) %>% .[which(is.na(.$`样本编号2`)==F),]
colnames(bpd_info) <- colname_df$X2;#bpd_info$diagnosis <- 1
rownames(bpd_info) <- str_c("case_",seq(1:nrow(bpd_info)))
no_bpd_info <- openxlsx::read.xlsx("/cluster/home/chenhy/project/BPD/data/BPD_临床信息汇总20200723_tochy.xlsx",sheet=2) %>% .[which(is.na(.$`样本编号2`)==F),]
colnames(no_bpd_info) <- colname_df$X4;#no_bpd_info$diagnosis <- 0
rownames(no_bpd_info) <- str_c("case_",seq(1:nrow(no_bpd_info)))
#summary sample info
sample_info <- rbind(select(bpd_info,Household_number,Sample_number2),select(no_bpd_info,Household_number,Sample_number2))
sample_info %>% dplyr::group_by(Household_number) %>% summarise(n=n()) %>% filter(n>1) %>% dplyr::count(n)
#output bpd_info
candidate_var <- c(openxlsx::read.xlsx("/cluster/home/chenhy/project/BPD/data/BPD_临床信息汇总20200723_tochy.xlsx",sheet=5,colNames=F)[,1])
bpd_df <- rbind(select(bpd_info,all_of(candidate_var)),select(no_bpd_info,all_of(candidate_var)))
for (i in colnames(bpd_df)) {
  bpd_df[,i] <- str_trim(bpd_df[,i])
}
bpd_df$Nervous_system_infections <- NA;bpd_df$urinary_tract_infections <- NA;bpd_df$digestive_tract_infections <- NA
for (i in 1:nrow(bpd_df)) {
  bpd_df[i,"Nervous_system_infections"] <- ifelse(str_detect(bpd_df[i,"Other_systemic_infections"],"1"),1,0)
  bpd_df[i,"urinary_tract_infections"] <- ifelse(str_detect(bpd_df[i,"Other_systemic_infections"],"2"),1,0)
  bpd_df[i,"digestive_tract_infections"] <- ifelse(str_detect(bpd_df[i,"Other_systemic_infections"],"3"),1,0)
}
bpd_df <- select(bpd_df,-Other_systemic_infections)
bpd_df$gestational_age_week <- as.numeric(str_split(bpd_df$gestational_age,"\\+",simplify = T)[,1])*7
bpd_df$gestational_age_day <- as.numeric(str_split(bpd_df$gestational_age,"\\+",simplify = T)[,2])
bpd_df$gestational_age <- apply(select(bpd_df,gestational_age_week,gestational_age_day),1,function(x){sum(x,na.rm = T)})
bpd_df <- select(bpd_df,-gestational_age_week,-gestational_age_day)
numeric_var <- c("gestational_age",colnames(bpd_df)[str_detect(colnames(bpd_df),"Number*")])
factor_var <- colnames(bpd_df)[-which(colnames(bpd_df) %in% numeric_var)]
factor_order_var <- c("BPD_scale","Ventilation_at_36_weeks")
id <- rownames(bpd_df)
bpd_df$genders <- ifelse(bpd_df$genders=="男",1,0)
bpd_df <- bpd_df %>% dplyr::mutate_at(colnames(bpd_df),~as.numeric(.x))
rownames(bpd_df) <- id
# bpd_df <- bpd_df %>% dplyr::mutate_at(numeric_var,~as.numeric(.x))
# bpd_df <- bpd_df %>% dplyr::mutate_at(factor_var,~as.factor(.x))
# bpd_df <- bpd_df %>% dplyr::mutate_at(factor_order_var,~factor(.x,ordered = T))
str(bpd_df)
save(bpd_df,file="/cluster/home/chenhy/project/BPD/data/bpd_info.RData")
