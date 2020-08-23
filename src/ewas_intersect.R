project_path <- "/cluster/home/chenhy/project/BPD/"
gene_ttest_paired <- read.csv(sprintf("%s/result/gene_ttest_paired.csv",project_path))
gene_ttest_bootstrap <- read.csv(sprintf("%s/result/gene_ttest.csv",project_path))
inter_Dis_gene <- intersect(gene_ttest_paired$Dis_gene,gene_ttest_bootstrap$Dis_gene)
inter_gene_ttest_bootstrap <- gene_ttest_bootstrap[which(gene_ttest_bootstrap$Dis_gene %in% inter_Dis_gene),]
inter_gene_ttest_paired <- gene_ttest_paired[which(gene_ttest_paired$Dis_gene %in% inter_Dis_gene),]
write.csv(inter_gene_ttest_bootstrap,file = sprintf("%s/result/inter_gene_ttest_bootstrap.csv",project_path),row.names = F)
write.csv(inter_gene_ttest_paired,file = sprintf("%s/result/inter_gene_ttest_paired.csv",project_path),row.names = F)

variant_paired <- read.csv(sprintf("%s/result/variant_fisher_paired.csv",project_path))
variant_bootstrap <- read.csv(sprintf("%s/result/variant_fisher.csv",project_path))
inter_Dis_variant <- intersect(variant_paired$Dis_variant,variant_bootstrap$Dis_variant)
inter_variant_bootstrap <- variant_bootstrap[which(variant_bootstrap$Dis_variant %in% inter_Dis_variant),]
inter_variant_paired <- variant_paired[which(variant_paired$Dis_variant %in% inter_Dis_variant),]
write.csv(inter_variant_bootstrap,file = sprintf("%s/result/inter_variant_bootstrap.csv",project_path),row.names = F)
write.csv(inter_variant_paired,file = sprintf("%s/result/inter_variant_paired.csv",project_path),row.names = F)
