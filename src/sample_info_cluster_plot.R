library(tidyverse)
library(stringr)
library(DMwR)
load("/cluster/home/chenhy/project/BPD/data/bpd_info.RData")
bpd_df <- bpd_df[-manyNAs(bpd_df[,1:ncol(bpd_df)],15/(ncol(bpd_df)-1)),]
rownames(bpd_df) <- bpd_df$id;bpd_df <- select(bpd_df,-id)
bpd_df <- t(bpd_df)
