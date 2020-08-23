################################
library(openxlsx) ## for input into excel
library(plyr) ## for xbind.fill data
library(dplyr) ## for rename colname
library(magrittr) ## for %>% operating
library(stringr) ## for character operating
library(reshape) ## for reshape dataframe
library(car) ## for leveneTest 


##' Geting file path and input data
#'
#' \code{read_file} is a function for geting file path and checking the existence of file path. 
#' 
#'
#' @param sample_info data.frame, is a subset of Blood_SampleList.
#' @return Return a list.
read_file <- function(sample_info) {
  tryCatch({
    existence_sample <- c()
    exist_path <- c()
    no_file <- c()
    # judge if exists file
    for (i in 1:nrow(sample_info)) {
      tmp <- str_c("/cluster/home/luyulan/project/Combined_Pipeline/data/dump_", sample_info[i, "Test_type"], 
                   "/", sample_info[i, "V10"], "/", sample_info[i, "file_name"], "/result_", sample_info[i, "file_name"], 
                   "/", sample_info[i, "Sample_name"], "_Essential_slt_report_AI.xlsx")
      if (file.exists(tmp) == 1) {
        existence_sample <- c(existence_sample, sample_info[i, "uniq_id"])
        exist_path <- c(exist_path, tmp)
      }else{
        no_file <- c(no_file,tmp)
      }
      rm(tmp)
    }
    save(no_file,file = sprintf("%s/no_file_path.RData",project_path))
    exist_path_list <- list()
    for (i in 1:length(exist_path)) {
      exist_path_list[[i]] <- exist_path[i]
    }
    # 读取文件
    sample_list <- list()
    sample_list <- lapply(exist_path_list, function(x) {
      want_var <- c("Dis_gene", "Pathogenicity_Score-Gene")
      tryCatch({
        tmp <- openxlsx::read.xlsx(x, 1)
        tmp2 <- if (all(want_var %in% colnames(tmp))) {
          dplyr::select(tmp, starts_with("Report"), Dis_gene, `Pathogenicity_Score-Gene`)
        }
        return(tmp2)
      }, error = function(e) {
        cat(sprintf("#Cannot open file! -> %s\n", x))
      })
    })
    names(sample_list) <- existence_sample
    return(sample_list)
  }, error = function(e) {
    cat("#Cannot open file!\n")
  })
}


##' Merge all data in sample list
#'
#' \code{merge_file} is a function to merge all dataframe in a list. 
#' 
#'
#' @param sample_list list, produced by function read_file.
#' @return Return a dataframe. 
merge_file <- function(sample_list){
  sample_file <- data.frame()
  for (i in 1:length(sample_list)) {
    if (is.data.frame(sample_list[[i]]) == 1) {
      sample_list[[i]]$uniq_id <- names(sample_list[i])
      sample_file <- plyr::rbind.fill(sample_file,sample_list[[i]])
    }
  }
  #sample_list <- purrr::map(sample_list, ~purrr::compact(.)) %>% purrr::keep(~length(.) != 0)  #Remove Null value in list
  #sample_file <- do.call("rbind.fill", sample_list)
  return(sample_file)
}

##' Calculate p value by compare with case sample many times.
#'
#' \code{cal_p} is a function to make t-test and fisher-test.
#' 
#'
#' @param x dataframe, contains all uniq_id of control.
#' @return Return a list. It's the final results including t-test and fisher-test. 
cal_p <- function(x) {
  #x <- ctr_mat[,1]
  tryCatch({
    ctr_sub_list <- ctr_all_file[x]
    ctr_file <- merge_file(ctr_sub_list)
    ctr_file$diagnosis <- "control"
    sample_file <- plyr::rbind.fill(case_file, ctr_file)
    #sample_file <- sample_file[which(sample_file$Dis_gene %in% use_gene$Dis_gene), ]
    want_var <- c("uniq_id", "Dis_gene", "Pathogenicity_Score-Gene")
    sample_file <- sample_file[!duplicated(sample_file[, want_var]), ]  #remove the duplicated sample_name and Dis_gene
    
    mat_case <- cast(sample_file[which(sample_file$diagnosis == "case"), want_var], Dis_gene ~ uniq_id)
    mat_control <- cast(sample_file[which(sample_file$diagnosis == "control"), want_var], Dis_gene ~ uniq_id)
    
    score_count1 <- ncol(mat_case) - apply(mat_case, 1, function(x) sum(is.na(x)))
    score_count1 <- score_count1[-which(score_count1 < 4)]
    mat_case <- mat_case[which(mat_case$Dis_gene %in% names(score_count1)), ]
    
    score_count2 <- ncol(mat_control) - apply(mat_control, 1, function(x) sum(is.na(x)))
    score_count2 <- score_count2[-which(score_count2 < 4)]
    mat_control <- mat_control[which(mat_control$Dis_gene %in% names(score_count2)), ]
    mat <- merge(mat_case, mat_control, by = "Dis_gene", all = F)
    
    case_name <- sample_file[which(sample_file$diagnosis == "case"), "uniq_id"]
    case_name <- case_name[!duplicated(case_name)]
    ctr_name <- sample_file[which(sample_file$diagnosis == "control"), "uniq_id"]
    ctr_name <- ctr_name[!duplicated(ctr_name)]
    # calculate P value
    rownames(mat) <- mat[,1];mat <- mat[,-1];mat <- as.matrix(mat)
    ttest_p_value <- unlist(apply(mat,1,function(w){
      tryCatch({
        case_tmp <- rep("case",length(w[case_name]))
        ctr_tmp <- rep("ctr",length(w[ctr_name]))
        test.df <- data.frame(w=w,
                              g=c(case_tmp,ctr_tmp))
        res.levene.test <- leveneTest(w~g,data = test.df)
        if(res.levene.test$`Pr(>F)`[1] > 0.05){
          return(t.test(w[case_name], w[ctr_name], alternative="greater",mu=0,var.equal=T)[["p.value"]])
        }else{
          return(t.test(w[case_name], w[ctr_name], alternative="greater",mu=0,var.equal=F)[["p.value"]])
        }
      },error = function(e) {
        cat("#Scores of gene are essentially constant!\n")
      })
    }))
    gene_case_n <- score_count1-1;gene_case_n <- gene_case_n[names(ttest_p_value)]
    gene_case_n_not <- length(case_name)-gene_case_n
    gene_control_n <- score_count2-1;gene_control_n <-gene_control_n[names(ttest_p_value)]
    gene_control_n_not <- length(ctr_name)-gene_control_n
    res <- data.frame(ttest_p_value=ttest_p_value,gene_case_n=gene_case_n,gene_case_n_not=gene_case_n_not,gene_control_n=gene_control_n,gene_control_n_not=gene_control_n_not)
    res$Dis_gene <- rownames(res)
    res <- merge(sample_file, res, by = "Dis_gene")
    
    chi_table <- res %>% dplyr::group_by(Dis_gene, diagnosis, Report_Variant) %>% 
      dplyr::summarize(count = n()) %>% .[-which(.$Report_Variant == "."), ]
    res$gene_report_variant <- str_c(res$Dis_gene, ":", res$Report_Variant)
    chi_table$gene_report_variant <- str_c(chi_table$Dis_gene, ":", chi_table$Report_Variant)
    
    case_ctr_sum <- as.data.frame(table(res[!duplicated(res$uniq_id), ]$diagnosis)) #统计case和control各有多少个样本
    for (i in unique(chi_table$gene_report_variant)) {
      #i <- unique(chi_table$gene_report_variant)[1]
      tryCatch({
        case_count <- chi_table[which(chi_table$gene_report_variant == i & chi_table$diagnosis == "case"), 
                                "count"]
        ctr_count <- chi_table[which(chi_table$gene_report_variant == i & chi_table$diagnosis == "control"), 
                               "count"]
        
        tmp_a <- ifelse(length(case_count$count) > 0, case_count$count, 0)
        tmp_b <- ifelse(length(ctr_count$count) > 0, ctr_count$count, 0)
        tmp_c <- case_ctr_sum[which(case_ctr_sum$Var1 == "case"), "Freq"] - tmp_a
        tmp_d <- case_ctr_sum[which(case_ctr_sum$Var1 == "control"), "Freq"] - tmp_b
        
        chi_table[which(chi_table$gene_report_variant == i),"case_variant_n"] <- tmp_a
        chi_table[which(chi_table$gene_report_variant == i),"case_variant_n_not"] <- length(case_name)-tmp_a
        chi_table[which(chi_table$gene_report_variant == i),"ctr_variant_n"] <- tmp_b
        chi_table[which(chi_table$gene_report_variant == i),"ctr_variant_n_not"] <- length(ctr_name)-tmp_b
        
        chi_table[which(chi_table$gene_report_variant == i), "fisher_p_value"] <- fisher.test(matrix(c(tmp_a, 
                                                                                                       tmp_b, tmp_c, tmp_d), nrow = 2))[["p.value"]]
      }, error = function(e) {
        cat(sprintf("#fisher.test of %s is failed!\n", i))
      })
    }
    chi_table <- dplyr::select(chi_table,-Dis_gene,-diagnosis)
    res <- merge(res, dplyr::select(chi_table ,gene_report_variant:fisher_p_value), by = "gene_report_variant")
    res <- res[which(res$ttest_p_value < 0.05 | res$fisher_p_value <0.1),]
    return(res)
  }, error = function(e) {
    cat("#Cannot return object res!\n")
  })
}


#### Step1: find essential file 1.1. find sample_info info line from data/Blood_SampleList_2020-02-27.txt
project_path <- "/cluster/home/chenhy/project/BPD/"
blood_samplelist <- read.table(file = sprintf("%s/data/Blood_SampleList_2020-02-27.txt", project_path), header = F, 
                               colClasses = "character")
colnames(blood_samplelist) <- c("file_name", "Sample_name", "FD_id", "V4", "Test_type", "Family", "Gender", "Patch", 
                                "V9", "V10", "V11", "V12")
blood_samplelist$uniq_id <- str_c(blood_samplelist$Sample_name, "_", blood_samplelist$Test_type, "_", blood_samplelist$Patch)

if(TRUE){
  if(file.exists(str_c(project_path,"src/ewas/case_file_panel.RData"))){
    load(str_c(project_path,"src/ewas/case_file_panel.RData"))
  }else{
    case_info <- openxlsx::read.xlsx(sprintf("%s/data/BPD_临床信息汇总DD20200809_tochy.xlsx", project_path),sheet=1)
    case_info <- blood_samplelist[which(blood_samplelist$Sample_name %in% unique(case_info$`样本编号2`)), ]
    dup_case_info <- case_info[case_info$Sample_name %in% case_info[duplicated(case_info$Sample_name),]$Sample_name,] %>% 
      .[order(.$Sample_name),] %>% 
      .[which(.$Test_type != "WGS" & .$Patch != "S5_trio"),] %>% .[which(.$Test_type !="WES"),]
    uniq_case_info <- case_info[!(case_info$Sample_name %in% case_info[duplicated(case_info$Sample_name),]$Sample_name),] %>% .[which(.$Test_type !="WES"),]
    case_info <- rbind(dup_case_info,uniq_case_info)
    save(case_info,file = sprintf("%s/src/ewas/case_info_panel.RData", project_path))
    case_list <- read_file(case_info)
    case_file <- merge_file(case_list)
    case_file$diagnosis <- "case"
    save(case_file,file = str_c(project_path,"src/ewas/case_file_panel.RData"))
  }
}
case_file$sample_id <- str_split(case_file$uniq_id,"_",simplify = T)[,1]
case_file <- case_file[which(case_file$sample_id %in% use_sample_id),]

# Paired samplesas control data
if(TRUE){
  if(file.exists(str_c(project_path,"src/ewas/ctr_all_file_panel.RData"))){
    load(str_c(project_path,"src/ewas/ctr_all_file_panel.RData"))
    load(str_c(project_path,"src/ewas/ctr_mat_panel.RData"))
    load(str_c(project_path,"src/ewas/ctr_list_panel.RData"))
  }else{
    ctr_info <- openxlsx::read.xlsx(sprintf("%s/data/BPD_临床信息汇总DD20200809_tochy.xlsx", project_path),sheet=2)
    ctr_info <- blood_samplelist[which(blood_samplelist$Sample_name %in% unique(ctr_info$`样本编号2`)), ]
    dup_ctr_info <- ctr_info[ctr_info$Sample_name %in% ctr_info[duplicated(ctr_info$Sample_name),]$Sample_name,] %>% 
      .[order(.$Sample_name),] %>% 
      .[-which(.$Test_type == "WGS"),] %>%  .[-which(.$uniq_id == "17F08903_YM_115th"),] %>% .[which(.$Test_type !="WES"),]
    uniq_ctr_info <- ctr_info[!(ctr_info$Sample_name %in% ctr_info[duplicated(ctr_info$Sample_name),]$Sample_name),] %>% .[which(.$Test_type !="WES"),]
    ctr_list <- rbind(dup_ctr_info,uniq_ctr_info)
    save(ctr_list,file = "ctr_list_panel.RData")
    ctr_mat <- ctr_list$uniq_id
    ctr_mat <- matrix(ctr_mat, nrow = length(ctr_mat), ncol = 1)
    save(ctr_mat,file = "ctr_mat_panel.RData")
    #Read all control data to ctr_all_file.
    ctr_all_file <- read_file(ctr_list)
    save(ctr_all_file,file = "ctr_all_file_panel.RData")
  }
}


## Step2: find re-current genes in the file with high Pathogenicity_Score-Gene
#use_gene <- read.table(file = "/cluster/apps/refseq/Gene/Clear_seq/gene.lst", col.names = "Dis_gene")
#########===============The following is a test==================================##########
# Use function cal_p no paralle system.time({ invisible( gene_list <- apply(ctr_mat[,1:4],2,cal_p) ) })
gene_list <- apply(ctr_mat,2,cal_p)[[1]]
gene_level <- gene_list[which(gene_list$ttest_p_value < 0.05),c("Dis_gene","Pathogenicity_Score-Gene","ttest_p_value")] %>% .[!duplicated(.),]
variant_level <- gene_list[which(gene_list$fisher_p_value < 0.05),c("gene_report_variant","case_variant_n" ,"case_variant_n_not","ctr_variant_n","ctr_variant_n_not","fisher_p_value")] %>% .[!duplicated(.),]
