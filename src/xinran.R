library(NetBID2)
library(tidyverse)
library(stringr)
#################
#d1 <- read.xlsx('2019BPD及对照组信息汇总-DD3.0.xlsx',1)
#d2 <- read.xlsx('2019BPD及对照组信息汇总-DD3.0.xlsx',2)
#t1 <- unique(read.delim('tableExport.csv',sep='\t',stringsAsFactors = F,fileEncoding = 'gbk'))
#t1 <- unique(t1[which(t1$Test_type %in% c('panel','WES')),])
#rownames(t1) <- t1$hospital_ID
#View(t1[d2$住院号,1])
##
load("/cluster/home/chenhy/project/BPD/src/ewas/case_info_panel.RData")
load("/cluster/home/chenhy/project/BPD/src/ewas/ctr_list_panel.RData")
ctr_list <- ctr_list[-which(ctr_list$Sample_name %in% "18F04208"),]
library(readxl)
colname_df <- openxlsx::read.xlsx("/cluster/home/chenhy/project/BPD/data/BPD_临床信息汇总DD20200809_tochy.xlsx",sheet=4,colNames=F)
d1 <- read_xlsx("/cluster/home/chenhy/project/BPD/data/BPD_临床信息汇总DD20200809_tochy.xlsx",1) %>% .[which(.$`样本编号2` %in% case_info$Sample_name),];d1 <- as.data.frame(d1)
colnames(d1) <- colname_df$X2
d2 <- read_xlsx("/cluster/home/chenhy/project/BPD/data/BPD_临床信息汇总DD20200809_tochy.xlsx",2) %>% .[which(.$`样本编号2` %in% ctr_list$Sample_name),];d2 <- as.data.frame(d2)
colnames(d2) <- colname_df$X4
d1 <- d1[which(is.na(d1$`Sample_number2`)==F),];rownames(d1) <- d1$`Sample_number2`;
d2 <- d2[which(is.na(d2$`Sample_number2`)==F),];rownames(d2) <- d2$`Sample_number2`;

#summary sample info
# d2 <- d2 %>% dplyr::mutate_at(colnames(d2),~as.factor(.x))
# tmp <- as.data.frame(summary(d2)) %>% na.omit()
# table(str_split(d2$`胎龄`,"\\+",simplify = T)[,1] %>% as.numeric(.) <32)
# sample_info <- rbind(select(d1,Household_number,Sample_number2),select(d2,Household_number,Sample_number2))
# sample_info %>% dplyr::group_by(Household_number) %>% summarise(n=n()) %>% filter(n>1) %>% dplyr::count(n)

##
d1$age <- as.numeric(gsub('(.*)\\+(.*)','\\1',d1$gestational_age))*7+
  ifelse(as.numeric(gsub('(.*)\\+(.*)','\\2',d1$gestational_age))<6,as.numeric(gsub('(.*)\\+(.*)','\\2',d1$gestational_age)),0)
d2$age <- as.numeric(gsub('(.*)\\+(.*)','\\1',d2$gestational_age))*7+
  ifelse(as.numeric(gsub('(.*)\\+(.*)','\\2',d2$gestational_age))<6,as.numeric(gsub('(.*)\\+(.*)','\\2',d2$gestational_age)),0)
d1$gender <- ifelse(d1$genders=='男',1,0)
d2$gender <- ifelse(d2$genders=='男',1,0)
##
fc <- intersect(colnames(d1),colnames(d2))
fc <- unique(c('age','gender',setdiff(fc,colnames(d1)[1:12])))
m1 <- d1[,fc];y1 <- d1$BPD_scale
M2 <- d2[,fc];y2 <- rep(0,nrow(d2))
m1 <- apply(m1,2,as.numeric)
M2 <- apply(M2,2,as.numeric)
rownames(m1) <- rownames(d1)
rownames(M2) <- rownames(d2)
X <- rbind(m1,M2)
X <- as.data.frame(X)
X$Nervous_system_infections <- NA;X$urinary_tract_infections <- NA;X$digestive_tract_infections <- NA
for (i in 1:nrow(X)) {
  X[i,"Nervous_system_infections"] <- ifelse(str_detect(X[i,"Other_systemic_infections"],"1"),1,0)
  X[i,"urinary_tract_infections"] <- ifelse(str_detect(X[i,"Other_systemic_infections"],"2"),1,0)
  X[i,"digestive_tract_infections"] <- ifelse(str_detect(X[i,"Other_systemic_infections"],"3"),1,0)
}
X <- select(X,-Other_systemic_infections)
X <- as.matrix(X)
Y <- c(y1,y2); names(Y) <- c(rownames(d1),rownames(d2))
naX1 <- apply(X,1,function(x){length(which(is.na(x)==T))})
naX2 <- apply(X,2,function(x){length(which(is.na(x)==T))})
##
u1 <- which(naX2<3) ## NA
X <- X[,u1]; 
#
naX1 <- apply(X,1,function(x){length(which(is.na(x)==T))})
u1 <- which(naX1<=0.2*ncol(X)) ## NA
X <- X[u1,]; Y <- Y[u1]
X <- as.data.frame(X)
X$age_category <- as.integer(cut(X$age,breaks = c(169,210,245,280)))
X$BPD <- as.numeric(Y[rownames(X)])-1
#############
X1 <- as.data.frame(X)
for(i in setdiff(colnames(X1),'age')) X1[,i] <- as.factor(X1[,i]) 
Y1 <- as.factor(sign(Y))
Y2 <- as.factor(Y)
X1 <- X1[,setdiff(colnames(X1),'age')]
################################# adjust for gender + CA
# remove the column with one levels
rm_cons_var <- function(df){
  rm_var <- c();two_levels <<- c()
  for (i in colnames(df)) {
    if(length(levels(df[,i]))==1){rm_var <- c(rm_var,i)}
    if(length(levels(df[,i]))==2){two_levels <<- c(two_levels,i)}
  }
  df <- df[,setdiff(colnames(df),rm_var)]
  return(df)
}

## 统计样本临床特征
library(tableone)
sample_info_stat <- X1
sample_info_stat <- within(sample_info_stat,{
  group <- NA
  group[BPD == -1] <- "no BPD"
  group[BPD %in% c(0,1,2)] <- "BPD"
})
sample_info_stat <- rm_cons_var(sample_info_stat)
if(T){
  sample_info_stat$BPD <- factor(sample_info_stat$BPD,ordered = TRUE,labels = c("no BPD","mild BPD","moderate BPD","severe BPD"))
  sample_info_stat$age_category <- factor(sample_info_stat$age_category,ordered = TRUE,labels = c("[169,210)","[210,245)","[245,280]"))
  sample_info_stat <- sample_info_stat %>% dplyr::mutate_at(two_levels,~factor(.x,ordered = TRUE,labels = c("No","Yes")))
  
  cat_table_overall <- CreateCatTable(vars = setdiff(colnames(X1),c("BPD","group")),
                                      strata = "BPD",
                                      data = sample_info_stat,includeNA = TRUE)
  des_frame1 <- as.data.frame(print(cat_table_overall,showAllLevels = TRUE,quote=F))
  des_frame1 <- des_frame1 %>% dplyr::rename(p1=p)
  des_frame1 <- des_frame1 %>% select(-test)
  cat_table_overall <- CreateCatTable(vars = setdiff(colnames(X1),c("BPD","group")),
                                      strata = "group",
                                      data = sample_info_stat,includeNA = TRUE)
  des_frame2 <- as.data.frame(print(cat_table_overall,showAllLevels = TRUE,quote=F))
  des_frame2 <- des_frame2 %>% dplyr::rename(p2=p)
  des_frame2 <- des_frame2 %>% select(-test,-level)
  des_frame <- cbind(des_frame1,des_frame2)
  select_var_df <- des_frame[which(str_detect(rownames(des_frame),"\\.$")),]
  select_var_df$p1 <- as.numeric(as.character(select_var_df$p1))
  select_var_df$p2 <- as.numeric(as.character(select_var_df$p2))
  select_var_df[which(is.na(select_var_df$p1)==1),"p1"] <- 0.001
  select_var_df[which(is.na(select_var_df$p2)==1),"p2"] <- 0.001
  u1 <- select_var_df[which(select_var_df$p1 < 0.05 | select_var_df$p2 < 0.05),] %>% 
    rownames() %>% str_replace_all(.,"\\.","")
  
  cat_table_overall <- CreateCatTable(vars = u1,
                                      strata = "BPD",
                                      data = sample_info_stat,includeNA = TRUE)
  des_frame1 <- as.data.frame(print(cat_table_overall,showAllLevels = TRUE,quote=F))
  des_frame1 <- des_frame1 %>% dplyr::rename(p1=p)
  des_frame1 <- des_frame1 %>% select(-test)
  cat_table_overall <- CreateCatTable(vars = u1,
                                      strata = "group",
                                      data = sample_info_stat,includeNA = TRUE)
  des_frame2 <- as.data.frame(print(cat_table_overall,showAllLevels = TRUE,quote=F))
  des_frame2 <- des_frame2 %>% dplyr::rename(p2=p)
  des_frame2 <- des_frame2 %>% select(-test,-level)
  des_frame <- cbind(des_frame1,des_frame2)
  write.csv(des_frame,file = "/cluster/home/chenhy/project/BPD/result/clinical_info_description.csv")
}


X1_cl <- X1[,u1];rowname.X1_cl <- rownames(X1_cl)
X1_cl <- rm_cons_var(X1_cl)
X1_cl <- X1_cl %>% mutate_at(colnames(X1_cl),~as.integer(.x))
X1_cl <- as.matrix(X1_cl)
rownames(X1_cl) <- rowname.X1_cl

#确定聚类个数
if(T){
  ## 样本聚类
  hc1 <- hclust(dist(X1_cl,method ='canberra'),method='ward.D') 
  pv <- list()
  for(j in 6:10){
    s <- cutree(hc1,k=j) #xinran's method
    print(table(list(s,Y1[names(s)])))
    pv[[j]] <- chisq.test(table(list(s,Y1[names(s)])))$p.value
  }
  use_k_sp <- c(6:10)[which.min(unlist(pv))]
  
  ## 临床信息聚类
  hc2 <- hclust(dist(t(X1_cl),method ='canberra'),method='ward.D')
  use_k_cl <- 4
}

##下面是个test
if(F){
  #Use k-modes to cluster the clinical info
  library(klaR)
  pv <- list()
  use_k <- c()
  for (i in 1:100) {
    for(j in 1:10){
      km <- kmodes(X1_cl, modes=j)
      s <- km$cluster;names(s) <- rowname.X1_cl
      #s <- cutree(hc,k=j) #xinran's method
      print(table(list(s,Y1[names(s)])))
      pv[[j]] <- chisq.test(table(list(s,Y1[names(s)])))$p.value
    }
    use_k[i] <- c(1:10)[which.min(unlist(pv))]
  }
}
#ct <- kmodes(X1_cl, modes=use_k)$cluster;names(ct) <- rowname.X1_cl
ct <- cutree(hc1,k=use_k_sp)
##
draw_mat <- X[,colnames(X1_cl)]
draw_mat <- t(draw_mat[names(ct[order(ct,decreasing = T)]),hc2$order])
draw_mat <- rbind(draw_mat[setdiff(rownames(draw_mat),'BPD'),],'BPD'=as.numeric(Y1[colnames(draw_mat)])-1,'BPD_detail'=as.numeric(Y2[colnames(draw_mat)])-1)
#draw_mat <- rbind(draw_mat,'BPD_detail'=as.numeric(Y2[colnames(draw_mat)])-1)
ct <- ct[colnames(draw_mat)]
sample_class <- paste0('C',ct)
names(sample_class) <- names(ct)

if(T){
  library(dendextend)
  dend_row = hclust(dist(X1_cl,method ='canberra'),method='ward.D')
  dend_row1 = color_branches(dend_row, k = 6)
  dend_col = hclust(dist(t(X1_cl),method ='canberra'),method='ward.D')
  dend_col1 = color_branches(dend_col, k = 3,col =  c("#FC4E07","#E7B800","#2E9FDF"))
  #kmode_cl <- kmodes(X1_cl, modes=6)
  #热图注释
  df1 <- data.frame(type = paste0("C", cutree(dend_row,k = 6)))
  col1 <- colorRampPalette(brewer.pal(6, "Set3"))(6)
  ha1 <-  HeatmapAnnotation(df = df1,
                            col = list(type = c("C1" = col1[1], "C2" = col1[2], "C3" = col1[3], 
                                                "C4" = col1[4],"C5" = col1[5],"C6" = col1[6])),
                            which = "row")
  
  df2 <- data.frame(type = paste0("S", cutree(dend_col,k = 3)))
  col2 <- c("#2E9FDF","#FC4E07","#E7B800")
  ha2 <-  HeatmapAnnotation(df = df2,
                            col = list(type = c("S1" = col2[1], "S2" = col2[2], "S3" = col2[3])),
                            which = "column")
  
  draw(ha2)
  pdf(file="/cluster/home/chenhy/project/BPD/result/plot/clinical_cha_pheatmap.pdf",width = 9,height = 15)
  Heatmap(X1_cl, col =c("#E6E6FA", "#9370DB"),
          #left_annotation = ha1,
          top_annotation = ha2,
          name = "foo",#同时绘制多个热图，这名字也可以作为唯一的id。稍后，我们可以使用这个名称到相应的热图中添加更多的图形 (具体看[Heatmap Decoration]vignette).
          heatmap_legend_param = list(title = "legend"),
          column_title = "Clinical character", 
          row_title = "Sample cluster",
          show_column_names = T,
          show_row_names = F,
          column_title_side = "top",
          #row_names_gp = gpar(fontsize = 20) #设置字体大小
          column_names_gp = gpar(fontsize = 9),
          #column_names_side = "top",
          #column_title_gp = gpar(fontsize = 20, fontface = "bold"),
          show_column_dend = T,#聚类树显示
          row_dend_side = "left",#聚类树位置
          column_dend_side = "top",
          row_dend_width = unit(1.5, "cm"),
          column_dend_height = unit(1.5, "cm"),
          clustering_distance_rows = function(m) dist(m,method ='canberra'),#预定义选项可用的值是dist()函数中支持的方法和pearson, spearman and kendall. NA 值在预定义选项的聚类中是被忽略的但会抛出 warnings (see example in Colors section).
          clustering_method_rows = "ward.D", #可以通过' clustering_method_rows '和' clustering_method_columns '指定实现分层聚类的方法。可用的方法是' hclust() '函数中支持的方法
          clustering_method_columns = "ward.D",
          #cluster_rows = as.dendrogram(diana(mat)),#默认情况下，聚类由' hclust() '执行。但是您也可以通过将“cluster_rows”或“cluster_columns”指定为“hclust”或“dendrogram”对象来利用其他方法生成的聚类结果。
          #cluster_columns = as.dendrogram(agnes(t(mat)))
          #cluster_rows = dend_row,#给树状图上颜色
          cluster_columns = dend_col1,
          #split = 6,
          split = paste0("C", cutree(dend_row,k = 6)), #按行分割热图
          gap = unit(0.5,"mm")#设置间隔高度
  )
  dev.off()
}
#
library(Cairo)
if(T){
  CairoPDF("/cluster/home/chenhy/project/BPD/result/plot/f1.pdf",width=8,height=12,family='GB1')
  par(mar=c(15,6,4,6))
  image(draw_mat,col=(colorRampPalette(c('white','red'))(5)),breaks=c(-1,0,1,2,3,4),bty='o',xaxt='n',yaxt='n')
  pp <- par()$usr ## get to know the xy position of the figure
  x_line_pos <- seq(pp[1],pp[2],length.out=nrow(draw_mat)+1) ## x position for vertical line
  y_line_pos <- seq(pp[3],pp[4],length.out=ncol(draw_mat)+1) ## y position for each horizon line
  x_mid_pos  <- x_line_pos[1:(length(x_line_pos)-1)]/2+x_line_pos[2:length(x_line_pos)]/2 ## x position for column mid position
  y_mid_pos  <- y_line_pos[1:(length(y_line_pos)-1)]/2+y_line_pos[2:length(y_line_pos)]/2 ## y position for row mid position
  text(y=min(y_line_pos)-1e-2,x=x_mid_pos,rownames
       (draw_mat),xpd=T,cex=0.6,adj=1,srt=90,family='GB1')
 
  all_c <- unique(sample_class)
  col1 <- colorRampPalette(brewer.pal(8, "Set3"))(length(all_c)); names(col1) <- all_c
  gen <- sample_class[colnames(draw_mat)] ## gene category
  w1 <- cumsum(table(gen)[unique(gen)])
  segments(x0=min(x_line_pos),x1=max(x_line_pos),y0=y_line_pos[c(1,w1+1)],y1=y_line_pos[c(1,w1+1)])
  dw <- (pp[2]-pp[1])/20
  rect(xleft=pp[1]-dw,xright=pp[1]-dw*0.1,ybottom=y_line_pos[1:(length(y_line_pos)-1)],
       ytop=y_line_pos[2:(length(y_line_pos))],xpd=T,col=col1[gen],border=NA)
  text(x=pp[1]-2*dw,y=y_line_pos[1+w1]/2+y_line_pos[1+c(0,w1[1:(length(w1)-1)])]/2,all_c,srt=90,xpd=T,adj=0.5)
  dev.off()
}
###############################
table(list(ct[names(Y1)],Y1))
table(list(ct[names(Y1)],Y2))


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

#######################get genetic data
load("/cluster/home/chenhy/project/BPD/src/ewas/case_file_panel.RData")
load("/cluster/home/chenhy/project/BPD/src/ewas/ctr_all_file_panel.RData")
ctr_file <- merge_file(ctr_all_file)
case_file$V1 <- str_split(case_file$uniq_id,"_",simplify = T)[,1]
case_file$Dis_variant <- str_c(case_file$Report_Gene,":",case_file$Report_Variant)
ctr_file$V1 <- str_split(ctr_file$uniq_id,"_",simplify = T)[,1]
ctr_file$Dis_variant <- str_c(ctr_file$Report_Gene,":",ctr_file$Report_Variant)
gene_ttest <- read.csv("/cluster/home/chenhy/project/BPD/result/ewas/gene_ttest_paired_panel.csv",stringsAsFactors = F) %>% .[order(.$p),]
variant_test <- read.csv("/cluster/home/chenhy/project/BPD/result/ewas/variant_fisher_paired_panel.csv",stringsAsFactors = F) %>% .[order(.$p),]
#Dis_gene <- intersect(gene_ttest$Dis_gene,variant_test$Dis_gene)
length(unique(case_file$V1));length(unique(ctr_file$V1))
table(unique(case_file$V1) %in% rownames(X1));table(unique(ctr_file$V1) %in% rownames(X1))
Dis_gene <- unique(gene_ttest$Dis_gene) %>% head(30)
Dis_variant <- unique(variant_test$Dis_variant) %>% head(30)
if(F){
  case_file <- case_file[,c("V1","Dis_gene","Pathogenicity_Score-Gene")] #case_file$Dis_gene %in% unique(db1$Dis_gene) & 
  ctr_file <- ctr_file[,c("V1","Dis_gene","Pathogenicity_Score-Gene")]
  db1 <- rbind(case_file,ctr_file)
  db1 <- dplyr::rename(db1,V2=Dis_gene,V3=`Pathogenicity_Score-Gene`) %>% .[!duplicated(.[,c("V1","V2")]),];
  ####################### load genetic information
  #db1 <- read.delim('E:/写写文章/GTLC/GTLC/CCGT_result_201912/data/sample2gene2level.txt',stringsAsFactors = F,header=F)
  db2 <- db1[which(db1$V1 %in% names(Y1)),]
  u_g <- names(which(table(db2$V2)>=3));u_s <- names(Y1)
  db2 <- db2[which(db2$V2 %in% u_g),]
  sample2gene <- matrix(0,nrow=length(u_s),ncol=length(u_g))
  rownames(sample2gene) <- u_s;colnames(sample2gene) <- u_g
  for(i in 1:nrow(db2)){
    sample2gene[db2$V1[i],db2$V2[i]] <- db2$V3[i]
  }
  apply(sample2gene,2,function(x){sum(is.na(x))})
  ## get candidate genes
  get_candidate <- function(g1,g0,pv_thre=0.1){
    pv_g <- apply(sample2gene,2,function(x){
      x1 <- x[which(u_s %in% g1)]
      x0 <- x[which(u_s %in% g0)]
      tryCatch({
        t.test(x1,x0)$p.value*sign(mean(x1)-mean(x0))
      },error = function(e) {
        cat("#Scores of gene are essentially constant!\n")
      })
    })
    
    tryCatch({
      names(pv_g) <- u_g
      pv_g <- unlist(pv_g)
      pv_g <- sort(pv_g)
    },error = function(e) {
      cat("#pv_g is NULL!\n")
    })
    return(pv_g)
    #pv_g[which(abs(pv_g)<pv_thre)]
  }
  check_candidate <- function(g1,g0,use_gene){
    x <- sample2gene[,use_gene]
    x1 <- x[which(u_s %in% g1)]
    x0 <- x[which(u_s %in% g0)]
    return(list(db_g=db_g,x1,x0,
                tryCatch({
                  pv=t.test(x1,x0)$p.value*sign(mean(x1)-mean(x0))
                },error = function(e) {
                  cat("#Scores of gene are essentially constant!\n")
                })
    ))
  }
  ## BPD Vs. non-BPD
  cd <- list()
  g1 <- names(Y1)[which(Y1==1)];g0 <- names(Y1)[which(Y1==0)];
  cd[['BPD vs N-BPD']] <- get_candidate(g1,g0,pv_thre=0.05)
  g1 <- names(Y1)[which(Y2==3)];g0 <- names(Y1)[which(Y1==0)];
  cd[['BPD3 vs N-BPD']] <- get_candidate(g1,g0,pv_thre=0.05)
  g1 <- names(Y1)[which(Y2==2)];g0 <- names(Y1)[which(Y1==0)];
  cd[['BPD2 vs N-BPD']] <- get_candidate(g1,g0,pv_thre=0.05)
  g1 <- names(Y1)[which(Y2==1)];g0 <- names(Y1)[which(Y1==0)];
  cd[['BPD1 vs N-BPD']] <- get_candidate(g1,g0,pv_thre=0.05)
  ##
  for(i in 1:6){
    g1 <- names(ct)[which(ct==i)];g0 <- names(ct)[which(ct!=i)];
    cd[[sprintf('Class %s vs others',i)]] <- get_candidate(g1,g0,pv_thre=0.05)
  }
  for(i in 3:6){
    g1 <- names(ct)[which(ct==i)];g0 <- names(Y1)[which(Y1==0)];
    cd[[sprintf('Class %s vs N-BPD',i)]] <- get_candidate(g1,g0,pv_thre=0.05)
  }
  ## candidate genes
  cd1 <- lapply(cd,function(x){x[which(x>0 & x<0.1)]}) ## candidate
  # cd2 <- lapply(cd,function(x){sort(x[which(x<0 & x> -0.01)],decreasing = T)}) # need to remove
  ### cluster sample by gene information
  M1 <- sample2gene
  ## feature selection
  pv1 <- apply(M1,2,function(x){
    t.test(x~Y1,alternative='less')$p.value
  })
  ##
  u1 <- which(pv1<=0.01) ## sig
  M2 <- M1[names(Y1),u1]
  #M2 <- M2[,which(colnames(M2) %in% Dis_gene)]
}

if(T){
  case_file <- case_file[,c("V1","Dis_gene","Pathogenicity_Score-Gene")] #which(case_file$V1 %in% rownames(X1))
  ctr_file <- ctr_file[,c("V1","Dis_gene","Pathogenicity_Score-Gene")] #which(ctr_file$Dis_gene %in% Dis_gene)
  db1 <- rbind(case_file,ctr_file)
  db1 <- dplyr::rename(db1,V2=Dis_gene,V3=`Pathogenicity_Score-Gene`) %>% .[!duplicated(.),];
  Y1 <- Y1[names(Y1) %in% unique(db1$V1)]
  Y2 <- Y2[names(Y2) %in% unique(db1$V1)]
  M2_gene <- matrix(0,nrow = length(Y1),ncol = length(Dis_gene))
  rownames(M2_gene) <- names(Y1);colnames(M2_gene) <- Dis_gene
  for (i in rownames(M2_gene)) {
    for (j in colnames(M2_gene)) {
      tmp <- db1[which(db1$V1==i & db1$V2==j),"V3"]
      M2_gene[i,j] <- ifelse(length(tmp)==0,0,tmp)
    }
  }
}


if(F){
  case_file <- case_file[,c("V1","Dis_variant","Pathogenicity_Score-Gene")] #which(case_file$V1 %in% rownames(X1))
  ctr_file <- ctr_file[,c("V1","Dis_variant","Pathogenicity_Score-Gene")] #which(ctr_file$Dis_variant %in% Dis_variant)
  db1 <- rbind(case_file,ctr_file)
  db1 <- dplyr::rename(db1,V2=Dis_variant,V3=`Pathogenicity_Score-Gene`) %>% .[!duplicated(.),];
  Y1 <- Y1[names(Y1) %in% unique(db1$V1)]
  Y2 <- Y2[names(Y2) %in% unique(db1$V1)]
  M2_variant <- matrix(0,nrow = length(Y1),ncol = length(Dis_variant))
  rownames(M2_variant) <- names(Y1);colnames(M2_variant) <- Dis_variant
  for (i in rownames(M2_variant)) {
    for (j in colnames(M2_variant)) {
      tmp <- db1[which(db1$V1==i & db1$V2==j),"V3"]
      M2_variant[i,j] <- ifelse(length(tmp)==0,0,tmp)
    }
  }
}

if(T){
  M2 <- M2_gene
  #M2 <- M2_variant
  #heatmap(M2,scale='none')
  ##
  hc <- hclust(dist(M2),method ='ward.D')
  hc1 <- hclust(dist(t(M2)),method ='ward.D')
  pv <- list()
  for(i in 5:10){
    pv[[i]] <- chisq.test(table(list(cutree(hc,k=i)[names(Y1)],Y2[names(Y1)])))$p.value
  }
  use_k <- c(5:10)[which.min(unlist(pv))]
  ct1 <- cutree(hc,k=use_k)
  table(list(ct[names(Y1)],ct1[names(Y1)]))
  table(list(Y1,ct1))
  table(list(Y2,ct1))
  
}

source("/cluster/home/chenhy/project/BPD/src/draw_heatmap.R")
if(T){
  pdf('/cluster/home/chenhy/project/BPD/result/plot/f2.pdf',width=12,height=6)
  #ct2 <- sprintf('B%sC%s',Y1,ct1);names(ct2) <- names(ct1);ct2 <- sort(ct2)
  #draw_mat <- M2[names(ct2),hc1$order]
  #par(mar=c(4,3,2,10))
  #draw_heatmap(t(draw_mat),gene_category =NULL,phenotype_category=ct2,row_text=T,col_cex=0.3,
  #             min_val=0,max_val=1,legend_col=(brewer.pal(9,'Reds')))
  ct2 <- sprintf('B%sC%s',Y2[names(ct1)],ct1[names(ct1)]);names(ct2) <- names(ct1);ct2 <- sort(ct2)
  draw_mat <- t(M2[names(ct2),hc1$order])
  par(mar=c(4,3,2,12))
  tmp2 <- sprintf("BPD%s",Y2);names(tmp2) <- names(Y2)
  draw_heatmap(draw_mat,gene_category =NULL,phenotype_category=tmp2[names(ct2)],row_text=T,col_cex=0.9,
               min_val=0,max_val=1,legend_col=c('white',brewer.pal(9,'Reds')))
  pp <- par()$usr ## get to know the xy position of the figure
  draw_mat <- t(draw_mat)
  x_line_pos <- seq(pp[1],pp[2],length.out=nrow(draw_mat)+1) ## x position for vertical line
  y_line_pos <- seq(pp[3],pp[4],length.out=ncol(draw_mat)+1) ## y position for each horizon line
  x_mid_pos  <- x_line_pos[1:(length(x_line_pos)-1)]/2+x_line_pos[2:length(x_line_pos)]/2 ## x position for column mid position
  y_mid_pos  <- y_line_pos[1:(length(y_line_pos)-1)]/2+y_line_pos[2:length(y_line_pos)]/2 ## y position for row mid position
  col1 <- colorRampPalette(brewer.pal(8, "Set3"))(length(all_c)); names(col1) <- all_c
  cc <- col1[paste0('C',ct[names(ct2)])]
  rect(xleft=x_line_pos[1:length(cc)],xright=x_line_pos[2:(length(cc)+1)],ytop=min(y_line_pos),
       ybottom=min(y_line_pos)-(y_line_pos[2]-y_line_pos[1])/2,col = cc,xpd=T,border = NA)
  text(x=x_mid_pos,y=min(y_line_pos),sprintf('C%s',ct[names(ct2)]),cex=0.15,xpd=T,pos=1)
  dev.off()
}

##
if(T){
  pdf('/cluster/home/chenhy/project/BPD/result/plot/f3.pdf',width=10,height=5)
  ct2 <- sprintf('%sC%s',ifelse(Y1==1,'Y','N'),ct1);names(ct2) <- names(ct1);ct2 <- sort(ct2)
  draw_mat <- M2[names(ct2),hc1$order]
  par(mar=c(4,3,2,12))
  draw_heatmap(t(draw_mat),gene_category =NULL,phenotype_category=ct2,row_text=T,col_cex=0.35,
               min_val=0,max_val=1,legend_col=c('white',brewer.pal(9,'Reds')))
  x_line_pos <- seq(pp[1],pp[2],length.out=nrow(draw_mat)+1) ## x position for vertical line
  y_line_pos <- seq(pp[3],pp[4],length.out=ncol(draw_mat)+1) ## y position for each horizon line
  x_mid_pos  <- x_line_pos[1:(length(x_line_pos)-1)]/2+x_line_pos[2:length(x_line_pos)]/2 ## x position for column mid position
  y_mid_pos  <- y_line_pos[1:(length(y_line_pos)-1)]/2+y_line_pos[2:length(y_line_pos)]/2 ## y position for row mid position
  col1 <- colorRampPalette(brewer.pal(8, "Set3"))(length(all_c)); names(col1) <- all_c
  cc <- col1[paste0('C',ct[names(ct2)])]
  rect(xleft=x_line_pos[1:length(cc)],xright=x_line_pos[2:(length(cc)+1)],ytop=min(y_line_pos)-(y_line_pos[2]-y_line_pos[1])/4,
       ybottom=min(y_line_pos)-(y_line_pos[2]-y_line_pos[1])/2,col = cc,xpd=T,border = NA)
  text(x=x_mid_pos,y=min(y_line_pos),sprintf('C%s',ct[names(ct2)]),cex=0.12,xpd=T,pos=1)
  text(x=x_mid_pos,y=min(y_line_pos)-(y_line_pos[2]-y_line_pos[1])/8,sprintf('B%s',Y2[names(ct2)]),cex=0.1,xpd=T)
  dev.off()
}





##############################################################
# Model Based Clustering for Mixed Data
library(clustMD)
M2 <- as.data.frame(M2);#M2.rowname <- rownames(M2)
sample_info_stat <- X1
# Start categorical variables at 1 rather than 0
sample_info_stat <- sample_info_stat %>% mutate_at(colnames(sample_info_stat),~as.numeric(.x));
rownames(sample_info_stat) <- rownames(X1)
choose_var <- apply(sample_info_stat,2,function(x) sum(is.na(x))) < 12
sample_info_stat <- sample_info_stat[,choose_var]
sample_info_stat <- na.omit(sample_info_stat)
bpd_scale <- sample_info_stat$BPD
sample_info_stat <- dplyr::select(sample_info_stat,-BPD)
if(T){
  tmp <- sample_info_stat
  tmp <- tmp %>% mutate_at(colnames(tmp),~as.factor(.x))
  binary_var <- c();ordinal_var <- c()
  for (i in colnames(tmp)) {
    if(length(levels(tmp[,i]))==2){binary_var <- c(binary_var,i)}
    if(length(levels(tmp[,i]))>2){ordinal_var <- c(ordinal_var,i)}
  }
}
sample_info_stat <- sample_info_stat[,c(binary_var,ordinal_var)]

# Standardise continuous variables
M3 <- scale(M2)
sampleinfo_gene <- merge(M3,sample_info_stat,by=0,all.y=T);sampleinfo_gene.rowname <- sampleinfo_gene$Row.names
sampleinfo_gene <- sampleinfo_gene[,-1];rownames(sampleinfo_gene) <- sampleinfo_gene.rowname
sampleinfo_gene$BPD <- bpd_scale

res <- clustMD(X = sampleinfo_gene, G = ncol(sample_info_stat), CnsIndx = ncol(M3), OrdIndx = ncol(sample_info_stat)+ncol(M3), 
               Nnorms = 20000, MaxIter = 500, model = "EII", store.params = FALSE, scale = TRUE,
               startCL = "kmeans", autoStop= TRUE, ma.band=30, stop.tol=0.0001)

res2 <- clustMDparallel(X = sampleinfo_gene, G = ncol(sample_info_stat), CnsIndx = ncol(M3), OrdIndx = ncol(sample_info_stat)+ncol(M3), 
                Nnorms = 20000, MaxIter = 500, model = c("EVI", "EII", "VII"), store.params = FALSE, scale = FALSE,
                startCL = "kmeans", autoStop= TRUE, ma.band=30, stop.tol=0.0001)


