library(party) 
set.seed(2)
data(iris)
#随机抽取训练集和测试集
index<-sample(1:nrow(iris),100)
train<-iris[index,]
test<-iris[-index,]
#建立C4.5决策树模型
ctree.model <- ctree(Species ~ ., data <- train)  
#输出决策树图
plot(ctree.model, type <- "simple")  
data1 <- matrix(sample(1:3,size=1000,replace = T),nrow=100)
data2 <- matrix(sample(4:6,size=1000,replace = T),nrow=100)
data <- rbind(data1,data2)
set.seed(1)
x <- rbind(matrix(rbinom(250, 2, 0.25), ncol = 5),
           matrix(rbinom(250, 2, 0.75), ncol = 5))
colnames(x) <- c("a", "b", "c", "d", "e")
## run algorithm on x:
(cl <- kmodes(x, 2))
## and visualize with some jitter:
plot(jitter(x), col = cl$cluster)
points(cl$modes, col = 1:5, pch = 8)


hc1 <- hclust(dist(X1_cl,method ='canberra'),method='ward.D') 
hc2 <- hclust(dist(t(X1_cl),method ='canberra'),method='ward.D')
o1 <- h1$order
o2 <- h2$order
use_mat1 <- X1_cl[o1,o2]

hcd <- as.dendrogram(h1)
plot(hcd,type = "rectangle", ylab = "",xlim=c(0.5,length(h1$order)+0.5),xaxs='i',yaxs='i',yaxt='n',
     ylim=c(-8,0.5+max(h1$height)))
hcd <- as.dendrogram(h2)
plot(hcd,type = "rectangle", ylab = "",xlim=c(0.5,length(h2$order)+0.5),xaxs='i',yaxs='i',yaxt='n',
     ylim=c(-8,0.5+max(h2$height)))



hcd <- as.dendrogram(hc2)
# vector of colors labelColors = c('red', 'blue', 'darkgreen', 'darkgrey',
# 'purple')
labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951")
# cut dendrogram in 4 clusters
clusMember = cutree(hc2, 4)
# function to get color labels
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}
# using dendrapply
clusDendro = dendrapply(hcd, colLab)
# make plot
plot(clusDendro,horiz = T,xlab="",ylab = "",type = "rectangle",xaxs='i',yaxs='i',yaxt='n',xaxt='n',axes = F, leaflab = "none")

plot(clusDendro,horiz = T,xlab="",ylab = "",type = "rectangle",xaxs='i',yaxs='i',yaxt='n',xaxt='n',
     xlim=c(0.5,length(h2$order)+0.5),ylim=c(-900,500+max(h2$height)))



library(ComplexHeatmap)
library(circlize)

set.seed(123) #郑宝童简书上，提供了一个关于set.seed()函数的解析

mat = cbind(rbind(matrix(rnorm(16, -1), 4), matrix(rnorm(32, 1), 8)),
            rbind(matrix(rnorm(24, 1), 4), matrix(rnorm(48, -1), 8)))

# permute the rows and columns 排列行和列
mat = mat[sample(nrow(mat), nrow(mat)), sample(ncol(mat), ncol(mat))]

rownames(mat) = paste0("R", 1:12)
colnames(mat) = paste0("C", 1:10)
Heatmap(mat)
mat2 = mat
mat2[1, 1] = 100000
Heatmap(mat2, col = colorRamp2(c(-3, 0, 3), c("green", "white", "red")), 
        cluster_rows = FALSE, cluster_columns = FALSE)


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


result <- get_dist(t(M2))
result_hc <- hclust(d = result, method = "ward.D")
tmp <- as.data.frame(cbind(M2, cluster = cutree(result_hc,k = 3)))
#fviz_dend(result_hc, cex = 0.6)蓝"#2E9FDF", 浅棕"#00AFBB", 橘黄"#E7B800", 红"#FC4E07"
fviz_dend(result_hc, k = 3, 
          cex = 0.5, 
          k_colors = c("#2E9FDF", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, 
          show_labels = F,
          rect = TRUE
)


result <- get_dist(M2)
result_hc <- hclust(d = result, method = "ward.D")
tmp <- as.data.frame(cbind(M2, cluster = cutree(result_hc,k = 3)))
#fviz_dend(result_hc, cex = 0.6)蓝"#2E9FDF", 浅棕"#00AFBB", 橘黄"#E7B800", 红"#FC4E07"
fviz_dend(result_hc, k = 4, 
          cex = 0.5, 
          k_colors = c("#2E9FDF", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, 
          show_labels = F,
          rect = TRUE
)
draw(ha1)


ha1 = HeatmapAnnotation(pt = anno_points(1:10, gp = gpar(col = rep(2:3, each = 5)), 
                                         height = unit(2, "cm")), show_annotation_name = FALSE)
ha2 = HeatmapAnnotation(ln = anno_lines(cbind(1:10, 10:1), gp = gpar(col = 4:5, lty = 1:2),
                                        height = unit(2, "cm")), show_annotation_name = FALSE)
m = matrix(rnorm(100), 10)
ht_list = Heatmap(m, name = "mat1", top_annotation = ha1) + 
  Heatmap(m, name = "mat2", top_annotation = ha2) +
  Heatmap(m[, 1], name = "mat3", 
          top_annotation = HeatmapAnnotation(
            summary = anno_summary(gp = gpar(fill = 2:3))
          ), width = unit(1, "cm"))
lgd_list = list(
  Legend(labels = c("red", "green"), title = "pt", type = "points", pch = 16, 
         legend_gp = gpar(col = 2:3)),
  Legend(labels = c("darkblue", "lightblue"), title = "ln", type = "lines", 
         legend_gp = gpar(col = 4:5, lty = 1:2)),
  Legend(labels = c("group1", "group2"), title = "km", type = "boxplot",
         legend_gp = gpar(fill = 2:3))
)
draw(ht_list, ht_gap = unit(7, "mm"), row_km = 2, annotation_legend_list = lgd_list)


for (i in candidate_driver) {
  i <- candidate_driver[1]
  for (j in 1:length(gse_list)) {
    group_vars <- c("disease status","ASD.CTL","diagnosis")
    assign(str_c("p",j),draw_boxplot(gse_list[[j]],group_vars,i))
  }
  ggarrange(p1,p2,p3,ncol=3,nrow = 1,align = "hv",hjust = 10,vjust = 10)
}

draw_network(gse_list[["GSE64018"]],"ITGA6")
setwd("/cluster/home/chenhy/")

library(qpdf)
setwd("/cluster/home/chenhy/project/autism/result/plot/TYMP/")
pdf_files <- list.files("/cluster/home/chenhy/project/autism/result/plot/TYMP/")
pdf_combine(pdf_files,output = "joined.pdf")




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
  for(j in 10:15){
    s <- cutree(hc1,k=j) #xinran's method
    print(table(list(s,Y1[names(s)])))
    pv[[j]] <- chisq.test(table(list(s,Y1[names(s)])))$p.value
  }
  use_k_sp <- c(10:15)[which.min(unlist(pv))]
  
  ## 临床信息聚类
  hc2 <- hclust(dist(t(X1_cl),method ='canberra'),method='ward.D')
  use_k_cl <- 3
}
ct <- cutree(hc1,k=use_k_sp)
draw_mat <- X[rownames(X1_cl),colnames(X1_cl)]
draw_mat <- t(draw_mat[names(ct[order(ct,decreasing = T)]),hc2$order])
#draw_mat <- rbind(draw_mat[setdiff(rownames(draw_mat),'BPD'),],'BPD'=as.numeric(Y1[colnames(draw_mat)])-1,'BPD_detail'=as.numeric(Y2[colnames(draw_mat)])-1)
# ct <- ct[colnames(draw_mat)]
# sample_class <- paste0('C',ct)
# names(sample_class) <- names(ct)
draw_mat <- t(draw_mat)
dend_row = hclust(dist(draw_mat,method ='canberra'),method='ward.D')
dend_row1 = color_branches(dend_row, k = use_k_sp)
dend_col = hclust(dist(t(draw_mat),method ='canberra'),method='ward.D')
dend_col1 = color_branches(dend_col, k = use_k_cl,col =  c("#FC4E07","#2E9FDF","#E7B800"))#蓝"#2E9FDF", 浅棕"#00AFBB", 橘黄"#E7B800", 红"#FC4E07"
## 构建热图注释
df1 <- data.frame(cluster = paste0("C", cutree(dend_row,k = use_k_sp)),
                  BPD_scale = paste0("BPD", Y),
                  BPD = paste0("BPD",Y1))
col1 <- colorRampPalette(brewer.pal(use_k_sp, "Set3"))(use_k_sp)
col2 <- colorRampPalette(brewer.pal(4, "Set1"))(4)
ha1 <-  HeatmapAnnotation(df = df1,
                          col = list(cluster = c("C1" = col1[1], "C2" = col1[2], "C3" = col1[3], 
                                                 "C4" = col1[4], "C5" = col1[5], "C6" = col1[6], 
                                                 "C7" = col1[7], "C8" = col1[8], "C9" = col1[9],
                                                 "C10" = col1[10], "C11" = col1[11], "C12" = col1[12], "C13" = col1[13]),
                                     BPD_scale = c("BPD1" = col2[1], "BPD2" = col2[2], "BPD3" = col2[3],"BPD0" = col2[4])),
                          which = "row")

df2 <- data.frame(type = paste0("S", cutree(dend_col,k = 3)))
col2 <- c("#2E9FDF","#FC4E07","#E7B800")
ha2 <-  HeatmapAnnotation(df = df2,
                          name = "type",
                          col = list(type = c("S1" = col2[1], "S2" = col2[2], "S3" = col2[3])),
                          which = "column",
                          annotation_legend_param = list(type=list(direction = "horizontal"))
)


#draw(ha2)
#pdf(file="/cluster/home/chenhy/project/BPD/result/plot/clinical_cha_pheatmap.pdf",width = 9,height = 15)
Heatmap(draw_mat, col =c("#E6E6FA", "#9370DB"),
        left_annotation = ha1,
        top_annotation = ha2,
        name = "foo",#同时绘制多个热图，这名字也可以作为唯一的id。稍后，我们可以使用这个名称到相应的热图中添加更多的图形 (具体看[Heatmap Decoration]vignette).
        km=2,
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
        row_dend_width = unit(2.5, "cm"),
        column_dend_height = unit(1, "cm"),
        clustering_distance_rows = function(m) dist(m,method ='canberra'),#预定义选项可用的值是dist()函数中支持的方法和pearson, spearman and kendall. NA 值在预定义选项的聚类中是被忽略的但会抛出 warnings (see example in Colors section).
        clustering_method_rows = "ward.D", #可以通过' clustering_method_rows '和' clustering_method_columns '指定实现分层聚类的方法。可用的方法是' hclust() '函数中支持的方法
        clustering_method_columns = "ward.D",
        #cluster_rows = as.dendrogram(diana(mat)),#默认情况下，聚类由' hclust() '执行。但是您也可以通过将“cluster_rows”或“cluster_columns”指定为“hclust”或“dendrogram”对象来利用其他方法生成的聚类结果。
        #cluster_columns = as.dendrogram(agnes(t(mat)))
        #cluster_rows = dend_row,#给树状图上颜色
        cluster_columns = dend_col1,
        #split = 6,
        #split = paste0("BPD",Y1), #按行分割热图
        gap = unit(0.5,"mm")#设置间隔高度
)








cluster_comp <- data.frame()
source("/cluster/home/chenhy/project/BPD/src/logistic_fun.R")
for (i in 1:use_k_sp) {
  i <- 1
  tmp_mat <- cl_mat[which(cl_mat$cluster %in% c(0,i)),]
  tmp_mat$group <- 0
  tmp_mat[which(tmp_mat$cluster == i),"group"] <- 1
  for (j in colnames(X1_cl)) {
    #j <- colnames(X1_cl)[43]
    res <- logitUniVar(tmp_mat, group = "group", var = j)
    if(res["p.val"] %>% as.numeric(.) < 0.05){
      log_p <- res["p.val"] %>% as.numeric(.) %>% -log10(.)
      OR <- res["OR"] %>% as.numeric(.)
      tmp.df <- data.frame(bpd = str_c("BPD",i),cl = j,log_p = log_p,OR = OR)
    }
    cluster_comp <- rbind(cluster_comp,tmp.df)
  }
  
}
cluster_comp <- cluster_comp[!duplicated(cluster_comp),]







# Load package
library(networkD3)

# Create fake data
src <- c("A", "A", "A", "A",
         "B", "B", "C", "C", "D")
target <- c("B", "C", "D", "J",
            "E", "F", "G", "H", "I")
networkData <- data.frame(src, target)

# Plot
simpleNetwork(network_df)


# Load data
data(MisLinks)
data(MisNodes)

network_df$source <- str_split(network_df$start,"D",simplify = T)[,2]
network_df$target <- as.factor(network_df$end) %>% as.numeric()
network_df$value <- network_df$weight
network_df$group <- network_df$source
network_df$name <- network_df$start
links <- select(network_df,source,target,value)
nodes <- rbind(select(network_df,name,group),select(network_df,end,group) %>% dplyr::rename(name = end)) %>% .[!duplicated(.),]
# Plot
forceNetwork(Links = links, Nodes = nodes,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.8)

forceNetwork(Links = MisLinks, Nodes = MisNodes,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.8)



forceNetwork(Links = links, Nodes = nodes,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", Nodesize = "size",
             fontSize = 30,zoom = T,opacity = 0.8)



dat <- cluster_comp %>% mutate(y=ifelse(OR>1,OR,-1/OR),color=ifelse(OR>1,"blue","red")) %>% select(bpd,cl,y,color) #%>% filter(bpd=="BPD1") 
ggplot(dat,aes(x = cl,y = y,fill=color))+
  geom_bar(stat ="identity",width = 0.6,position = "dodge")+
  facet_grid(.~bpd)+
  xlab("Clinical characteristics")+
  ylab("OR")+
  scale_y_continuous(breaks=c(-150,-100,-50,0,50,100,150,200,250,300),labels=c("1/150","1/100","1/50",0,50,100,150,200,250,300))+
  #geom_text(aes(label = dat$cl),position=position_dodge(width = 0.9),size = 5,vjust = 0)+
  theme(legend.position='none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_line(size = 0.08, colour = '#1391FF'), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white",color="darkblue"),axis.line = element_line(colour = "darkblue",size = 0.1)
        #axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank() #去掉纵坐标轴
        )+
  coord_flip()



library(factoextra)
hc1 <- hclust(dist(X1_cl,method ='canberra'),method='ward.D')
fviz_nbclust(X1_cl, FUN = kmeans, method = "wss",diss = dist(X1_cl, method = "canberra")) +
  geom_vline(xintercept = 3, linetype = 2)

library(igraph)
my_edges <- network_df[,1:2]
tmp <- graph.data.frame(my_edges)
tmp;summary(tmp)
plot(tmp, layout=layout.kamada.kawai)
subnet <- as_graphnel(tmp)
subnet


intersect(network_df[which(network_df$start == str_c("C",2)),"end"],variant_C_list[[2]])
