# library(VennDiagram)
# library(grid)
# library(futile.logger)


venn.plot_up <- venn.diagram(
  #数据列表
  x = list(
    twoSon_VariantGene = select_variant$gene,
    GSE106910_driver = str_split(rownames(all_driver),"_",simplify = T)[,1],
    paired_test_gene = paired_test_gene
  ),
  filename = NULL,    #保存路径
  #main="Overlap of up regulated genes",
  #sub = "Subtitle",
  height = 450, 
  width = 450,
  resolution =300, 
  #imagetype="png", 
  col = "transparent",      #指定图形的圆周边缘颜色  transparent 透明           
  fill = c(colors()[148], colors()[589], colors()[116]),  #填充颜色
  alpha = c(0.6, 0.6, 0.6),                                      #透明度
  label.col = c("orange", "white", "darkorchid4", "white",
                "white", "darkgreen", "white"),
  lwd = 0.5,
  cex = 1.5,    #每个区域label名称的大小
  fontfamily = "serif",  #字体
  fontface = "bold",     #字体格式
  cat.col = c("black"),  #分类颜色 
  cat.cex = 1,      #每个分类名称大小
  cat.pos = c(100, 260, 0),        #
  cat.dist = c(0.17, 0.17, 0.05),    #
  cat.fontfamily = "serif",     #分类字体
  rotation.degree =180,        #旋转角度
  margin = 0.2               #在网格单元中给出图周围空白量的编号
);
#可以不保存查看图片，但是效果不佳（命令如下，但是需要首先把filename设置为（filename=NULL））
grid.draw(venn.plot_up)
#dev.off()




