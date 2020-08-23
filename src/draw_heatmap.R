library(RColorBrewer)
library(gplots)
###
draw_heatmap <- function(draw_mat,gene_category,phenotype_category,row_text = FALSE,col_text = FALSE,
                         col1 = brewer.pal(8, "Set3"),col2 = brewer.pal(8, "Set2"),
                         col_cex=1,min_val=-2,max_val=2,legend_pos='bottomright',
                         legend_col=rev(brewer.pal(11,'RdBu'))){
  phe <- phenotype_category[colnames(draw_mat)]
  all_p <- unique(phe)
  col2 <- colorRampPalette(col2)(length(all_p)); names(col2) <- all_p
  if(is.null(gene_category)==F){
    gen <- gene_category[rownames(draw_mat)] ## gene category
    all_c <- unique(gen)
    col1 <- colorRampPalette(col1)(length(all_c)); names(col1) <- all_c
    mod_data <- do.call(rbind,lapply(all_c,function(x){
      use_gene <- names(which(gen==x))
      use_data <- draw_mat[use_gene,]
      if(length(use_gene)==1) return(use_data)
      ord <- hclust(dist(use_data))$order ## order of the genes
      use_data[ord,]
    }))
    gen <- gen[rownames(mod_data)]
    draw_mat <- t(mod_data)
  }else{
    draw_mat <- t(draw_mat)
  }
  ##
  bb <- seq(min_val,max_val,length.out = 101)
  cc <- colorRampPalette(legend_col)(length(bb)-1)
  #par(mar=c(6,6,4,8))
  image(draw_mat,col=cc,breaks=bb,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',xaxs='i',yaxs='i')
  pp <- par()$usr ## get to know the xy position of the figure
  x_line_pos <- seq(pp[1],pp[2],length.out=nrow(draw_mat)+1) ## x position for vertical line
  y_line_pos <- seq(pp[3],pp[4],length.out=ncol(draw_mat)+1) ## y position for each horizon line
  x_mid_pos  <- x_line_pos[1:(length(x_line_pos)-1)]/2+x_line_pos[2:length(x_line_pos)]/2 ## x position for column mid position
  y_mid_pos  <- y_line_pos[1:(length(y_line_pos)-1)]/2+y_line_pos[2:length(y_line_pos)]/2 ## y position for row mid position
  # add text
  if(row_text==TRUE) text(x=max(x_line_pos),y=y_mid_pos,colnames(draw_mat),xpd=T,pos=4,cex=0.65)
  if(col_text==TRUE) text(y=min(y_line_pos),x=x_mid_pos,rownames(draw_mat),xpd=T,cex=1.2,pos=1)
  # try to draw line to seperate groups
  # try to draw box to seperate groups
  if(is.null(gene_category)==FALSE){
    w1 <- cumsum(table(gen)[unique(gen)])
    segments(x0=min(x_line_pos),x1=max(x_line_pos),y0=y_line_pos[c(1,w1+1)],y1=y_line_pos[c(1,w1+1)])
    dw <- (pp[2]-pp[1])/20
    rect(xleft=pp[1]-dw,xright=pp[1]-dw*0.1,ybottom=y_line_pos[1:(length(y_line_pos)-1)],
         ytop=y_line_pos[2:(length(y_line_pos))],xpd=T,col=col1[gen],border=NA)
    text(x=pp[1]-2*dw,y=y_line_pos[1+w1]/2+y_line_pos[1+c(0,w1[1:(length(w1)-1)])]/2,all_c,srt=90,xpd=T,adj=0.5)
  }
  #
  if(is.null(phenotype_category)==FALSE){
    w2 <- cumsum(table(phe)[unique(phe)])
    segments(y0=min(y_line_pos),y1=max(y_line_pos),x0=x_line_pos[c(1,w2+1)],x1=x_line_pos[c(1,w2+1)])
    dw <- (pp[4]-pp[3])/30
    rect(ybottom=pp[4]+0.1*dw,ytop=pp[4]+dw,xleft=x_line_pos[1:(length(x_line_pos)-1)],
         xright=x_line_pos[2:(length(x_line_pos))],xpd=T,col=col2[phe],border=NA)
     text(y=pp[4]+1.5*dw,x=x_line_pos[1+w2]/2+x_line_pos[1+c(0,w2[1:length(w2)-1])]/2,unique(phe),srt=0,xpd=T,adj=0.5,cex=col_cex)
  }
  # legend
  if(legend_pos=='bottomright'){
    bb <- seq(min_val,max_val,length.out = 11)
    cc <- colorRampPalette(legend_col)(length(bb)-1)
    dw <- (pp[2]-pp[1])/4
    dh <- (pp[4]-pp[3])/20
    ss1 <- seq(max(x_line_pos),max(x_line_pos)+dw,length.out = length(bb))
    ss1_left <- ss1[1:(length(ss1)-1)]
    ss1_right <- ss1[2:length(ss1)]
    rect(xleft=ss1_left,xright=ss1_right,ybottom = min(y_line_pos)-dh,ytop=min(y_line_pos)-dh*0.6,col=cc,xpd=T,border=NA)
    text(x=ss1,y=min(y_line_pos)-1.1*dh,bb,xpd=T,srt=90,cex=0.5,adj=1)
    text(x=mean(ss1),y=min(y_line_pos)-dh*0.5,'Value',cex=0.6,xpd=T)
  }
  if(legend_pos=='topleft'){
    bb <- seq(min_val,max_val,length.out = 11)
    cc <- colorRampPalette(legend_col)(length(bb)-1)
    dw <- (pp[2]-pp[1])/3
    dh <- (pp[4]-pp[3])/10
    ss1 <- seq(min(x_line_pos),min(x_line_pos)+dw,length.out = length(bb))
    ss1_left <- ss1[1:(length(ss1)-1)]
    ss1_right <- ss1[2:length(ss1)]
    rect(xleft=ss1_left,xright=ss1_right,ybottom = max(y_line_pos)+1.2*dh,ytop=max(y_line_pos)+dh*1.6,col=cc,xpd=T,border=NA)
    text(x=ss1,y=max(y_line_pos)+1.1*dh,bb,xpd=T,srt=90,cex=0.5,adj=1,xpd=T)
    text(x=mean(ss1),y=max(y_line_pos)+dh*1.8,'Value',cex=1,xpd=T,xpd=T)
  }
}
##
list2mat <- function(input_list){
  all_x <- base::unique(unlist(input_list))
  all_y <- base::unique(names(input_list))
  mat1 <- matrix(0,nrow=base::length(all_x),ncol=base::length(all_y))
  rownames(mat1) <- all_x; colnames(mat1) <- all_y;
  for(i in names(input_list)){
    mat1[input_list[[i]],i] <- 1
  }
  return(mat1)
}
vec2list <- function(input_v,sep=NULL){
  if(is.null(sep)==TRUE){
    tmp2 <- list()
    input_vn <- names(input_v)
    input_v <- clean_charVector(input_v); names(input_v) <- input_vn
    for(i in 1:base::length(input_v)){
      if(input_v[i] %in% names(tmp2)){
        tmp2[[input_v[i]]] <- c(tmp2[[input_v[i]]],names(input_v)[i])
      }else{
        tmp2[[input_v[i]]] <- names(input_v)[i]
      }
    }
  }else{
    tmp1 <- stats::aggregate(names(input_v),list(input_v),function(x)base::paste(x,collapse=sep))
    tmp2 <- tmp1$x; names(tmp2) <- tmp1$Group.1
  }
  tmp2
}
clean_charVector <- function(x){
  x1 <- names(x)
  x <- as.character(x);
  x[which(x=='')] <- 'NULL';
  x[which(is.null(x)==TRUE)] <- 'NULL'
  x[which(is.na(x)==TRUE)] <- 'NA'
  names(x) <- x1
  x
}
strwidthMod <- function(s, units = "inch", cex = 1,ori=TRUE,mod=FALSE){
  if(ori==TRUE) return(strwidth(s=s,units=units,cex=cex))
  if(mod==TRUE){
    plot.new()
    rt <- strwidth(s,units=units)/strwidth('W',units=units); rt <- ceiling(rt)
    if(units=='user') r1 <- par.char2pos()[1]*cex*rt
    if(units=='inch') r1 <- par.char2inch()[1]*cex*rt
    dev.off(); return(r1)
  }else{
    if(units=='user') return(par.char2pos()[1]*cex*nchar(s))
    if(units=='inch'| units=='inches') return(par.char2inch()[1]*cex*nchar(s))
  }
}
par.char2inch <- function(){return(par()$cin)} ## letter W
strheightMod <- function(s, units = "inch", cex = 1,ori=TRUE,mod=FALSE){
  if(ori==TRUE) return(strheight(s=s,units=units,cex=cex))
  if(units=='user') return(par.char2pos()[2]*cex)
  if(units=='inch' | units=='inches') return(par.char2inch()[2]*cex)
}
par.char2pos <- function(){par()$cxy}

draw.funcEnrich.cluster.mod <- function (funcEnrich_res = NULL, top_number = 30, Pv_col = "Ori_P", 
          name_col = "#Name", item_col = "Intersected_items", 
          Pv_thre = 0.1, gs_cex = 0.7, gene_cex = 0.8, pv_cex = 0.7, 
          main = "", h = 0.95, cluster_gs = TRUE, cluster_gene = TRUE, 
          pdf_file = NULL, use_genes = NULL, return_mat = FALSE,
          col_mid=(grDevices::colorRampPalette(brewer.pal(9, "Reds")[3:9]))(100),
          pv_col = brewer.pal(9, "Set1")[1]) 
{
#  all_input_para <- c("funcEnrich_res", "Pv_col", 
#                      "item_col", "Pv_thre", "name_col", 
#                      "gs_cex", "gene_cex", "pv_cex", "main", 
#                      "h", "cluster_gs", "cluster_gene", 
#                      "return_mat")
#  check_res <- sapply(all_input_para, function(x) check_para(x, 
#                                                             envir = environment()))
#  if (min(check_res) == 0) {
#    message("Please check and re-try!")
#    return(FALSE)
#  }
#  check_res <- c(check_option("cluster_gs", c(TRUE, FALSE), 
#                              envir = environment()), check_option("cluster_gene", 
#                                                                   c(TRUE, FALSE), envir = environment()), check_option("return_mat", 
#                                                                                                                        c(TRUE, FALSE), envir = environment()))
#  if (min(check_res) == 0) {
#    message("Please check and re-try!")
#    return(FALSE)
#  }
  if (is.null(top_number) == TRUE) 
    top_number <- nrow(funcEnrich_res)
  funcEnrich_res <- funcEnrich_res[which(funcEnrich_res[, Pv_col] <= 
                                           Pv_thre), ]
  if (nrow(funcEnrich_res) > top_number) 
    funcEnrich_res <- funcEnrich_res[1:top_number, ]
  pv_val <- funcEnrich_res[, Pv_col]
  names(pv_val) <- rownames(funcEnrich_res)
  all_g2s <- lapply(funcEnrich_res[, item_col], function(x1) unlist(strsplit(x1, 
                                                                             ";")))
  names(all_g2s) <- funcEnrich_res[, name_col]
  mat1 <- t(list2mat(all_g2s))
  mat1 <- mat1[rev(funcEnrich_res[, name_col]), ]
  if (is.null(use_genes) == FALSE) 
    mat1 <- mat1[, base::intersect(colnames(mat1), use_genes)]
  if (ncol(mat1) == 0) {
    message("No genes left, please check and re-try!")
    return(FALSE)
  }
  h_gs <- hclust(dist(mat1, method = "binary"))
  h_gene <- hclust(dist(t(mat1), method = "binary"))
  gs_cluster <- cutree(h_gs, h = h)
  gene_cluster <- cutree(h_gene, h = h)
  if (cluster_gs == FALSE) {
    gs_cluster <- rep(1, length.out = nrow(mat1))
    names(gs_cluster) <- rownames(mat1)
  }
  if (cluster_gene == FALSE) {
    gene_cluster <- rep(1, length.out = ncol(mat1))
    names(gene_cluster) <- colnames(mat1)
  }
  cc1 <- (grDevices::colorRampPalette(brewer.pal(8, "Dark2")))(base::length(base::unique(gs_cluster)))
  cc2 <- (grDevices::colorRampPalette(brewer.pal(9, "Pastel1")))(base::length(base::unique(gene_cluster)))
  #cc3 <- (grDevices::colorRampPalette(brewer.pal(9, "Reds")[3:9]))(100)
  cc3 <- col_mid
  if (cluster_gs == TRUE) 
    gs_cluster <- gs_cluster[h_gs$order]
  tmp2 <- vec2list(gs_cluster, sep = NULL)
  tmp2 <- tmp2[rev(order(unlist(lapply(tmp2, function(x) base::min(funcEnrich_res[x, 
                                                                                  Pv_col])))))]
  mat1 <- mat1[unlist(tmp2), ]
  if (cluster_gene == TRUE) 
    mat1 <- mat1[, h_gene$order]
  gs_cluster <- gs_cluster[rownames(mat1)]
  gene_cluster <- gene_cluster[colnames(mat1)]
  pv <- pv_val[rownames(mat1)]
  pv1 <- format(pv, scientific = TRUE, digits = 3)
  plot_part <- function(ori = TRUE, before_off = FALSE) {
    gsWidth <- base::max(strwidthMod(rownames(mat1), units = "inch", 
                                     cex = gs_cex, ori = ori)) + par.char2inch()[1]
    gsHeight <- base::max(strheightMod(rownames(mat1), units = "inch", 
                                       cex = gs_cex)) * nrow(mat1) * 1.75
    geneWidth <- base::max(strheightMod(colnames(mat1), units = "inch", 
                                        cex = gene_cex)) * ncol(mat1) * 1.5
    geneHeight <- base::max(strwidthMod(colnames(mat1), units = "inch", 
                                        cex = gene_cex, ori = ori)) * 1.05 + par.char2inch()[2] * 
      1.1
    pvWidth <- base::max(strwidthMod(pv1, units = "inch", 
                                     cex = pv_cex, ori = ori)) + par.char2inch()[1]
    pvHeight <- base::max(strheightMod(pv1, units = "inch", 
                                       cex = pv_cex, ori = ori)) * nrow(mat1) * 1.75
    gsHeight <- base::max(gsHeight, pvHeight)
    mr <- 1/pvWidth
    geneWidth1 <- ceiling((geneWidth + 0.5) * mr)
    pvWidth1 <- 1
    gsWidth1 <- ceiling((gsWidth + 0.5) * mr)
    if (geneWidth1 + pvWidth1 + gsWidth1 > 200) {
      mr <- 180/(geneWidth1 + pvWidth1 + gsWidth1)
      geneWidth1 <- round(geneWidth1 * mr)
      pvWidth1 <- round(pvWidth1 * mr)
      gsWidth1 <- ceiling(gsWidth1 * mr)
    }
    ww <- (gsWidth + pvWidth + geneWidth) + 1
    hh <- geneHeight + gsHeight + 0.5
    if (before_off == TRUE) 
      dev.off()
    if (is.null(pdf_file) == FALSE) 
      pdf(file = pdf_file, width = ww, height = hh)
    graphics::layout(t(matrix(c(rep(1, geneWidth1), rep(2, 
                                                        pvWidth1), rep(3, gsWidth1)), byrow = TRUE)))
    par(mai = c(0.5, 0.5, geneHeight, 0))
    graphics::image(t(mat1), col = c("white", cc3[1]), 
                    xaxt = "n", yaxt = "n", bty = "n")
    pp <- par()$usr
    gs_cs <- cumsum(base::table(gs_cluster)[base::unique(gs_cluster)])
    gene_cs <- cumsum(base::table(gene_cluster)[base::unique(gene_cluster)])
    xx <- (pp[2] - pp[1])/base::length(gene_cluster)
    yy <- (pp[4] - pp[3])/base::length(gs_cluster)
    graphics::abline(h = gs_cs * yy + pp[3], col = "black", 
                     lwd = 0.25)
    graphics::abline(v = gene_cs * xx + pp[1], col = "black", 
                     lwd = 0.25)
    graphics::abline(v = pp[1:2], col = "black", lwd = 0.5)
    graphics::abline(h = pp[3:4], col = "black", lwd = 0.5)
    yy <- par.char2pos()[2] * 0.9
    if (cluster_gene == TRUE) {
      graphics::text(c(1:base::length(gene_cluster)) * 
                       xx + pp[1] - xx/2, pp[4] + yy, colnames(mat1), 
                     xpd = TRUE, adj = 0, cex = gene_cex, srt = 90)
      graphics::rect(xleft = c(1:base::length(gene_cluster)) * 
                       xx + pp[1] - xx, xright = c(1:base::length(gene_cluster)) * 
                       xx + pp[1], ybottom = pp[4], ytop = pp[4] + yy * 
                       0.8, col = cc2[gene_cluster[colnames(mat1)]], 
                     xpd = TRUE, border = NA)
    }
    else {
      graphics::text(c(1:base::length(gene_cluster)) * 
                       xx + pp[1] - xx/2, pp[4] + 0.5 * yy, colnames(mat1), 
                     xpd = TRUE, adj = 0, cex = gene_cex, srt = 90)
    }
    par(mai = c(0.5, 0, geneHeight, 0))
    graphics::plot(1, xaxt = "n", yaxt = "n", 
                   bty = "n", xlim = c(pp[1], pp[2]), ylim = c(pp[3], 
                                                               pp[4]), col = "white", xlab = "", 
                   ylab = "")
    pp <- par()$usr
    yy <- (pp[4] - pp[3])/base::length(gs_cluster)

    pv_c <- z2col(-qnorm(pv),red_col = pv_col)
    graphics::rect(xleft = pp[1], xright = pp[2], ybottom = c(1:base::length(gs_cluster)) * 
                     yy + pp[3] - yy, ytop = c(1:base::length(gs_cluster)) * 
                     yy + pp[3], col = pv_c, border = NA)
    graphics::text(0.5, c(1:base::length(gs_cluster)) * yy + 
                     pp[3] - yy/2, pv1, xpd = TRUE, adj = 0.5, cex = pv_cex)
    graphics::abline(v = pp[1], col = "black", lwd = 2)
    par(mai = c(0.5, 0, geneHeight, 0.5))
    zz <- base::min(c(xx, yy))
    graphics::plot(1, xaxt = "n", yaxt = "n", 
                   bty = "n", xlim = c(pp[1], pp[2]), ylim = c(pp[3], 
                                                               pp[4]), col = "white", xlab = "", 
                   ylab = "")
    pp <- par()$usr
    yy <- (pp[4] - pp[3])/base::length(gs_cluster)
    graphics::text(pp[1] + zz * 0.2, c(1:base::length(gs_cluster)) * 
                     yy + pp[3] - yy/2, rownames(mat1), xpd = TRUE, adj = 0, 
                   cex = gs_cex)
    graphics::abline(h = pp[3:4], col = "black", lwd = 0.5)
    graphics::abline(v = pp[1], col = "black", lwd = 0.5)
    graphics::abline(h = gs_cs * yy + pp[3], col = "black", 
                     lwd = 0.25)
  }
  if (is.null(pdf_file) == FALSE) {
    plot_part(ori = TRUE)
    plot_part(ori = TRUE, before_off = TRUE)
    while (!is.null(dev.list())) dev.off()
  }
  else {
    plot_part()
  }
  if (return_mat == TRUE) {
    return(mat1)
  }
  else {
    return(TRUE)
  }
}

