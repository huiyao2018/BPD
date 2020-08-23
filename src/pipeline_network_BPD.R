library(NetBID2)
db.preload()
gse <- 'GSE106910'
gpl <- 'GPL16791'
project_main_dir <- '/cluster/home/chenhy/project/BPD/data/' 
project_name <- sprintf('%s_%s',gse,gpl) 
network.par <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name = project_name)
##### Step1: load in gene expression datasets for network construction
tmp <- load.exp.GEO(network.par$out.dir.DATA,GSE='GSE106910',GPL='GPL16791',getGPL = F)
phe <- pData(tmp)
mat <- read.delim(sprintf('%s/GSE106910_combined_counts_geo.txt',network.par$out.dir.DATA),stringsAsFactors = FALSE,row.names=1)
phe <- phe[,c('title','bronchopulmonary dysplasia severity:ch1','cga at collection:ch1','cell type:ch1','passages in-vitro:ch1','tissue:ch1')]
colnames(phe) <- c('title','BPD','cga','cell_type','passages','tissue')
rownames(phe) <- gsub(" ",".",phe[,1])
phe <- phe[colnames(mat),]
dds <- DESeqDataSetFromMatrix(mat,phe,design=~BPD)
dds <- DESeq(dds)
vsd <- vst(dds)
mat <- assay(vsd) ##
##
choose1 <- apply(mat<= quantile(mat, probs = 0.05), 1, sum)<= ncol(mat) * 0.90
print(table(choose1))
mat <- mat[choose1,]
draw.emb.kmeans(mat,obs_label = get_obs_label(phe,'BPD'),plot_type = '2D.interactive')
##
rm_sample <- 'S1.rep2'
phe <- phe[setdiff(rownames(phe),'S1.rep2'),]
mat <- mat[,rownames(phe)]
##
eset <- generate.eset(mat,phe)
network.par$net.eset <- eset
NetBID.saveRData(network.par = network.par,step='exp-load')

##########
pData(eset)$group <- ifelse(pData(eset)$BPD=='No BPD','No BPD','BPD')
phe_info <- pData(eset)
#save(eset,file='BPD.RData')
#load("/cluster/home/chenhy/project/BPD/data/GSE106910_GPL16791/DATA/BPD.RData")
# Update eset
net_eset <- generate.eset(exp_mat=mat, phenotype_info=phe_info,
                          feature_info=fData(eset)[rownames(mat),], 
                          annotation_info=annotation(eset))
use_genes <- rownames(mat)
transfer_tab <- get_IDtransfer(from_type = 'ensembl_gene_id',
                               to_type='external_gene_name',
                               use_genes=use_genes,
                               dataset='hsapiens_gene_ensembl')
transfer_tab <- dplyr::rename(transfer_tab,gene = ensembl_gene_id)
net_eset@featureData@data <- merge(net_eset@featureData@data,transfer_tab,by='gene')
net_eset <- update_eset.feature(use_eset=net_eset,use_feature_info=fData(net_eset),from_feature='gene',to_feature='external_gene_name',merge_method='median')

network.par$net.eset <- net_eset
NetBID.saveRData(network.par = network.par,step='exp-QC')
# Load database
db.preload(use_level='gene',use_spe='human',update=FALSE)

# Converts gene ID into the corresponding TF/SIG list
use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)

# Select samples for analysis
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) # here is using all samples, users can modify
prj.name <- network.par$project.name # if use different samples, need to change the project name
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0.5,IQR.loose_thre = 0.1,
                 SJAR.project_name=prj.name,SJAR.main_dir=network.par$out.dir.SJAR)
if(F){
  phe <- pData(eset)
  draw.eset.QC(eset,outdir = network.par$out.dir.QC)
  ################
  denovo_gene <- read.xlsx('../result.xlsx',sheet=2)
  g1 <- unique(denovo_gene[which(denovo_gene$Effect=='MODERATE'),'Genes'])
  ##
  load('BPD_BPD_net_2019-10-11_2019-10-11_16-27/BPD_BPD_net_2019-10-11/DATA/BPD_BPD_net.RData')
  net1 <- analysis.par$merge.network$igraph_obj
  transfer_tab <- analysis.par$transfer_tab
  ms_tab <- analysis.par$final_ms_tab
  vn <- V(net1)$name
  vn <- get_name_transfertab(vn,transfer_tab,from_type = 'ensembl_gene_id',to_type='external_gene_name')
  vn1 <- ms_tab[vn,'gene_label']; vn1[which(is.na(vn1)==TRUE)] <- vn[which(is.na(vn1)==TRUE)]
  V(net1)$name <- vn1
  ## sp
  g2 <- intersect(vn1,c(g1,paste0(g1,'_TF'),paste0(g1,'_SIG')))
  sp1 <- all_shortest_paths(net1,from=c('HIF1AN_SIG','HIF1AN'),to=g2,mode='all')
  ##
  sig_d <- draw.volcanoPlot(ms_tab,logFC_thre = 0.2,Pv_thre = 1e-5,
                            label_col='gene_label',logFC_col='logFC.BPD Vs. No BPD_DA',
                            Pv_col='P.Value.BPD Vs. No BPD_DA')
  sig_d_up <- sig_d[which(sig_d$`logFC.BPD Vs. No BPD_DA`>0),]
  sig_d_down <- sig_d[which(sig_d$`logFC.BPD Vs. No BPD_DA`<0),]
  ##
  g3 <- intersect(c(sig_d_up$gene_label,gsub('(.*)_.*','\\1',sig_d_up$gene_label)),vn1)
  g4 <- intersect(c(sig_d_down$gene_label,gsub('(.*)_.*','\\1',sig_d_down$gene_label)),vn1)
  ##
  sp1 <- lapply(g3,function(x)shortest_paths(net1,from=x,to=g2,mode='all'))
  names(sp1) <- g3
  gg <- c(g2,g3)
  view(as.matrix(gg[!grepl('_',gg)]))
  ##
  gg <- c(g2,g4)
  view(as.matrix(gg[!grepl('_',gg)]))
  
  ##
  gs.preload()
  res1 <- funcEnrich.Fisher(g3,use_gs = c('H','BP','CP:BIOCARTA','CP:REACTOME','CP:KEGG'),Pv_adj = 'none');
  draw.funcEnrich.cluster(res1,pdf_file = 'up.pdf',top_number = 50)
  res2 <- funcEnrich.Fisher(g4,use_gs = c('H','BP','CP:BIOCARTA','CP:REACTOME','CP:KEGG'));
  draw.funcEnrich.cluster(res2,pdf_file = 'down.pdf',top_number = 50)
}

##
############### Step 3: Prepare files to run SJARACNe (sjaracne-prep) ###############
NetBID.loadRData(network.par = network.par,step='exp-QC')
######### Step 1: Load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###############
# Get the demo's constructed network data
network.dir <- "/cluster/home/chenhy/project/BPD/data/GSE106910_GPL16791/"
network.project.name <- 'BPD'
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)
analysis.par$cal.eset <- network.par$net.eset

# Save Step 1 network.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')

############### Step 2: Read in network files and calcualte driver activity (act-get) ###############

# Reload network.par RData from Step 1
NetBID.loadRData(analysis.par=analysis.par,step='exp-QC')
# Get network information
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)

# Merge network first
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)

# Get activity matrix
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')

# Create eset using activity matrix
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

# QC plot for activity eset
intgroups <- "group"
draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=intgroups,do.logtransform=FALSE,prefix='AC_')

# Save Step 2 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-get')

############### Step 3: Get differential expression (DE) / differential activity (DA) for drivers (act-DA) ###############

# Reload network.par RData from Step 2
analysis.par <- list()
analysis.par$out.dir.DATA <- '/cluster/home/chenhy/project/BPD/data/GSE106910_GPL16791//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='act-get')

# Create empty list to store comparison result
analysis.par$DE <- list()
analysis.par$DA <- list()

# First comparison: ASD vs. TD
comp_name <- 'BPD.Vs.noBPD' # Each comparison must has a name
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$group=='BPD')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$group=='No BPD')] # Control group
DE_gene_bid <- getDE.limma.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='BPD',G0_name='No BPD')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='BPD',G0_name='No BPD')
# Save comparison result to list element in analysis.par, with comparison name
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid
# Save Step 3 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-DA')

############### Step 4: Generate a master table for drivers (ms-tab) ###############

# Reload analysis.par RData from Step 3
NetBID.loadRData(analysis.par=analysis.par,step='act-DA')

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(analysis.par$DE) # Users can use index or name to get target ones
# Prepare the conversion table (OPTIONAL)
use_genes <- unique(c(analysis.par$merge.network$network_dat$source.symbol,analysis.par$merge.network$network_dat$target.symbol))
transfer_tab <- get_IDtransfer2symbol2type(from_type = 'external_gene_name',use_genes=use_genes)
analysis.par$transfer_tab <- transfer_tab
# Creat the final master table
analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp,DE=analysis.par$DE,DA=analysis.par$DA,
                                                  target_list=analysis.par$merge.network$target_list,
                                                  tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                                  main_id_type='external_gene_name')

# Path and file name of the output EXCEL file
out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)
# Save the final master table as EXCEL file
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file)
# Save Step 4 analysis.par as RData, ESSENTIAL
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')

all_driver <- draw.volcanoPlot(dat=analysis.par$final_ms_tab,label_col='originalID_label',
                                         logFC_col='logFC.BPD.Vs.noBPD_DA',
                                         Pv_col='P.Value.BPD.Vs.noBPD_DA',
                                         logFC_thre=0.1,
                                         Pv_thre=1e-3,
                                         main='Volcano Plot for BPD.Vs.noBPD',
                                         show_label=F,
                                         label_type = 'origin',
                                         label_cex = 0.5,
                                         show_plot = T);nrow(all_driver)/nrow(analysis.par$final_ms_tab)
all_gene <- draw.volcanoPlot(dat=analysis.par$DE$BPD.Vs.noBPD,label_col='ID',
                                       logFC_col='logFC',
                                       Pv_col='adj.P.Val',
                                       logFC_thre=0.3,
                                       Pv_thre=1e-2,
                                       main='Volcano Plot for BPD.Vs.noBPD',
                                       show_label=F,
                                       label_type = 'origin',
                                       label_cex = 0.5,
                                       show_plot = T);nrow(all_gene)/nrow(analysis.par$DE$BPD.Vs.noBPD)
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')






##ewas 差异基因富集分析
gs.preload()
ewas_res1 <- funcEnrich.Fisher(input_list= unique(gene_ttest$Dis_gene),bg_list= unique(case_file$Report_Gene),use_gs=c('H','CP:REACTOME','BP','CGP','CP:KEGG'),
                          Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
pdf("/cluster/home/chenhy/project/BPD/result/plot/funcEnrich_ewas.pdf",width = 15,height = 6)
draw.funcEnrich.cluster(funcEnrich_res= ewas_res1,top_number=20,gs_cex = 1.2,gene_cex=1.2,pv_cex=1,Pv_thre=0.1,
                        cluster_gs=TRUE,cluster_gene = TRUE,h=0.95)
dev.off()

library(stringr)
ms_tab <- openxlsx::read.xlsx("/cluster/home/chenhy/project/BPD/data/GSE106910_GPL16791/DATA/BPD_BPD_net_2019-10-11_ms_tab.xlsx",1)
gene_ttest <- read.csv("/cluster/home/chenhy/project/BPD/result/ewas/gene_ttest_paired_panel.csv",stringsAsFactors = F) %>% .[order(.$p),]
variant_test <- read.csv("/cluster/home/chenhy/project/BPD/result/ewas/variant_fisher_paired_panel.csv",stringsAsFactors = F) %>% .[order(.$p),]
Dis_gene <- unique(gene_ttest$Dis_gene)
Dis_variant <- unique(variant_test$Dis_variant)
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
                               logFC_col='logFC.BPD.Vs..No.BPD_DA',
                               Pv_col='P.Value.BPD.Vs..No.BPD_DA',
                               logFC_thre=0.12,
                               Pv_thre=1e-3,
                               main='Volcano Plot for BPD.Vs..No.BPD_DA',
                               show_label=FALSE,
                               label_type = 'origin',
                               label_cex = 0.5,
                               show_plot = T);nrow(sig_driver)/nrow(ms_tab)*100
candidate_driver <- intersect(str_split(sig_driver$gene_label,"_",simplify = T)[,1],Dis_gene)
out2excel(ms_tab[which(str_split(ms_tab$gene_label,"_",simplify = T)[,1] %in% candidate_driver),],out.xlsx = "/cluster/home/chenhy/project/BPD/result/candidate_driver_mstab.xlsx")


driver_list <- rownames(sig_driver)
use_driver <- driver_list[3]
exp_mat <- Biobase::exprs(analysis.par$cal.eset)
## expression,the rownames could match originalID
ac_mat <- Biobase::exprs(analysis.par$merge.ac.eset)
## activity,the rownames could match originalID_label
phe_info <- Biobase::pData(analysis.par$cal.eset)
use_obs_class <- get_obs_label(phe_info = phe_info,'subgroup')
draw.categoryValue(ac_val=ac_mat[use_driver,],
                   exp_val=exp_mat[ms_tab[use_driver,'originalID'],],
                   use_obs_class=use_obs_class,
                   class_order=c('WNT','SHH','G4'),
                   class_srt=30,
                   main_ac = ms_tab[use_driver,'gene_label'],
                   main_exp=ms_tab[use_driver,'geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G4'='green'))










