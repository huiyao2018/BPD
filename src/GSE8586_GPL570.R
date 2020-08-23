#load packages
library(tidyverse)
library(NetBID2)
# Define main working directory and project name
gse <- 'GSE8586'
gpl <- 'GPL570'
project_main_dir <- "/cluster/home/chenhy/project/BPD/data/"
project_name <- sprintf('%s_%s',gse,gpl) 
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

############### Step 1: Load in gene expression datasets for network construction (exp-load) ###############
net_eset <- load.exp.GEO(out.dir=network.par$out.dir.DATA,GSE=gse,GPL=gpl,getGPL=TRUE,update=FALSE)

# Remove duplicated data
net_eset@featureData@data <- net_eset@featureData@data[!duplicated(net_eset@featureData@data$ID),]
# ID conversion, or merge transcript level to expression level, use_feature_info can be other dataframe instead of fData; optional;
db.preload(use_level='transcript',use_spe='human',update=FALSE)
transfer_tab <- get_IDtransfer(from_type = 'affy_hg_u133_plus_2',to_type = 'external_gene_name')
transfer_tab <- dplyr::rename(transfer_tab,ID = affy_hg_u133_plus_2)
transfer_tab <- transfer_tab[!duplicated(transfer_tab$ID),] #删除一个探针对应多个基因名的探针
net_eset@featureData@data <- merge(net_eset@featureData@data,transfer_tab,by='ID')
net_eset <- update_eset.feature(use_eset=net_eset,use_feature_info=fData(net_eset),from_feature='ID',to_feature='external_gene_name',merge_method='median')

# Select phenotype columns or user added phenotype info
phe <- pData(net_eset)
phe <- within(phe,{
  `diagnosis:ch1` <- NA
  `diagnosis:ch1`[str_detect(title,"_bpd")] <- "case"
  `diagnosis:ch1`[str_detect(title,"_nobpd")] <- "control"
})
net_eset <- update_eset.phenotype(use_eset=net_eset,use_phenotype_info=phe,use_sample_col='geo_accession',use_col='GEO-auto')
network.par$net.eset <- net_eset
# Save Step 1 network.par as RData
NetBID.saveRData(network.par = network.par,step='exp-load')

############### Step 2: Normalization for the expression dataset (exp-QC) ###############

# Reload network.par RData from Step 1
NetBID.loadRData(network.par = network.par,step='exp-load')
# QC
intgroups <- get_int_group(network.par$net.eset)
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=intgroups,do.logtransform=FALSE,prefix='beforeQC_')

mat <- exprs(network.par$net.eset)
#查看行和列的缺失值
sample_na_count <- apply(mat,1,function(x){length(which(is.na(x)==TRUE))});print(table(sample_na_count==0))
gene_na_count <- apply(mat,2,function(x){length(which(is.na(x)==TRUE))});print(table(gene_na_count))
# Perform imputation
if(sum(sample_na_count)+sum(gene_na_count)>0) mat <- impute.knn(mat)$data
## Secondly, the log2 transformation.
med_val <- median(apply(mat,2,median)); print(med_val);if(med_val>16){mat <- log2(mat)}
## Thirdly, the quantile normalization across samples.
# Perform limma quantile normalization
mat <- normalizeQuantiles(mat)
## Fourthly, filter out genes with very low expression values (bottom 5%) in most samples (more than 90%).
# Filter out low-expression genes
choose1 <- apply(mat<= quantile(mat, probs = 0.05), 1, sum)<= ncol(mat) * 0.90; print(table(choose1));mat <- mat[choose1,]


# Update eset
net_eset <- generate.eset(exp_mat=mat, phenotype_info=phe,
                          feature_info=fData(network.par$net.eset)[rownames(mat),], 
                          annotation_info=annotation(network.par$net.eset))
# QC
draw.eset.QC(net_eset,outdir=network.par$out.dir.QC,intgroup=intgroups,do.logtransform=FALSE,prefix='afterQC_')

network.par$net.eset <- net_eset
NetBID.saveRData(network.par = network.par,step='exp-QC')

############### Step 3: Prepare files to run SJARACNe (sjaracne-prep) ###############
NetBID.loadRData(network.par = network.par,step='exp-QC')
######### Step 1: Load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###############
# Get the demo's constructed network data
network.dir <- "/cluster/home/chenhy/project/Brainspan/network/brainspan_nonorm"
network.project.name <- 'all_brain'
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
draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=intgroups,do.logtransform=FALSE,prefix='removebatch_AC_')

# Save Step 2 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-get')

############### Step 3: Get differential expression (DE) / differential activity (DA) for drivers (act-DA) ###############

# Reload network.par RData from Step 2
analysis.par <- list()
analysis.par$out.dir.DATA <- '/cluster/home/chenhy/project/neuron_drivers_landscape/data/Autism_GSE113834_GPL15207//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='act-get')

# Create empty list to store comparison result
analysis.par$DE <- list()
analysis.par$DA <- list()

# First comparison: ASD vs. TD
comp_name <- 'ASD.Vs.TD' # Each comparison must has a name
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$group=='ASD')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$group=='CTRL')] # Control group
DE_gene_bid <- getDE.limma.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='ASD',G0_name='CTRL')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='ASD',G0_name='CTRL')
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



