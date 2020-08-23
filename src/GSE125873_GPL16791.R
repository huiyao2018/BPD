#load packages
library(tidyverse)
library(NetBID2)
# Define main working directory and project name
gse <- 'GSE125873'
gpl <- 'GPL16791'
project_main_dir <- '/cluster/home/chenhy/project/BPD/data/' 
project_name <- sprintf('%s_%s',gse,gpl) 
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

############### Step 1: Load in gene expression datasets for network construction (exp-load) ###############

## Depending on how the files are loaded, some of the code may need to be changed
datMeta <- read.table("/cluster/home/chenhy/project/BPD/data/Autism_GSE125873_GPL16791/DATA/GSE125873/GSE125873_series_matrix.txt",sep="\t",quote = "\"",header = T)
countData <- read.table("/cluster/home/chenhy/project/neuron_drivers_landscape/src/Autism_GSE64018_GPL11154/GSE64018_countlevel_12asd_12ctl.txt",sep="\t",header=TRUE)
colnames(countData) <- rownames(datMeta)
condition <- factor(datMeta$characteristics_ch1.2)
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design= ~ condition)
dds <- estimateSizeFactors(dds)
norcounts <- counts(dds, normalized=T)
rld <- rlog(dds, blind=FALSE)
boxplot(assay(rld))
datFPM <- assay(rld)
#datFPM <- sweep(log2(datExpr + 1), 2, log2(apply(datExpr,2,sum)/10^6))
keep <- which(rowMeans(datFPM) > 0)## Keep those with > 0 FPM in all of samples
datExpr <- datFPM[keep,]

## ref method:https://www.cell.com/cell/fulltext/S0092-8674(14)01512-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867414015128%3Fshowall%3Dtrue#secsectitle0030
# covariate vars:Differential expression was assessed using a linear model with log2(FPKM) as the outcome, and diagnosis, age, sex, RIN, brain bank, a surrogate variable for sequencing depth, and a surrogate variable for sequencing 3′ bias and GC content as covariates. 
#colnames(datExpr) <- rownames(datMeta)
## Do a PCA of the sequencing statistics of the full sample
colnames(datMeta) <- gsub(":ch1","",colnames(datMeta))
datSeq <- datMeta[,c(62,63,66,68:76,79,84,85)] ## All 15 columns, "TotalReads.picard" to "PropExonicReads.HTSC"
for (i in 1:ncol(datSeq)) {
  datSeq[,i] <- as.numeric(datSeq[,i])
}
datSeqNorm <- t(scale(datSeq,scale=F))
PC.datSeq <- prcomp(datSeqNorm);
varexp <- (PC.datSeq$sdev)^2 / sum(PC.datSeq$sdev^2)
print(varexp[1:2])
topPC.datSeq <- PC.datSeq$rotation[,1:2]; ## Explains 99% of variance in datSeq
colnames(topPC.datSeq) <- c("SeqPC1 - Depth","SeqPC2 - GC/Length")

## Get the data
condition <- 2-as.numeric(as.factor(datMeta[,"diagnosis"]))
age <- as.numeric(datMeta[,"age"])
sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
region <- as.numeric(as.factor(datMeta[,"regionid"]))-1
RIN <- as.numeric(datMeta[,"rin"])
bank <- as.numeric(as.factor(datMeta[,"brainbank"]))-1
seqStatPC1 <- topPC.datSeq[,1]
seqStatPC2 <- topPC.datSeq[,2]

varnames <- c("condition","age","sex","RIN","bank","seqStatPC1","seqStatPC2")
Bmat <- SEmat <- Pmat <- matrix(NA,nrow=nrow(datExpr),ncol=7)
colnames(Bmat) <- paste("beta",varnames,sep=".")
colnames(SEmat) <- paste("SE",varnames,sep=".")
colnames(Pmat) <- paste("p",varnames,sep=".")
rownames(Bmat) <- rownames(SEmat) <- rownames(Pmat) <- rownames(datExpr)

## Adjusted values
datExpr.reg <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)
regvars <- data.frame(condition=condition,age=age,sex=sex,RIN=RIN,bank=bank,seqStatPC1=seqStatPC1,seqStatPC2=seqStatPC2)
coefmat <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(regvars)+1)

for (i in 1:nrow(datExpr)) {
  if (i %% 1000 == 0) {cat(paste("On gene ",i,"\n",sep=""))}
  thisExpr <- as.numeric(datExpr[i,])
  lm.out <- lm(thisExpr ~ condition + age + sex + RIN + bank + seqStatPC1 + seqStatPC2)
  
  coef <- coef(lm.out) ## Get the coefficients from the model
  coefmat[i,] <- coef
  datExpr.reg[i,] <- coef[1] + coef[2]*regvars[,"condition"] + lm.out$residuals 
  
  tabOut <- summary(lm.out)$coefficients
  # Bmat[i,] <- tabOut[-c(1),"Estimate"]
  # SEmat[i,] <- tabOut[-c(1),"Std. Error"]
  # Pmat[i,] <- tabOut[-c(1),"Pr(>|t|)"]
}
gse64018 <- datExpr.reg
#gse64018 <- read.table(str_c(project_main_dir,project_name,"/DATA/GSE64018_adjfpkm_12asd_12ctl.txt"),header = T,col.names = rownames(pData(net_eset)),stringsAsFactors = F)
feature_df <- as.data.frame(rownames(gse64018))
colnames(feature_df) <- "ID"
feature_df$ID <- as.character(feature_df$ID)
rownames(feature_df) <- feature_df[,1]
feature_df <- feature_df[!duplicated(feature_df$ID),]
net_eset <- generate.eset(exp_mat = gse64018,
                          phenotype_info = pData(net_eset),
                          feature_info = feature_df,
                          annotation_info = "GPL11154")
# Remove duplicated data
# ID conversion, or merge transcript level to expression level, use_feature_info can be other dataframe instead of fData; optional;
# net_eset@featureData@data <- bitr(net_eset@featureData@data$gene, fromType="ENSEMBL", toType=c("SYMBOL", "GENENAME"), OrgDb="org.Hs.eg.db")
db.preload(use_level = 'gene')
transfer_tab <- read.csv(file = "/cluster/home/chenhy/project/neuron_drivers_landscape/result/transfer_tab.csv") %>% .[which(.$gene_biotype=="protein_coding"),1:2]
colnames(transfer_tab) <- c("ENSEMBL","SYMBOL")
net_eset@featureData@data <- transfer_tab
net_eset <- update_eset.feature(use_eset=net_eset,use_feature_info=fData(net_eset),from_feature='ENSEMBL',to_feature='SYMBOL',merge_method='median')

# Select phenotype columns or user added phenotype info; optional
phe <- pData(net_eset)
phe <- within(phe,{
  `diagnosis:ch1` <- NA
  `diagnosis:ch1`[characteristics_ch1.2 == "diagnosis: ASD"] <- "ASD"
  `diagnosis:ch1`[characteristics_ch1.2 == "diagnosis: CTL"] <- "TD"
})
net_eset <- update_eset.phenotype(use_eset=net_eset,use_phenotype_info=phe,use_sample_col='geo_accession',use_col='GEO-auto')
# before_QC
intgroups <- c("diagnosis","brainbank","Sex")
draw.eset.QC(net_eset,outdir=network.par$out.dir.QC,intgroup=intgroups,do.logtransform=FALSE,prefix='beforeQC_')

network.par$net.eset <- net_eset
# Save Step 1 network.par as RData
NetBID.saveRData(network.par = network.par,step='exp-load')

############### Step 2: Normalization for the expression dataset (exp-QC) ###############

# Reload network.par RData from Step 1
NetBID.loadRData(network.par = network.par,step='exp-load')
mat <- exprs(network.par$net.eset)
#查看行和列的缺失值
sample_na_count <- apply(mat,1,function(x){length(which(is.na(x)==TRUE))})
print(table(sample_na_count))
gene_na_count <- apply(mat,2,function(x){length(which(is.na(x)==TRUE))})
print(table(gene_na_count))
# Perform imputation
if(sum(sample_na_count)+sum(gene_na_count)>0) mat <- impute.knn(mat)$data

## Fourthly, filter out genes with very low expression values (bottom 5%) in most samples (more than 90%).
# Filter out low-expression genes
choose1 <- apply(mat<= quantile(mat, probs = 0.05), 1, sum)<= ncol(mat) * 0.90
print(table(choose1))
mat <- mat[choose1,]
phe <- pData(net_eset)
# Update eset
net_eset <- generate.eset(exp_mat=mat, phenotype_info=phe,
                          feature_info=fData(network.par$net.eset)[rownames(mat),],
                          annotation_info=annotation(network.par$net.eset))
draw.eset.QC(net_eset,outdir=network.par$out.dir.QC,intgroup=intgroups,do.logtransform=FALSE,prefix='afterQC_')
network.par$net.eset <- net_eset
NetBID.saveRData(network.par = network.par,step='exp-QC')

#===========================================analysis=================================================
######### Step 1: Load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###############
# Get the demo's constructed network data
NetBID.loadRData(network.par = network.par,step='exp-QC')
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
draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=intgroups,do.logtransform=FALSE,prefix='AC_')

# Save Step 2 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-get')

############### Step 3: Get differential expression (DE) / differential activity (DA) for drivers (act-DA) ###############

# Reload network.par RData from Step 2
analysis.par <- list()
analysis.par$out.dir.DATA <- '/cluster/home/chenhy/project/neuron_drivers_landscape/data/Autism_GSE64018_GPL11154//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='act-get')

# Create empty list to store comparison result
analysis.par$DE <- list()
analysis.par$DA <- list()

# First comparison: ASD vs. TD
comp_name <- 'ASD.Vs.TD' # Each comparison must has a name
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`diagnosis`=='ASD')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`diagnosis`=='TD')] # Control group
DE_gene_bid <- getDE.limma.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='ASD',G0_name='TD')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='ASD',G0_name='TD')
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



