#####################################
# Dimensionality reduction techniques smart-seq2
#
#
# Mouse E5.0 analysis
#####################################
#set working directory
path_data="C:/"
setwd(path_data)

set.seed(1234567) 
library(corrplot)
library(RColorBrewer)
library(ggplot2)
library(mclust)
library(pheatmap)
library("scatterplot3d")
library(umap)
library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
source("http://bioconductor.org/biocLite.R")
library("KEGGREST")
library(circlize)

#Load in featurecounts file and QC and trim samples
ALL_DATA_plate1 <- read.table('C:/featurecountsplate1.txt', head=TRUE)
sample_key<-read.csv('C:/samplekey.csv',stringsAsFactors = F)
qc_mouse<-read.table('C:/QC.txt', head=TRUE)
ALL_DATA_counts<-ALL_DATA_plate1[,c(7:ncol(ALL_DATA))]
gene_ids<-ALL_DATA_plate1[,1]
alldata_ids<-colnames(ALL_DATA_counts)
rownames(ALL_DATA_counts)<-gene_ids
pass_qc<-unlist(lapply(qc_mouse$File[qc_mouse$percent>.60],grep,alldata_ids))
ALL_DATA_trimed<-ALL_DATA_counts[,pass_qc]
colnames(ALL_DATA_trimed)<-colnames(ALL_DATA_counts[,pass_qc])
sample_key_trim<-sample_key[match(colnames(ALL_DATA_trimed),sample_key[,1]),]
colnames(ALL_DATA_trimed)<-sample_key_trim[,2]

ALL_DATA_2<-read.table('C:/featurecountsplate2.txt', head=TRUE)
mouse_2_qc <- read.table('C:/Users/Dylan/Documents/cambridge/forDylan/mouse_analysis/QC_antoniaplate23.txt', head=F)
ALL_DATA_2_genes<-ALL_DATA_2[,1]
ALL_DATA_2<-ALL_DATA_2[,-1]
mouse_2_qc<-mouse_2_qc[match(colnames(ALL_DATA_23),mouse_2_qc[,1]),]
samples_DATA_3<-ALL_DATA_2[,mouse_2_qc[,5]==1]

#Generate seurat objects and plot E5.0 clustering 

#plate 1
mouse_data <- CreateSeuratObject(counts = ALL_DATA_trimed, assay = "RNA",min.cells = 0, min.features = 0)
mouse_data$species <- "1) Mouse"
mouse_data$divergence1 <- "Mouse"
mouse_data <- subset(mouse_data, subset = nFeature_RNA > 0)
mouse_data <- NormalizeData(mouse_data, verbose = FALSE)
mouse_data <- FindVariableFeatures(mouse_data, selection.method = "vst", nfeatures = 20000)
mouse_data <- ScaleData(mouse_data, verbose = FALSE)
mouse_data <- RunPCA(mouse_data, npcs = 20, verbose = FALSE)
mouse_data <- RunUMAP(mouse_data, reduction = "pca", dims = 1:20)
mouse_data <- RunTSNE(mouse_data, reduction = "pca", dims = 1:20)
Idents(mouse_data)<-sample_key_trim[,3]
DimPlot(mouse_data, reduction = "umap",dims = c(1, 2),pt.size=10) 

#plate 2
mouse_data3 <- CreateSeuratObject(counts = samples_DATA_3, assay = "RNA",min.cells = 0, min.features = 0)
mouse_data3$species <- "1) Mouse"
mouse_data3$divergence1 <- "Mouse"
mouse_data3 <- subset(mouse_data3, subset = nFeature_RNA > 0)
mouse_data3 <- NormalizeData(mouse_data3, verbose = FALSE)
mouse_data3 <- FindVariableFeatures(mouse_data3, selection.method = "vst", nfeatures = 20000)
mouse_data3 <- ScaleData(mouse_data3, verbose = FALSE)

mouse_data3 <- RunPCA(mouse_data3, npcs = 20, verbose = FALSE)
mouse_data3 <- RunUMAP(mouse_data3, reduction = "pca", dims = 1:20)
mouse_data3 <- RunTSNE(mouse_data3, reduction = "pca", dims = 1:20)
Idents(mouse_data3)<-mouse23_samples_only[which(mouse23_samples_only[,4] == "WT2"|mouse23_samples_only[,4] == "DEL4"|mouse23_samples_only[,4] == "DEL5"|mouse23_samples_only[,4] == "DEL6"),4]
DimPlot(mouse_data3, reduction = "umap",dims = c(1, 2),pt.size=10) 

### CCA 
mammal.anchors <- FindIntegrationAnchors(object.list = list(mouse_data,mouse_data3), dims = 1:20, anchor.features = 20000,k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20, k.weight = 50)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined, reduction = "umap",dims = c(1, 2),pt.size=10) 

#separate Wild type and KO embryos as determined by read alignment
WT_samples=subset(x=mammal.combined,idents=c("WT1","WT2","DEL1"))
KO_samples=subset(x=mammal.combined,idents=c("DEL2","DEL3","DEL4","DEL5","DEL6"))

#perform unbiased clustering and label clusters by gene expression
mammal.combined <- FindNeighbors(mammal.combined, reduction = "umap", dims = 1:2)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
clusters<-Idents(mammal.combined)
clusters<-as.character(clusters)
clusters[clusters==0]<-rep("E5.0_KO_VE",length(clusters[clusters==0]))
clusters[clusters==1]<-rep("E5.0_TE",length(clusters[clusters==1]))
clusters[clusters==2]<-rep("E5.0_EPI",length(clusters[clusters==2]))
clusters[clusters==3]<-rep("E5.0_EPI",length(clusters[clusters==3]))
clusters[clusters==4]<-rep("E5.0_WT_VE",length(clusters[clusters==4]))
Idents(mammal.combined)<-clusters
mammal.combined_invitro<-mammal.combined 
DimPlot(mammal.combined, reduction = "pca",dims = c(1, 2), cols=cols ,split.by="species",  label = FALSE,label.size = 4,pt.size =3) 
DimPlot(mammal.combined, reduction = "umap",dims = c(1, 2), split.by="species",  label = FALSE,label.size = 5,pt.size =5) 
FeaturePlot(mammal.combined,features = c("Pou5f1","Otx2","Cdx2","Gata6"),pt.size=5)


ALL_DATA_mohammed <- read.table("GSE121650_rna_counts.tsv",header = T)
ALL_DATA_mohammed<-ALL_DATA_mohammed[!duplicated(ALL_DATA_mohammed[,1]),]
mohamed_genes<-ALL_DATA_mohammed$ensembl_id
ALL_DATA_mohammed<-ALL_DATA_mohammed[,-1]
rownames(ALL_DATA_mohammed)<-mohamed_genes
mohammed_key <- read.table("GSE121650_sample_metadata.txt",header = T)
mohammed_key_trim<-mohammed_key[match(colnames(ALL_DATA_mohammed),mohammed_key[,1]),]
ALL_DATA_mohammed<-ALL_DATA_mohammed[,mohammed_key_trim$pass_rnaQC]
mohammed_key_trim<-mohammed_key_trim[mohammed_key_trim$pass_rnaQC,]
ALL_DATA_mohammed_trim<-ALL_DATA_mohammed[,-which(is.na(mohammed_key_trim$lineage10x_2))]
mohammed_key_trim<-mohammed_key_trim[-which(is.na(mohammed_key_trim$lineage10x_2)),]
ALL_DATA_mohammed_trim<-ALL_DATA_mohammed_trim[,-grep("7.5",mohammed_key_trim$stage)]
mohammed_key_trim<-mohammed_key_trim[-grep("7.5",mohammed_key_trim$stage),]

#Generate seurat objects
mohammed_key_trim_fullID <- paste(mohammed_key_trim$stage, mohammed_key_trim$lineage10x_2, sep="_")

mouse_data_mohammed <- CreateSeuratObject(counts = ALL_DATA_mohammed_trim, assay = "RNA",min.cells = 0, min.features = 0)
mouse_data_mohammed$species <- "invivo"
mouse_data_mohammed$divergence1 <- "Mouse"
mouse_data_mohammed <- subset(mouse_data_mohammed, subset = nFeature_RNA > 0)
mouse_data_mohammed <- NormalizeData(mouse_data_mohammed, verbose = FALSE)
mouse_data_mohammed <- FindVariableFeatures(mouse_data_mohammed, selection.method = "vst", nfeatures = 4000)
mouse_data_mohammed <- ScaleData(mouse_data_mohammed, verbose = FALSE)

mouse_data_mohammed <- RunPCA(mouse_data_mohammed, npcs = 20, verbose = FALSE)
mouse_data_mohammed <- RunUMAP(mouse_data_mohammed, reduction = "pca", dims = 1:20)
mouse_data_mohammed <- RunTSNE(mouse_data_mohammed, reduction = "pca", dims = 1:20)
Idents(mouse_data_mohammed)<-mohammed_key_trim_fullID
DimPlot(mouse_data_mohammed, reduction = "umap",dims = c(1, 2),pt.size=2)# + NoLegend()

###merge Mohammed et al with E5.0
mammal.anchors <- FindIntegrationAnchors(object.list = list(WT_samples,mouse_data_mohammed), dims = 1:20, anchor.features = 20000,k.filter = 50)
mammal.combined_all <- IntegrateData(anchorset = mammal.anchors, dims = 1:20, k.weight = 40)
DefaultAssay(mammal.combined_all) <- "integrated"
mammal.combined_all <- ScaleData(mammal.combined_all, verbose = FALSE)
mammal.combined_all <- RunPCA(mammal.combined_all, npcs = 20, verbose = FALSE)
mammal.combined_all <- RunUMAP(mammal.combined_all, reduction = "pca", dims = 1:20)
mammal.combined_all <- RunTSNE(mammal.combined_all, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined_all, reduction = "umap",dims = c(1, 2), split.by = 'species',pt.size=2) 

#perform unbiased clustering and label clusters by gene expression and ensure matches previosu annotations
mammal.combined_all_tissue <- FindNeighbors(mammal.combined_all_tissue, reduction = "umap", dims = 1:2)
mammal.combined_all_tissue <- FindClusters(mammal.combined_all_tissue, resolution = 0.2)
Idents(WT_samples)=Idents(mammal.combined_all_tissue)[1:length(Idents(WT_samples))]
clusters<-Idents(WT_samples)
clusters<-as.character(clusters)
clusters[clusters==1]<-rep("E5.0_VE",length(clusters[clusters==1]))
clusters[clusters==3]<-rep("E5.0_EPI",length(clusters[clusters==3]))
clusters[clusters==6]<-rep("E5.0_TE",length(clusters[clusters==6]))
Idents(WT_samples)<-clusters


#regenerate combined dataset with only WT samples
mammal.anchors <- FindIntegrationAnchors(object.list = list(WT_samples,mouse_data_mohammed), dims = 1:20, anchor.features = 4000,k.filter = 40)
mammal.combined_all <- IntegrateData(anchorset = mammal.anchors, dims = 1:20, k.weight = 40)
DefaultAssay(mammal.combined_all) <- "integrated"
mammal.combined_all <- ScaleData(mammal.combined_all,  verbose = FALSE)
mammal.combined_all <- ScaleData(mammal.combined_all, assay="RNA", features=rownames(mammal.combined_all),verbose = FALSE)
mammal.combined_all <- RunPCA(mammal.combined_all, npcs = 20, verbose = FALSE)
mammal.combined_all <- RunUMAP(mammal.combined_all, reduction = "pca", dims = 1:20)
mammal.combined_all <- RunTSNE(mammal.combined_all, reduction = "pca", dims = 1:20)

#set colours
cols <-  c(
  ## preimplantation
  "E4.5_Epiblast" = "#46E0F9",
  "E4.5_Primitive_endoderm" = "#FFC306",
  "E5.0_EPI" = "#0080F4",
  "E5.0_TE" = "#DF70F8",
  "E5.0_KO_TE" = "black",
  "E5.0_VE" = "#FF5A08",
  "E5.5_Epiblast" = "#003AF2",
  "E5.5_Visceral_endoderm" = "#AA0400",
  "E6.5_Epiblast" = "#4A23E8",
  "E6.5_ExE_ectoderm" = "#6BF400",
  "E6.5_Mesoderm" = "#00BF9F",
  "E6.5_Primitive_Streak" = "#00C12D",
  "E6.5_Visceral_endoderm" = "#5E0B0B"
  
)
DimPlot(mammal.combined_all, reduction = "umap",dims = c(1, 2), split.by = 'species',pt.size=2,col=cols)
DimPlot(mammal.combined_all, reduction = "umap",dims = c(1, 2),pt.size=2,cols=cols)
ggsave("Wildtype_all_relabbeled5_pca.pdf",width =16, height = 12)

###Identify developmentally relevant genes
invivo_markers<-FindAllMarkers(mammal.combined_all,assay="integrated")
celltypes<-unique(invivo_markers$cluster)
heatmap_invivo_genes<-c()
for (i in 1:length(celltypes)){
  temp_cell<-celltypes[i]
  temp_data<-invivo_markers[invivo_markers$cluster==temp_cell,]
  temp_data<-temp_data[temp_data$avg_log2FC >0.2,]
  heatmap_invivo_genes<-c(heatmap_invivo_genes,temp_data$gene[1:10])
  
}
levels(mammal.combined_all)<-levels(mammal.combined_all)[c(5,3,7,12,4,1,5,8,10,11,2,9)]
new_order<-levels(mammal.combined_all)[c(11,1,3,4,6,2,7,8,9,10,5,12)]
mammal.combined_all@active.ident <- factor(x = mammal.combined_all@active.ident, levels = new_order)

DoHeatmap(object = mammal.combined_all,features=heatmap_invivo_genes,group.by='ident') +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("unbiasedgenes_heatmap.pdf",width =16, height = 12)

#generate gene violin plots
dev_cols<-c("#E0E836","#FF9100","#F8766D","#FF0000","#46E0F9","#0080F4","#003AF2","#4A23E8","#00BF9F","#00C12D","#DF70F8","#6BF400")
specific_genes<-c("Pou5f1","Sox2","Otx2","Gata6","Sox7","Cer1","Hand1","Tfap2c",'Fgfr2')
for (i in 1:length(specific_genes)){
  VlnPlot(object=mammal.combined_all,assay="RNA",features=specific_genes[i],cols =dev_cols)+
    geom_boxplot()
  ggsave(paste("C:/",specific_genes[i],".pdf",sep=""),width =10, height = 7)
}


#generate just VE  violin plots
ve_mammal<-subset(x=mammal.combined_all,idents=c("E4.5_Primitive_endoderm","E5.0_VE","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm"))

DimPlot(ve_mammal, reduction = "umap",dims = c(1, 2),pt.size=3, cols=c("#E0E836","#FF9100","#F8766D","#FF0000")) 
ggsave("VE_umap.pdf",width =10, height = 8)
specific_genes<-c("Pou5f1","Cdh2","Sox17")
for (i in 1:length(specific_genes)){
  VlnPlot(object=ve_mammal,assay="RNA",features=specific_genes[i],cols =dev_cols,,pt.size=1.5)
  ggsave(paste("C:/",specific_genes[i],"_justVE.pdf",sep=""),width =10, height = 7)
}


# Pairwise differential gene expresion analysis 
WT_VE_4.5PEmarkers=FindMarkers(mammal.combined_all,"E5.0_WT_VE","E4.5_Primitive_endoderm")
WT_VE_5.5VEmarkers=FindMarkers(mammal.combined_all,"E5.0_WT_VE","E5.5_Visceral_endoderm")
WT_EPI_4.5EPImarkers=FindMarkers(mammal.combined_all,"E5.0_EPI","E4.5_Epiblast")
WT_EPI_5.5EPImarkers=FindMarkers(mammal.combined_all,"E5.0_EPI","E5.5_Epiblast")
WT_TE_5.5EPImarkers=FindMarkers(mammal.combined_all,"E5.0_TE","E6.5_ExE_ectoderm")

# Generate violin plots
WT_VE_4.5PEmarkers_wilcox=FindMarkers(mammal.combined_all,"E5.0_VE","E5.5_Visceral_endoderm",logfc.threshold = 0)
WT_VE_4.5PEmarkers_wilcoxmut <- mutate(WT_VE_4.5PEmarkers_wilcox, sig=ifelse(WT_VE_4.5PEmarkers_wilcox$avg_log2FC>0.3 | WT_VE_4.5PEmarkers_wilcox$avg_log2FC<(-0.3)  , yes = "Sig", no = "Nonsig")) 
WT_VE_4.5PEmarkers_wilcoxmut <- cbind(gene=rownames(WT_VE_4.5PEmarkers_wilcoxmut ), WT_VE_4.5PEmarkers_wilcoxmut ) 
WT_VE_4.5PEmarkers_wilcox_filter<-WT_VE_4.5PEmarkers_wilcoxmut
WT_VE_4.5PEmarkers_wilcox_filter[-log2(WT_VE_4.5PEmarkers_wilcox_filter$p_val)>30,2]<-0.000000001
WT_VE_4.5PEmarkers_wilcox_filter[WT_VE_4.5PEmarkers_wilcox_filter$avg_log2FC>1.5,3]<-1.5

volc = ggplot(WT_VE_4.5PEmarkers_wilcox_filter,aes(avg_log2FC, -log2(p_val)))+
  geom_point(aes(col=sig))+        scale_color_manual(values=c("grey", "red")) +
  geom_vline(xintercept=c(-0.3, 0.3), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  xlim(c(-1.5,1.5))+ylim(c(0,30))+
  theme(panel.background = element_blank())
volc+geom_text_repel(data=WT_VE_4.5PEmarkers_wilcox_filter[WT_VE_4.5PEmarkers_wilcox_filter$sig=="Sig",], aes(label=gene))
ggsave("WT_VE_4.5PEmarkers.pdf",width =10, height = 10)


#incorporate KO samples
KO_samples_lintest<-KO_samples
KO_samples_lintest <- FindNeighbors(KO_samples_lintest, reduction = "umap", dims = 1:2)
KO_samples_lintest <- FindClusters(KO_samples_lintest, resolution = 0.2)

#unbiased clustering and labelling of clusters by embryo annotation and gene expression
Idents(KO_samples)=Idents(KO_samples_lintest)
clusters<-Idents(KO_samples)
clusters<-as.character(clusters)
clusters[clusters==0]<-rep("E5.0_KO_VE",length(clusters[clusters==0]))
clusters[clusters==2]<-rep("E5.0_KO_EPI",length(clusters[clusters==2]))
clusters[clusters==1]<-rep("E5.0_KO_EPI",length(clusters[clusters==1]))
clusters[clusters==3]<-rep("E5.0_KO_TE",length(clusters[clusters==3]))
Idents(KO_samples)<-clusters
invitro_combined_IDENTS<- c(Idents(WT_samples),Idents(KO_samples))

#Merge with KO samples and remerge datasets using CCA
mammal.anchors <- FindIntegrationAnchors(object.list = list(mammal.combined_invitro,mouse_data_mohammed), dims = 1:20, anchor.features = 20000,k.filter = 50)
mammal.combined_all <- IntegrateData(anchorset = mammal.anchors, dims = 1:20, k.weight = 40)
DefaultAssay(mammal.combined_all) <- "integrated"
mammal.combined_all <- ScaleData(mammal.combined_all, verbose = FALSE)
mammal.combined_all <- RunPCA(mammal.combined_all, npcs = 20, verbose = FALSE)
mammal.combined_all <- RunUMAP(mammal.combined_all, reduction = "pca", dims = 1:20)
mammal.combined_all <- RunTSNE(mammal.combined_all, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined_all, reduction = "umap",dims = c(1, 2), split.by = 'species',pt.size=2) 


#reorder samples
lineage_order<-c("E4.5_Primitive_endoderm",
                 "E5.0_WT_VE",
                 "E5.0_KO_VE",
                 "E5.5_Visceral_endoderm",
                 "E6.5_Visceral_endoderm",
                 "E4.5_Epiblast",
                 "E5.0_EPI",
                 "E5.5_Epiblast",
                 "E6.5_Epiblast",
                 "E6.5_Mesoderm",
                 "E6.5_Primitive_Streak",
                 "E5.0_TE",
                 "E6.5_ExE_ectoderm")
new_order_ko<-levels(mammal.combined_all)[match(lineage_order, levels(mammal.combined_all))]
mammal.combined_all@active.ident <- factor(x = mammal.combined_all@active.ident, levels = new_order_ko)
dev_cols_KO<-c("#E0E836","#FF9100","#6BF400","#F8766D","#FF0000","#46E0F9","#0080F4","#003AF2","#4A23E8","#00BF9F","#00C12D","#DF70F8","#6BF400")
DimPlot(mammal.combined_all, reduction = "umap",dims = c(1, 2),cols=dev_cols_KO ,pt.size=2) 

#Pairwise differential gene expression analysis with KO embryos
KO_VE_5.0VEmarkers=FindMarkers(object=mammal.combined_all,ident.1="E5.0_KO_VE",ident.2="E5.0_WT_VE",logfc.threshold = 0)
KO_VE_4.5PEmarkers=FindMarkers(object=mammal.combined_all,"E5.0_KO_VE","E4.5_Primitive_endoderm",logfc.threshold = 0)
KO_VE_5.5VEmarkers=FindMarkers(object=mammal.combined_all,"E5.0_KO_VE","E5.5_Visceral_endoderm",logfc.threshold = 0)
write.csv(KO_VE_4.5PEmarkers,file = "KO_VE_4.5PEmarkers.csv",col.names=TRUE,row.names=TRUE)
write.csv(KO_VE_5.0VEmarkers,file = "KO_VE_5.0VEmarkers.csv",col.names=TRUE,row.names=TRUE)
write.csv(KO_VE_5.5VEmarkers,file = "KO_VE_5.5VEmarkers.csv",col.names=TRUE,row.names=TRUE)

#volcano plot of DE KO genes

KO_VE_4.5PEmarkers_wilcoxmut <- KO_VE_4.5PEmarkers
KO_VE_4.5PEmarkers_wilcoxmut$sig<- "NS"
KO_VE_4.5PEmarkers_wilcoxmut[KO_VE_4.5PEmarkers_wilcoxmut$avg_log2FC>0.3 & KO_VE_4.5PEmarkers_wilcoxmut$p_val<0.05,6]<- "UP"
KO_VE_4.5PEmarkers_wilcoxmut[KO_VE_4.5PEmarkers_wilcoxmut$avg_log2FC<(-0.3) & KO_VE_4.5PEmarkers_wilcoxmut$p_val<0.05,6]<- "down"
KO_VE_4.5PEmarkers_wilcoxmut <- cbind(gene=rownames(KO_VE_4.5PEmarkers_wilcoxmut ), KO_VE_4.5PEmarkers_wilcoxmut ) 
KO_VE_4.5PEmarkers_wilcoxmut[-log2(KO_VE_4.5PEmarkers_wilcoxmut$p_val)>15,2]<-0.00003051757

volc = ggplot(KO_VE_4.5PEmarkers_wilcoxmut,aes(avg_log2FC, -log2(p_val)))+
  geom_point(aes(col=sig))+        scale_color_manual(values=c("#E0E836","grey", "#6BF400")) +
  theme(panel.background = element_blank())
volc+geom_text_repel(data=KO_VE_4.5PEmarkers_wilcoxmut[KO_VE_4.5PEmarkers_wilcoxmut$sig=="UP" | KO_VE_4.5PEmarkers_wilcoxmut$sig=="down",], aes(label=gene))
ggsave("KO_VE_4.5PEmarkers_wilcoxmut.pdf",width =10, height = 10)


KO_VE_5.0PEmarkers_wilcoxmut <- KO_VE_5.0VEmarkers
KO_VE_5.0PEmarkers_wilcoxmut$sig<- "NS"
KO_VE_5.0PEmarkers_wilcoxmut[KO_VE_5.0PEmarkers_wilcoxmut$avg_log2FC>0.3 & KO_VE_5.0PEmarkers_wilcoxmut$p_val<0.05,6]<- "UP"
KO_VE_5.0PEmarkers_wilcoxmut[KO_VE_5.0PEmarkers_wilcoxmut$avg_log2FC<(-0.3) & KO_VE_5.0PEmarkers_wilcoxmut$p_val<0.05,6]<- "down"
KO_VE_5.0PEmarkers_wilcoxmut <- cbind(gene=rownames(KO_VE_5.0PEmarkers_wilcoxmut ), KO_VE_5.0PEmarkers_wilcoxmut ) 
volc = ggplot(KO_VE_5.0PEmarkers_wilcoxmut,aes(avg_log2FC, -log2(p_val)))+
  geom_point(aes(col=sig))+        scale_color_manual(values=c("#FF9100","grey", "#6BF400")) +
  theme(panel.background = element_blank())
volc+geom_text_repel(data=KO_VE_5.0PEmarkers_wilcoxmut[KO_VE_5.0PEmarkers_wilcoxmut$sig=="UP" | KO_VE_5.0PEmarkers_wilcoxmut$sig=="down",], aes(label=gene))+geom_text(data=KO_VE_5.0PEmarkers_wilcoxmut[KO_VE_5.0PEmarkers_wilcoxmut$gene=="Pou5f1",], aes(label=gene))
ggsave("KO_VE_5.0PEmarkers_wilcoxmut.pdf",width =10, height = 10)


KO_VE_5.5VEmarkers_wilcoxmut <- KO_VE_5.5VEmarkers
KO_VE_5.5VEmarkers_wilcoxmut$sig<- "NS"
KO_VE_5.5VEmarkers_wilcoxmut[KO_VE_5.5VEmarkers_wilcoxmut$avg_log2FC>0.3 & KO_VE_5.5VEmarkers_wilcoxmut$p_val<0.05,6]<- "UP"
KO_VE_5.5VEmarkers_wilcoxmut[KO_VE_5.5VEmarkers_wilcoxmut$avg_log2FC<(-0.3) & KO_VE_5.5VEmarkers_wilcoxmut$p_val<0.05,6]<- "down"
KO_VE_5.5VEmarkers_wilcoxmut <- cbind(gene=rownames(KO_VE_5.5VEmarkers_wilcoxmut ), KO_VE_5.5VEmarkers_wilcoxmut ) 
volc = ggplot(KO_VE_5.5VEmarkers_wilcoxmut,aes(avg_log2FC, -log2(p_val)))+
  geom_point(aes(col=sig))+        scale_color_manual(values=c("#F8766D","grey", "#6BF400")) +
  theme(panel.background = element_blank())
volc+geom_text_repel(data=KO_VE_5.5VEmarkers_wilcoxmut[KO_VE_5.5VEmarkers_wilcoxmut$sig=="UP" | KO_VE_5.5VEmarkers_wilcoxmut$sig=="down",], aes(label=gene))

ggsave("KO_VE_5.5VEmarkers_wilcoxmut.pdf",width =10, height = 10)


###generate violin plots 

specific_genes<-c('Cited1',"Gpc3","Apom","Plet1","Pou5f1","Sox17","Pdgfra","Lama1","Cer1","Otx2","Apoa1","Dab2",'Fabp3',"Glipr2","Ttr","Lamb1","Lamc1")
for (i in 1:length(specific_genes)){
  VlnPlot(object=ve_mammal_KO,assay="RNA",features=specific_genes[i],cols =dev_cols_KO,pt.size=1.5)
  ggsave(paste("C:/Users/Dylan/Documents/cambridge/forDylan/mouse_analysis/violinsKO/",specific_genes[i],"_justVE.pdf",sep=""),width =10, height = 7)
}


for (i in 1:length(specific_genes)){
  FeaturePlot(object=ve_mammal_KO,features=specific_genes[i],pt.size=4)
  ggsave(paste("C:/Users/Dylan/Documents/cambridge/forDylan/mouse_analysis/umapgenes_KO/",specific_genes[i],"_justVE.pdf",sep=""),width =10, height = 7)
}

##WT signalling

### Download KEGG datasets 
ve_mammal<-subset(x=mammal.combined_all,idents=c("E4.5_Primitive_endoderm","E5.0_WT_VE","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm"))
ve_mammal<-ScaleData(ve_mammal, assay = 'RNA', verbose = FALSE)
ensembl = useDataset("mmusculus_gene_ensembl", mart = ensembl)

mmu_list <- getGenesets(org = "mmu", db = "kegg", cache = TRUE, return.type="list")
pathway_names<-c("WNT","TGFB","MAPK","PI3K-AKT","JAK-STAT","Pluripotency","HIPPO","MTOR","NOTCH","Apoptosis","Oxphos","Glycolosis")
pathway_list<-list(mmu_list$mmu04310,
                   mmu_list$mmu04350,
                   mmu_list$mmu04010,
                   mmu_list$mmu04151,
                   mmu_list$mmu04630,
                   mmu_list$mmu04550,
                   mmu_list$mmu04392,
                   mmu_list$mmu04150,
                   mmu_list$mmu04330,
                   mmu_list$mmu04210,
                   mmu_list$mmu00190,
                   mmu_list$mmu00010)
pathway_list_genesymbol<-pathway_list

for (i in 1:length(pathway_list)){
  pathway_list_genesymbol[i]<-getBM(attributes = c("external_gene_name"), 
                                    filters = "entrezgene_id", 
                                    values = pathway_list[i], 
                                    mart = ensembl)
}

#Calculate module scores for KEGG pathways

ve_mammal<-AddModuleScore(ve_mammal,features=pathway_list_genesymbol,assay="integrated")

#Generate heatmaps for signalling pathways

#BMP signalling
DoHeatmap(object = ve_mammal,assay="RNA",features=pathway_list_genesymbol[3],group.by='ident') +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("bmpnodal_wt_heatmap.pdf",width =16, height = 12)

#heatmap just VE lineages
specific_genes<-c('Smad1',"Smad2","Smad3","Smad4","Smad5","Dkk1","Bmp2","BMP4","Nodal","Pitx2")
DoHeatmap(object = ve_mammal,assay="RNA",features=specific_genes,group.by='ident') +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("bmpnodal_wt_heatmap.pdf",width =16, height = 12)

#heatmap all lineages
mammal.combined_all<-ScaleData(mammal.combined_all, assay = 'RNA', verbose = FALSE)
DoHeatmap(object = mammal.combined_all,assay="RNA",features=specific_genes,group.by='ident') +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("bmpnodal_all_heatmap.pdf",width =16, height = 12)

ve_mammal_KO<-subset(x=mammal.combined_all,idents=c("E4.5_Primitive_endoderm","E5.0_WT_VE","E5.0_KO_VE","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm"))
ve_mammal_KO<-ScaleData(ve_mammal_KO, assay = 'RNA', verbose = FALSE)

#SLC genes obtained from gene ontology
slcs<-read.csv('slcs.csv',stringsAsFactors = F)
specific_genes<-c("Cubn","Slc16a1","Apoa1","Apob","Soat2",'Slc16a1',"Slc16a10","Pla2g12b","Fasn","Soat2","Apoa1","Apoa4","Apob")
for (i in 1:length(specific_genes)){
  VlnPlot(object=ve_mammal_KO,assay="RNA",features=specific_genes[i],cols =dev_cols_KO,pt.size=1.5)
  ggsave(paste("C:/",specific_genes[i],"_justVE.pdf",sep=""),width =10, height = 7)
}

#KO signalling
mammal.combined_all<-FindVariableFeatures(mammal.combined_all)
expressed_genes<-VariableFeatures(mammal.combined_all)
EG<-apply(mammal.combined_all@assays$RNA@data[,Idents(mammal.combined_all)=="E5.0_WT_VE"],1,mean)
SLCEG<-EG[grep("Slc",names(EG))]
expressedSLC<-names(SLCEG[SLCEG>0.5])

DoHeatmap(object = ve_mammal_KO,features=expressedSLC,group.by='ident',disp.min=-1,disp.max=1) +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))

# Signalling module scores for WT and KO samples
mammal.combined_all<-AddModuleScore(mammal.combined_all,features=pathway_list_genesymbol,assay="RNA")
sig_heatmap<-FetchData(object = mammal.combined_all, vars = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6", "Cluster7", "Cluster8", "Cluster9", "Cluster10", "Cluster11", "Cluster12", "Cluster13", "Cluster14", "Cluster15"))
sig_heatmap_av<-matrix(,nrow=length(unique(Idents(mammal.combined_all))),ncol=length(pathway_list_genesymbol))
tissues<-lineage_order
for (j in 1:length(tissues)){
  temp_data<-sig_heatmap[Idents(mammal.combined_all)==tissues[j],]
  temp_data_av<-colMeans(temp_data)
  sig_heatmap_av[j,]<-temp_data_av
  
}
rownames(sig_heatmap_av)<-tissues
colnames(sig_heatmap_av)<-pathway_names
sig_heatmap_av<-t(sig_heatmap_av)
sig_heatmap_av_scale = t(scale(t(sig_heatmap_av)))
col_fun = colorRamp2(c(-2, 0, 2), c("#00B0F0", "white", "#FF0000"))

pdf("Signalling_alllin2_heatmap.pdf", height=10,width=12, useDingbats = FALSE)
Heatmap(sig_heatmap_av_scale,col_fun,cluster_columns = FALSE)
dev.off()

#EMT and basement membrane violin plots

specific_genes<-EMT_list$basement
for (i in 1:length(specific_genes)){
  VlnPlot(object=ve_mammal_KO,assay="RNA",features=specific_genes[i],cols =dev_cols_KO,pt.size=1.5)
  ggsave(paste("C:/",specific_genes[i],"_justVE.pdf",sep=""),width =10, height = 7)
}

specific_genes<-c("Lama1","Lamb1","Lamc1","Timp1","Timp3","Col4a1","Col4a2","Hspg2","Loxl2","Itgb1","Lad1","P3h1","Anxa2","Hspg2","Fn1","Nid2","Dag1")
DoHeatmap(object = ve_mammal_KO,features=specific_genes,group.by='ident') +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("basementKO_heatmap.pdf",width =10, height = 7)

#Junctions and ECM

EMT_list<-read.csv('EMT_GO.csv',stringsAsFactors = F)
pathway_names<-c("Focal adhesion","Adherens junctions","Tight junctions","Gap junctions","ECM RECEPTOR","EMT","Basement membrane")
pathway_list<-list(mmu_list$mmu04510,
                   mmu_list$mmu04520,
                   mmu_list$mmu04530,
                   mmu_list$mmu04540,
                   mmu_list$mmu04512,
                   EMT_list[,1],
                   EMT_list[lapply(EMT_list[,2],nchar)>0,2])
pathway_list_genesymbol<-pathway_list


for (i in 1:(length(pathway_list)-2)){
  pathway_list_genesymbol[i]<-getBM(attributes = c("external_gene_name"), 
                                    filters = "entrezgene_id", 
                                    values = pathway_list[i], 
                                    mart = ensembl)
}

#Just VE lineages

ve_mammal<-subset(x=mammal.combined_all,idents=c("E4.5_Primitive_endoderm","E5.0_WT_VE","E5.0_KO_VE","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm"))
ve_mammal<-AddModuleScore(ve_mammal,features=pathway_list_genesymbol,assay="integrated")
sig_heatmap<-FetchData(object = ve_mammal, vars = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6", "Cluster7"))
sig_heatmap_av<-matrix(,nrow=length(unique(Idents(ve_mammal))),ncol=length(pathway_list_genesymbol))
tissues<-unique(Idents(ve_mammal))
tissues<-c("E4.5_Primitive_endoderm","E5.0_WT_VE","E5.0_KO_VE","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm")
for (j in 1:length(tissues)){
  temp_data<-sig_heatmap[Idents(ve_mammal)==tissues[j],]
  temp_data_av<-colMeans(temp_data)
  sig_heatmap_av[j,]<-temp_data_av
  
}
rownames(sig_heatmap_av)<-tissues
colnames(sig_heatmap_av)<-pathway_names
sig_heatmap_av<-t(sig_heatmap_av)
sig_heatmap_av_scale = t(scale(t(sig_heatmap_av)))

col_fun = colorRamp2(c(-1, 0, 1), c("#00B0F0", "white", "#FF0000"))
pdf("junctions_basment_ve__heatmap.pdf", height=10,width=12, useDingbats = FALSE)

Heatmap(sig_heatmap_av_scale,col_fun,cluster_columns = FALSE)
dev.off()

#FGF signalling
specific_genes<-c('Fgfr1',"Fgfr2","Pdgfra","Fgf2","Fgf3","Fgf4","Fgf5","Fgf6","Fgf8","Fgf10","Mapk1")
for (i in 1:length(specific_genes)){
  VlnPlot(object=ve_mammal_KO,assay="RNA",features=specific_genes[i],cols =dev_cols_KO,pt.size=1.5)
  ggsave(paste("C:/Users/Dylan/Documents/cambridge/forDylan/mouse_analysis/violinsKO_fgf/",specific_genes[i],"_justVE.pdf",sep=""),width =10, height = 7)
}

DoHeatmap(object = mammal.combined_all,features=specific_genes,group.by='ident') +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("FGF_all_heatmap.pdf",width =16, height = 12)

## EMT signalling
EMT_list<-read.csv('EMT_GO.csv',stringsAsFactors = F)
pathway_names<-c("Focal adhesion","Adherens junctions","Tight junctions","Gap junctions","ECM RECEPTOR","EMT","Basement membrane")
pathway_list<-list(mmu_list$mmu04510,
                   mmu_list$mmu04520,
                   mmu_list$mmu04530,
                   mmu_list$mmu04540,
                   mmu_list$mmu04512,
                   EMT_list[,1],
                   EMT_list[lapply(EMT_list[,2],nchar)>0,2])
pathway_list_genesymbol<-pathway_list



for (i in 1:(length(pathway_list)-2)){
  pathway_list_genesymbol[i]<-getBM(attributes = c("external_gene_name"), 
                                    filters = "entrezgene_id", 
                                    values = pathway_list[i], 
                                    mart = ensembl)
}

ve_mammal<-subset(x=mammal.combined_all,idents=c("E4.5_Primitive_endoderm","E5.0_VE","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm"))
ve_mammal<-AddModuleScore(ve_mammal,features=pathway_list_genesymbol,assay="integrated")
sig_heatmap<-FetchData(object = ve_mammal, vars = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6", "Cluster7"))
sig_heatmap_av<-matrix(,nrow=length(unique(Idents(ve_mammal))),ncol=length(pathway_list_genesymbol))
tissues<-unique(Idents(ve_mammal))
tissues<-c("E4.5_Primitive_endoderm","E5.0_VE","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm")
for (j in 1:length(tissues)){
  temp_data<-sig_heatmap[Idents(ve_mammal)==tissues[j],]
  temp_data_av<-colMeans(temp_data)
  sig_heatmap_av[j,]<-temp_data_av
  
}

rownames(sig_heatmap_av)<-tissues
colnames(sig_heatmap_av)<-pathway_names
sig_heatmap_av<-t(sig_heatmap_av)
sig_heatmap_av_scale = t(scale(t(sig_heatmap_av)))
col_fun = colorRamp2(c(-1, 0, 1), c("#00B0F0", "white", "#FF0000"))
pdf("Signalling_ve__heatmap.pdf", height=10,width=12, useDingbats = FALSE)

#TGFB
specific_genes<-c('Smad1',"Smad2","Smad3","Smad4","Smad5","Dkk1","Bmp2","Nodal","Pitx2")
for (i in 1:length(specific_genes)){
  VlnPlot(object=ve_mammal,assay="RNA",features=specific_genes[i],cols =dev_cols,pt.size=1.5)
  ggsave(paste("C:/Users/Dylan/Documents/cambridge/forDylan/mouse_analysis/violinsWT_BMPNODAL/",specific_genes[i],"_justVE.pdf",sep=""),width =10, height = 7)
}

#Transporters
transporter_list<-read.csv('transporter_terms.csv',stringsAsFactors = F)
pathway_names<-c("protein_transporter_activity","lipid_transport_activity",
                 "Active transmembrane transporter activity","Ion transport")
pathway_list<-list(transporter_list[,1],
                   transporter_list[,2],
                   transporter_list[,5],
                   transporter_list[,7])
pathway_list_genesymbol<-pathway_list
mammal.combined_all<-AddModuleScore(mammal.combined_all,features=pathway_list_genesymbol,assay="RNA")
sig_heatmap<-FetchData(object = mammal.combined_all, vars = c("Cluster1", "Cluster2", "Cluster3", "Cluster4"))
sig_heatmap_av<-matrix(,nrow=length(unique(Idents(mammal.combined_all))),ncol=length(pathway_list_genesymbol))
tissues<-lineage_order

for (j in 1:length(tissues)){
  temp_data<-sig_heatmap[Idents(mammal.combined_all)==tissues[j],]
  temp_data_av<-colMeans(temp_data)
  sig_heatmap_av[j,]<-temp_data_av
  
}


#Correlation analysis

marmoset.combined.av<-matrix(,nrow=nrow(mammal.combined_all[["integrated"]]@scale.data),ncol=length(unique(Idents(mammal.combined_all))))
lineage_order<-c("E4.5_Primitive_endoderm",
                 "E5.0_WT_VE",
                 "E5.0_KO_VE",
                 "E5.5_Visceral_endoderm",
                 "E6.5_Visceral_endoderm",
                 "E4.5_Epiblast",
                 "E5.0_EPI",
                 "E5.5_Epiblast",
                 "E6.5_Epiblast",
                 "E6.5_Mesoderm",
                 "E6.5_Primitive_Streak",
                 "E5.0_TE",
                 "E6.5_ExE_ectoderm")

for (i in 1:length(lineage_order)){
  tis_temp<-lineage_order[i]
  marmoset.combined.av[,i]<-rowMeans(mammal.combined_all[["integrated"]]@scale.data[,which(Idents(mammal.combined_all)==tis_temp)],na.rm=TRUE)
  
  
}
colnames(marmoset.combined.av)<-lineage_order
rownames(marmoset.combined.av)<-rownames(mammal.combined_all[["integrated"]]@scale.data)
M<-cor(as.matrix(marmoset.combined.av))
M_2<-(M+1)/2
pdf("Correlation_all__heatmap.pdf", height=10,width=12, useDingbats = FALSE)

ve_mammal_KO<-subset(x=mammal.combined_all,idents=c("E4.5_Primitive_endoderm","E5.0_WT_VE","E5.0_KO_VE","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm"))

sig_heatmap<-FetchData(object = ve_mammal_KO, vars = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6", "Cluster7",
                                                       "Cluster8","Cluster9","Cluster10","Cluster11","Cluster12"))
sig_heatmap_av<-matrix(,nrow=length(unique(Idents(ve_mammal_KO))),ncol=length(pathway_list_genesymbol))
tissues<-c("E4.5_Primitive_endoderm",
           "E5.0_WT_VE",
           "E5.0_KO_VE",
           "E5.5_Visceral_endoderm",
           "E6.5_Visceral_endoderm")
for (j in 1:length(tissues)){
  temp_data<-sig_heatmap[Idents(ve_mammal_KO)==tissues[j],]
  temp_data_av<-colMeans(temp_data)
  sig_heatmap_av[j,]<-temp_data_av
  
}

rownames(sig_heatmap_av)<-tissues
colnames(sig_heatmap_av)<-pathway_names
sig_heatmap_av<-t(sig_heatmap_av)

sig_heatmap_av_scale = t(scale(t(sig_heatmap_av)))

col_fun = colorRamp2(c(-2, 0, 2), c("#00B0F0", "white", "#FF0000"))

pdf("transport_heatmap.pdf", height=10,width=12, useDingbats = FALSE)
Heatmap(sig_heatmap_av_scale,col_fun,cluster_columns = FALSE)
dev.off()

#Transciption factor analysis

tf_list <-read.csv('Mus_musculus_TF.csv',stringsAsFactors = F)
ko_counts<-ve_mammal_KO[["integrated"]]@scale.data
tf_KO<-ko_counts[which(rownames(ko_counts) %in% tf_list$Symbol),]
tissues_ve<-unique(Idents(ve_mammal_KO))
ko_avTF<-matrix(,nrow=nrow(tf_KO),ncol=length(tissues_ve))
for (g in 1:length(tissues_ve)){
  ko_avTF[,g]<-rowMeans(tf_KO[,Idents(ve_mammal_KO)== tissues_ve[g]])
}
rownames(ko_avTF)<-rownames(tf_KO)
colnames(ko_avTF)<-tissues_ve  
write.csv(ko_avTF,file = "TF_VEtissueaverages.csv",col.names=TRUE,row.names=TRUE)