library(readr)
library(readxl)
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(sctransform)
library(harmony)
library(ComplexHeatmap)
library(clustree)

########## Initialize Functions ##########
analyze_seurat <- function(seurat_object){
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object, approx = F)
  seurat_object <- RunHarmony(seurat_object, group.by.vars="PatientID", max.iter.cluster=500)
  seurat_object <- FindNeighbors(seurat_object, reduction = "harmony")
  seurat_object <- FindClusters(seurat_object, resolution = seq(0,1,0.1), random.seed=123)
  seurat_object <- FindClusters(seurat_object, resolution=0.5, random.seed=123, graph.name = 'PSI_snn')
  seurat_object <- RunUMAP(seurat_object, reduction = "harmony", seed.use = 123, dims=1:10, verbose=FALSE, a=1, b=0.75)
  seurat_object <- RunTSNE(seurat_object, reduction = "harmony")
}

groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}


cluster_membership <- function(x,y){
  a <- table(x,y)
  b <- t(a/rowSums(a))
  p <- pheatmap::pheatmap(mat = b, color = colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100), 
                          cluster_cols = FALSE, cellheight=15, cellwidth = 15)
  return(p)
}


cluster_membership_row <- function(x,y){
  a <- t(table(x,y))
  b <- a/rowSums(a)
  p <- pheatmap::pheatmap(mat = b, color = colorRampPalette(brewer.pal(n = 7, name = "Greens"))(100), 
                          cluster_cols = FALSE, cellheight=15, cellwidth = 15)
  return(p)
}

ColAssign <- function(Var,palettes="Classic 20"){
  require(ggthemes);require(RColorBrewer)
  pal <- tableau_color_pal(palette = palettes,direction = 1,type="regular")
  if (length(Var) > 20) {
    palOut <- colorRampPalette(pal(20))(length(Var))
    names(palOut) <- Var
  } else if (length(Var) == 20) {
    palOut <- pal(20)
    names(palOut) <- Var
  } else if (length(Var) < 20) {
    palOut <- pal(20)
    palOut <- setdiff(palOut,c("#7f7f7f","#c7c7c7"))# remove grey colors
    #palOut <- sample(palOut)
    palOut <- c(palOut,c("#7f7f7f","#c7c7c7"))
    palOut <- palOut[1:length(Var)]
    names(palOut) <- Var
  }
  return(palOut)
}

########## Load & Process Data ##########
sc_data <- read.csv("/MDA Blood + BM_CD4+, CD8+_Raw Data_Table.csv")

clinical_data <- read_excel("isoplexis_metadata_noID.xlsx", sheet = 1)


names(sc_data) <- gsub("\\.","",names(sc_data))            # remove the '.' from column names
colnames(sc_data)[colnames(sc_data)=='X'] <- 'Cell_ID'     # rename X as Cell_ID
sc_data = subset(sc_data, select = -c(Stimulation))        # remove stimulation column - not necessary
colnames(sc_data)[colnames(sc_data) == "InternalDonor"] <- "PatientID"  # change to PatientID for primary key

PatientID <- c(sc_data$PatientID)
sc_data$CellID <- paste0("cell",1:111420)

# Create df for metadata
metadata <- sc_data[,c("CellID", "CellSubset")]

# Remove variables from count dataframe
drop <- c("CellSubset")
sc_data = sc_data[,!(names(sc_data)%in% drop)]

# move cellid to front and create as row names
sc_data <- sc_data %>% 
  select(CellID, everything())
sc_data_new <- sc_data[,-1]
rownames(sc_data_new) <- sc_data[,1]

# Create matrix for PSI counts
sc_data_new <- sc_data_new[,-1]
sc_matrix <- t(sc_data_new)

# Create seurat object
object <- CreateSeuratObject(counts = sc_matrix, assay = "PSI")


# Create metadata dataframe and add to Seurat object
# Use the provided key to get sample and time information
clinical_data <- clinical_data %>%
  mutate(Sample = recode(Sample, Blood = 'PB', BM = 'BM'))
#subset clinical data that is desired for the object
clinical_sub <- clinical_data[,c(1:6)] #add any other data that is desired
clinical_sub$Time[clinical_sub$Time == "2 Month"] <- "Two Month"

metadata_clin <- left_join(metadata, clinical_sub, by=c("Patient"="PatientID", "Time", "Sample", "CellSubset"))

metadata_clin <- metadata_clin %>% select(CellID, everything())
metadata_new <- metadata_clin[,-1]
rownames(metadata_new) <- metadata_clin[,1]
object <- AddMetaData(object, metadata_new)
DefaultAssay(object) <- "PSI"

########## Seurat Pipeline No Harmony ##########
# Run standard Seurat pipeline
# This set of code is without harmony to see what the object looks like w/o batch correction
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, approx = F)
object <- FindNeighbors(object)
object <- FindClusters(object, resolution = seq(0,1,0.1), random.seed=123)
object <- FindClusters(object, resolution=0.5, random.seed=123, graph.name = 'PSI_snn')
object <- RunUMAP(object, seed.use = 123, dims=1:10, verbose=FALSE)
object <- RunTSNE(object)

pdf("Documents/Isoplexis/QC/all_samples_no_harmony_umap.pdf", width=10, height=5)
DimPlot(object, group.by = "Patient") #split is patient 1-10, and Patient 11-21
dev.off()

########## Seurat Pipeline Harmony ##########
# run harmony & normalize
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, approx = F)
object <- RunHarmony(object, group.by.vars="PatientID", max.iter.cluster=500)
object <- FindNeighbors(object, reduction = "harmony")
object <- FindClusters(object, resolution = seq(0,1,0.1), random.seed=123)
object <- FindClusters(object, resolution=0.5, random.seed=123, graph.name = 'PSI_snn')
object <- RunUMAP(object, reduction = "harmony", seed.use = 123, dims=1:10, verbose=FALSE)
object <- RunTSNE(object, reduction = "harmony")

pdf("Documents/Isoplexis/QC/all_samples_harmony_umap.pdf", width=10, height=5)
DimPlot(object, group.by = "Patient") 
dev.off()

pdf("Documents/Isoplexis/QC/all_samples_cellsubset.pdf", width=6, height=5)
DimPlot(object, group.by = "CellSubset") 
dev.off()

pdf("Documents/Isoplexis/QC/all_samples_sample.pdf", width=6, height=5)
DimPlot(object, group.by = "Sample") 
dev.off()

pdf("Documents/Isoplexis/QC/all_samples_split_sample.pdf", width=10, height=5)
DimPlot(object, split.by = "Sample") 
dev.off()

pdf("Documents/Isoplexis/QC/all_samples_cellsubset_sample.pdf", width=10, height=5)
DimPlot(object, group.by = "CellSubset", split.by = "Sample") 
dev.off()

########## Remove Outliers ##########
# 935 cells w/ Count > 5000, 2805 cells w/ Count > 4000, 10286 cells w/ Count > 3000
# 5262 cells / Count > 3500
all_cells <- subset(object, subset=nCount_PSI < 3000)
all_cells <- NormalizeData(all_cells)
all_cells <- FindVariableFeatures(all_cells)
all_cells <- ScaleData(all_cells)
all_cells <- RunPCA(all_cells, approx = F)
all_cells <- RunHarmony(all_cells, group.by.vars="PatientID", max.iter.cluster=500)
all_cells <- FindNeighbors(all_cells, reduction = "harmony")
all_cells <- FindClusters(all_cells, resolution = seq(0,1,0.1), random.seed=123)
all_cells <- FindClusters(all_cells, resolution=0.5, random.seed=123, graph.name = 'PSI_snn')
all_cells <- RunUMAP(all_cells, reduction = "harmony", seed.use = 123, dims=1:10, verbose=FALSE, a=1,b=0.75)
all_cells <- RunTSNE(all_cells, reduction = "harmony")


rm(object)
########## Fix Data/Add Scores ##########


effector <- c("GranzymeB", "IFNg", "MIP1a", "Perforin", "TNFa", "TNFb")
stimulatory <- c("GMCSF", "IL2", "IL5", "IL7", "IL8", "IL9", "IL12", "IL15", "IL21")
chemoattractive <- c("CCL11", "IP10", "MIP1b", "RANTES")
regulatory <- c("IL4", "IL10", "IL13", "IL22", "TGFb1", "sCD137", "sCD40L")
inflammatory <- c("IL1b", "IL6", "IL17A", "IL17F", "MCP1", "MCP4")

cyto_groups <- list("effector"=effector, 
                    "stimulatory"=stimulatory, 
                    "chemoattractive"=chemoattractive, 
                    "regulatory"=regulatory, 
                    "inflammatory"=inflammatory)

all_cells$effector_score <- rowSums(FetchData(all_cells, effector), na.rm = TRUE)
all_cells$stimulatory_score <- rowSums(FetchData(all_cells, stimulatory), na.rm = TRUE)
all_cells$chemoattractive_score <- rowSums(FetchData(all_cells, chemoattractive), na.rm = TRUE)
all_cells$regulatory_score <- rowSums(FetchData(all_cells, regulatory), na.rm = TRUE)
all_cells$inflammatory_score <- rowSums(FetchData(all_cells, inflammatory), na.rm = TRUE)

all_cells$monocytic <- "No"
all_cells$monocytic[all_cells$Patient %in% c("MDADD0007", "MDADD0008")] <- "Mono"


saveRDS(all_cells, "Documents/Isoplexis/Objects/scIso_filtered.rds")

########## Find Appropriate Resolution for All Cells ##########
clustree <- clustree(all_cells) + labs(title="Clustree")
plot_umap <- function(res){
  res_name = paste0("PSI_snn_res.", res)
  umap <- DimPlot(all_cells, group.by = res_name,label=T, raster=F) + 
    labs(title= paste0("UMAP for ", res_name))
  return(umap)
}
umap_plots <- map(seq(0,1,0.1), plot_umap)
pdf("Documents/Isoplexis/UMAPs/all_cells.pdf",width = 8,height=7)
print(ggarrange(plotlist=umap_plots, ncol=1,nrow=1))
print(ggarrange(clustree))
dev.off()

DimPlot(all_cells, group.by="CellSubset",cols=c("pink","seagreen"),raster=F, shuffle=T)+
  NoAxes()+ggtitle("All Cells by Cell Subset")+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 8),
        legend.key.height = unit(0.25,"cm"),
        legend.key.width = unit(0.25,"cm"))+NoAxes()
ggsave("Documents/Isoplexis/UMAPs/all_cells_subset.pdf",height=5,width=5,unit="in")

DimPlot(all_cells, group.by="PSI_snn_res.0.3",cols=ColAssign(unique(all_cells$PSI_snn_res.0.3)),raster=F, shuffle=T, label=T)+
  NoAxes()+ggtitle("All Cells at Resolution 0.3")+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 8),
        legend.key.height = unit(0.25,"cm"),
        legend.key.width = unit(0.25,"cm"))+NoAxes()
ggsave("Documents/Isoplexis/UMAPs/all_cells_res0.3.pdf",height=5,width=5,unit="in")

pdf("Documents/Isoplexis/ClusterMembership/cm_all_cells_subset.pdf", width=3, height=1.5)
print(cluster_membership(all_cells$PSI_snn_res.0.3, all_cells$CellSubset))
dev.off()

########## Heatmap ##########
cytokines <- unique(rownames(iso_cd4))
ckgroups <- c("Chemoattractive","Stimulatory","Effector","Effector","Regulatory","Stimulatory","Regulatory","Stimulatory",
              "Inflammatory","Inflammatory","Inflammatory","Stimulatory","Stimulatory","Regulatory","Regulatory","Stimulatory",
              "Inflammatory","Stimulatory","Stimulatory","Stimulatory","Chemoattractive","Inflammatory","Inflammatory","Effector",
              "Chemoattractive","Effector","Chemoattractive","Regulatory","Regulatory","Regulatory","Effector","Effector")
names(ckgroups) <- cytokines

rightannodf = data.frame(Groups = ckgroups)
rightanno = rowAnnotation(df = rightannodf, simple_anno_size = unit(2, "mm"),
  col = list(Groups = c("Chemoattractive" = "purple", "Stimulatory" = "deepskyblue", "Effector" = "darkseagreen3", "Regulatory" = "gold2", 
                        "Inflammatory" = "indianred2")))


all_cells_counts <- GetAssayData(object = all_cells, slot = "counts")
counts_avg <- groupMeans(all_cells_counts, groups = all_cells$PSI_snn_res.0.3)
counts_avg <- t(scale(t(counts_avg)))
celltype_cluster <- FetchData(all_cells, vars=c('PSI_snn_res.0.3', 'CellSubset'))
celltype_cluster <- celltype_cluster %>% group_by(PSI_snn_res.0.3, CellSubset) %>% dplyr::summarize(n=dplyr::n()) %>% 
  pivot_wider(names_from = CellSubset, values_from = n) %>% data.frame()
rownames(celltype_cluster) <- celltype_cluster$PSI_snn_res.0.3
celltype_cluster <- celltype_cluster[,-1]
celltype_cluster <- celltype_cluster/rowSums(celltype_cluster)
colnames(celltype_cluster) <- c("CD4", "CD8")
celltype_cluster <- celltype_cluster[colnames(counts_avg),]
ha <- HeatmapAnnotation(CellSubset = anno_barplot(celltype_cluster, gp = gpar(fill=c("pink","seagreen"))))
lgd_anno <- Legend(labels = c("CD4", "CD8"), legend_gp = gpar(fill = c("pink", "seagreen")), title = "Cell Subset")

pdf("Documents/Isoplexis/Heatmaps/avg_exp_heatmap_res0.3.pdf", height=8, width=5)
ht <- Heatmap(counts_avg,  width = ncol(counts_avg)*unit(5, "mm"), height = nrow(counts_avg)*unit(5, "mm"),
        heatmap_legend_param = list(title="Z-score"), top_annotation = ha, cluster_columns = FALSE,
        right_annotation = rightanno)
draw(ht, annotation_legend_list = list(lgd_anno))
dev.off()


all_cells_counts <- GetAssayData(object = all_cells, slot = "counts")
counts_avg <- groupMeans(all_cells_counts, groups = all_cells$PSI_snn_res.0.5)
counts_avg <- t(scale(t(counts_avg)))
celltype_cluster <- FetchData(all_cells, vars=c('PSI_snn_res.0.5', 'CellSubset'))
celltype_cluster <- celltype_cluster %>% group_by(PSI_snn_res.0.5, CellSubset) %>% dplyr::summarize(n=dplyr::n()) %>% 
  pivot_wider(names_from = CellSubset, values_from = n) %>% data.frame()
rownames(celltype_cluster) <- celltype_cluster$PSI_snn_res.0.5
celltype_cluster <- celltype_cluster[,-1]
celltype_cluster <- celltype_cluster/rowSums(celltype_cluster)
colnames(celltype_cluster) <- c("CD4", "CD8")
celltype_cluster <- celltype_cluster[colnames(counts_avg),]
ha <- HeatmapAnnotation(CellSubset = anno_barplot(celltype_cluster,gp = gpar(fill=c("pink","seagreen"))))
lgd_anno <- Legend(labels = c("CD4", "CD8"), legend_gp = gpar(fill = c("pink", "seagreen")),title = "Cell Subset")

pdf("Documents/Isoplexis/Heatmaps/all_cells_heatmap_res0.5.pdf", height=8, width=6)
ht <- Heatmap(counts_avg,  width = ncol(counts_avg)*unit(5, "mm"), height = nrow(counts_avg)*unit(5, "mm"),
              heatmap_legend_param = list(title="Z-score"), top_annotation = ha, cluster_columns = FALSE,
              right_annotation = rightanno)
draw(ht, annotation_legend_list = list(lgd_anno))
dev.off()


########## Subsets ##########
#Split objects based on the sampling site and run Seurat pipeline with harmony
iso_bm <- subset(all_cells, subset=Sample %in% "BM")
iso_pb <- subset(all_cells, subset=Sample %in% "PB")

iso_bm <- FindVariableFeatures(iso_bm)
iso_bm <- analyze_seurat(iso_bm)

iso_pb <- FindVariableFeatures(iso_pb)
iso_pb <- analyze_seurat(iso_pb)


iso_cd8 <- subset(all_cells, subset = CellSubset %in% "CD8+")
iso_cd4 <- subset(all_cells, subset = CellSubset %in% "CD4+")

iso_cd8 <- FindVariableFeatures(iso_cd8)
iso_cd8 <- analyze_seurat(iso_cd8)

iso_cd4 <- FindVariableFeatures(iso_cd4)
iso_cd4 <- analyze_seurat(iso_cd4)

saveRDS(iso_cd4, "Documents/Isoplexis/Objects/scIso_cd4.rds")
saveRDS(iso_cd8, "Documents/Isoplexis/Objects/scIso_cd8.rds")

cd4_base <- subset(iso_cd4, subset = Time %in% c("Baseline"))
cd8_base <- subset(iso_cd8, subset = Time %in% c("Baseline"))

saveRDS(cd4_base, "Documents/Isoplexis/Objects/cd4_base.rds")
saveRDS(cd8_base, "Documents/Isoplexis/Objects/cd8_base.rds")

DimPlot(iso_cd4, label=T)
DimPlot(iso_cd4, group.by="Sample")

iso_cd4_bm <- subset(iso_cd4, subset = Sample %in% "BM")
iso_cd4_pb <- subset(iso_cd4, subset = Sample %in% "PB")

iso_cd8_bm <- subset(iso_cd8, subset = Sample %in% "BM")
iso_cd8_pb <- subset(iso_cd8, subset = Sample %in% "PB")

iso_cd4_bm <- analyze_seurat(iso_cd4_bm)
iso_cd4_pb <- analyze_seurat(iso_cd4_pb)
iso_cd8_bm <- analyze_seurat(iso_cd8_bm)
iso_cd8_pb <- analyze_seurat(iso_cd8_pb)

########## CD4 ##########

clustree <- clustree(iso_cd4) + labs(title="Clustree")
plot_umap <- function(res){
  res_name = paste0("PSI_snn_res.", res)
  umap <- DimPlot(iso_cd4, group.by = res_name,label=T, raster=F) + 
    labs(title= paste0("UMAP for ", res_name))
  return(umap)
}
umap_plots <- map(seq(0,1,0.1), plot_umap)
pdf("Documents/Isoplexis/UMAPs/cd4.pdf",width = 8,height=7)
print(ggarrange(plotlist=umap_plots, ncol=1,nrow=1))
print(ggarrange(clustree))
dev.off()

DimPlot(iso_cd4, group.by="Time", cols=c("seagreen", "pink"), raster=F, shuffle=T) +
  NoAxes() + ggtitle (NULL) +
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 14),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"))+NoAxes() +
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("Documents/Isoplexis/UMAPs/cd4_time_new2.pdf",height=5,width=6,unit="in")

DimPlot(iso_cd4, group.by="Sample", cols=c("#B266FF","#FF9933"), raster=F, shuffle=F) +
  NoAxes()+ggtitle(NULL)+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 14),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"))+NoAxes()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("Documents/Isoplexis/UMAPs/cd4_sample_new2.pdf",height=5,width=6,unit="in")

DimPlot(iso_cd4, group.by="Sample", cols=c("#B266FF","#FF9933"), raster=F, shuffle=F, split.by = "Sample") +
  NoAxes()+ggtitle(NULL)+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 14),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"))+NoAxes()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("Documents/Isoplexis/UMAPs/cd4_sample_split.pdf",height=5,width=9,unit="in")

DimPlot(iso_cd4, group.by="Resp", cols=c("brown4", "steelblue1"), raster=F, shuffle=T) +
  NoAxes() + ggtitle (NULL) +
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 14),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"))+NoAxes() +
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("Documents/Isoplexis/UMAPs/cd4_resp.pdf",height=5,width=6,unit="in")


DimPlot(iso_cd4, group.by="PSI_snn_res.0.3",raster=F, shuffle=T, label = T, cols=ColAssign(unique(iso_cd4$PSI_snn_res.0.3)))+
  NoAxes()+ggtitle(NULL)+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 14),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"))+
  NoAxes()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("Documents/Isoplexis/UMAPs/cd4_res0.3_new.pdf",height=5,width=5,unit="in")


pdf("Documents/Isoplexis/ClusterMembership/cm_cd4_res0.3_time.pdf", width=5, height=1.5)
print(cluster_membership(iso_cd4$PSI_snn_res.0.3, iso_cd4$Time))
dev.off()

pdf("Documents/Isoplexis/ClusterMembership/cm_cd4_res0.3_sample.pdf", width=5, height=1.5)
print(cluster_membership(iso_cd4$PSI_snn_res.0.3, iso_cd4$Sample))
dev.off()

pdf("Documents/Isoplexis/ClusterMembership/cm_cd4_res0.3_Resp.pdf", width=5, height=1.5)
print(cluster_membership(iso_cd4$PSI_snn_res.0.3, iso_cd4$Resp))
dev.off()

cd4_counts <- GetAssayData(object = iso_cd4, slot = "counts")
counts_avg <- groupMeans(cd4_counts, groups = iso_cd4$PSI_snn_res.0.3)
counts_avg <- t(scale(t(counts_avg)))
time_cluster <- FetchData(iso_cd4, vars=c('PSI_snn_res.0.3', 'Time'))
time_cluster <- time_cluster %>% group_by(PSI_snn_res.0.3, Time) %>% dplyr::summarize(n=dplyr::n()) %>% 
  pivot_wider(names_from = Time, values_from = n) %>% data.frame() %>% mutate_all(~replace(., is.na(.), 0)) 
rownames(time_cluster) <- time_cluster$PSI_snn_res.0.3
time_cluster <- time_cluster[,-1]
time_cluster <- time_cluster/rowSums(time_cluster)
colnames(time_cluster) <- c("Baseline", "Two Month")
time_cluster <- time_cluster[colnames(counts_avg),]

sample_cluster <- FetchData(iso_cd4, vars=c('PSI_snn_res.0.3', 'Sample'))
sample_cluster <- sample_cluster %>% group_by(PSI_snn_res.0.3, Sample) %>% dplyr::summarize(n=dplyr::n()) %>% 
  pivot_wider(names_from = Sample, values_from = n) %>% data.frame() %>% mutate_all(~replace(., is.na(.), 0)) 
rownames(sample_cluster) <- sample_cluster$PSI_snn_res.0.3
sample_cluster <- sample_cluster[,-1]
sample_cluster <- sample_cluster/rowSums(sample_cluster)
sample_cluster <- sample_cluster[colnames(counts_avg),]
ha <- HeatmapAnnotation(Time = anno_barplot(time_cluster, gp = gpar(fill=c("seagreen", "pink"))),
                        Sample = anno_barplot(sample_cluster, gp=gpar(fill=c("#B266FF","#FF9933"))))


lgd_anno <- Legend(labels = c("Baseline", "Two Month"), legend_gp = gpar(fill = c("seagreen", "pink")),title = "Time Point")
lgd_anno2 <- Legend(labels = c("BM", "PB"), legend_gp = gpar(fill = c("#B266FF","#FF9933")),title = "Sample Location")

pdf("~/Documents/Isoplexis/Heatmaps/cd4_heatmap_res0.3_clust.pdf", height=8, width=7)
ht <- Heatmap(counts_avg,  width = ncol(counts_avg)*unit(5, "mm"), height = nrow(counts_avg)*unit(5, "mm"),
              heatmap_legend_param = list(title="Z-score"), cluster_columns = TRUE, top_annotation = ha,
              right_annotation = rightanno, show_column_dend = FALSE)
draw(ht, annotation_legend_list = list(lgd_anno, lgd_anno2))
dev.off()


counts_avg <- groupMeans(cd4_counts, groups = iso_cd4$PSI_snn_res.0.3)
counts_avg <- t(scale(t(counts_avg)))
time_cluster <- FetchData(iso_cd4, vars=c('PSI_snn_res.0.3', 'Time'))
time_cluster <- time_cluster %>% group_by(PSI_snn_res.0.3, Time) %>% dplyr::summarize(n=dplyr::n()) %>% 
  pivot_wider(names_from = Time, values_from = n) %>% data.frame() %>% mutate_all(~replace(., is.na(.), 0)) 
rownames(time_cluster) <- time_cluster$PSI_snn_res.0.3
time_cluster <- time_cluster[,-1]
time_cluster <- time_cluster/rowSums(time_cluster)
colnames(time_cluster) <- c("2Month", "Baseline")
time_cluster <- time_cluster[colnames(counts_avg),]

sample_cluster <- FetchData(iso_cd4, vars=c('PSI_snn_res.0.3', 'Sample'))
sample_cluster <- sample_cluster %>% group_by(PSI_snn_res.0.3, Sample) %>% dplyr::summarize(n=dplyr::n()) %>% 
  pivot_wider(names_from = Sample, values_from = n) %>% data.frame() %>% mutate_all(~replace(., is.na(.), 0)) 
rownames(sample_cluster) <- sample_cluster$PSI_snn_res.0.3
sample_cluster <- sample_cluster[,-1]
sample_cluster <- sample_cluster/rowSums(sample_cluster)
sample_cluster <- sample_cluster[colnames(counts_avg),]
ha <- HeatmapAnnotation(Time = anno_barplot(time_cluster, gp = gpar(fill=c("pink","seagreen"))),
                        Sample = anno_barplot(sample_cluster, gp=gpar(fill=c("#B266FF","#FF9933"))))

pdf("Documents/Isoplexis/Heatmaps/cd4_heatmap_res0.3.pdf", height=8, width=7)
ht <- Heatmap(counts_avg,  width = ncol(counts_avg)*unit(5, "mm"), height = nrow(counts_avg)*unit(5, "mm"),
              heatmap_legend_param = list(title="Z-score"), cluster_columns = FALSE, top_annotation = ha,
              right_annotation = rightanno)
draw(ht, annotation_legend_list = list(lgd_anno, lgd_anno2))
dev.off()

########## CD8 ##########

clustree <- clustree(iso_cd8) + labs(title="Clustree")
plot_umap <- function(res){
  res_name = paste0("PSI_snn_res.", res)
  umap <- DimPlot(iso_cd8, group.by = res_name,label=T, raster=F) + 
    labs(title= paste0("UMAP for ", res_name))
  return(umap)
}
umap_plots <- map(seq(0,1,0.1), plot_umap)
pdf("Documents/Isoplexis/UMAPs/cd8.pdf",width = 8,height=7)
print(ggarrange(plotlist=umap_plots, ncol=1,nrow=1))
print(ggarrange(clustree))
dev.off()

DimPlot(iso_cd8, group.by="Time", cols=c("seagreen", "pink"), raster=F, shuffle=T) +
  NoAxes() + ggtitle (NULL) +
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 14),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"))+NoAxes() +
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("Documents/Isoplexis/UMAPs/cd8_time_new2.pdf",height=5,width=6,unit="in")

DimPlot(iso_cd8, group.by="Sample", cols=c("#B266FF","#FF9933"), raster=F, shuffle=F) +
  NoAxes()+ggtitle(NULL)+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 14),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"))+NoAxes()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("Documents/Isoplexis/UMAPs/cd8_sample_new2.pdf",height=5,width=6,unit="in")

DimPlot(iso_cd8, group.by="Sample", cols=c("#B266FF","#FF9933"), raster=F, shuffle=F, split.by = "Sample") +
  NoAxes()+ggtitle(NULL)+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 14),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"))+NoAxes()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("Documents/Isoplexis/UMAPs/cd8_sample_split.pdf",height=5,width=9,unit="in")

DimPlot(iso_cd8, group.by="PSI_snn_res.0.4",raster=F, shuffle=T, label = T, cols=ColAssign(unique(iso_cd8$PSI_snn_res.0.4)))+
  NoAxes()+ggtitle(NULL)+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 14),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"))+
  NoAxes()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("Documents/Isoplexis/UMAPs/cd8_res0.4_new.pdf",height=5,width=5,unit="in")

DimPlot(iso_cd8, group.by="Resp", cols=c("brown4", "steelblue1"), raster=F, shuffle=T) +
  NoAxes() + ggtitle (NULL) +
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 14),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"))+NoAxes() +
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("Documents/Isoplexis/UMAPs/cd8_resp.pdf",height=5,width=6,unit="in")


pdf("Documents/Isoplexis/ClusterMembership/cm_cd8_res0.4_time.pdf", width=5, height=1.5)
print(cluster_membership(iso_cd8$PSI_snn_res.0.4, iso_cd8$Time))
dev.off()

pdf("Documents/Isoplexis/ClusterMembership/cm_cd8_res0.4_sample.pdf", width=5, height=1.5)
print(cluster_membership(iso_cd8$PSI_snn_res.0.4, iso_cd8$Sample))
dev.off()

pdf("Documents/Isoplexis/ClusterMembership/cm_cd8_res0.4_Resp.pdf", width=5, height=1.5)
print(cluster_membership(iso_cd8$PSI_snn_res.0.4, iso_cd8$Resp))
dev.off()


cd8_counts <- GetAssayData(object = iso_cd8, slot = "counts")
counts_avg <- groupMeans(cd8_counts, groups = iso_cd8$PSI_snn_res.0.4)
counts_avg <- t(scale(t(counts_avg)))
time_cluster <- FetchData(iso_cd8, vars=c('PSI_snn_res.0.4', 'Time'))
time_cluster <- time_cluster %>% group_by(PSI_snn_res.0.4, Time) %>% dplyr::summarize(n=dplyr::n()) %>% 
  pivot_wider(names_from = Time, values_from = n) %>% data.frame() %>% mutate_all(~replace(., is.na(.), 0)) 
rownames(time_cluster) <- time_cluster$PSI_snn_res.0.4
time_cluster <- time_cluster[,-1]
time_cluster <- time_cluster/rowSums(time_cluster)
colnames(time_cluster) <- c("Baseline", "Two Month")
time_cluster <- time_cluster[colnames(counts_avg),]

sample_cluster <- FetchData(iso_cd8, vars=c('PSI_snn_res.0.4', 'Sample'))
sample_cluster <- sample_cluster %>% group_by(PSI_snn_res.0.4, Sample) %>% dplyr::summarize(n=dplyr::n()) %>% 
  pivot_wider(names_from = Sample, values_from = n) %>% data.frame() %>% mutate_all(~replace(., is.na(.), 0)) 
rownames(sample_cluster) <- sample_cluster$PSI_snn_res.0.4
sample_cluster <- sample_cluster[,-1]
sample_cluster <- sample_cluster/rowSums(sample_cluster)
sample_cluster <- sample_cluster[colnames(counts_avg),]
ha <- HeatmapAnnotation(Time = anno_barplot(time_cluster, gp = gpar(fill=c("seagreen", "pink"))),
                        Sample = anno_barplot(sample_cluster, gp=gpar(fill=c("#B266FF","#FF9933"))))


pdf("~/Documents/Isoplexis/Heatmaps/cd8_heatmap_res0.4_clust.pdf", height=8, width=7)
ht <- Heatmap(counts_avg,  width = ncol(counts_avg)*unit(5, "mm"), height = nrow(counts_avg)*unit(5, "mm"),
              heatmap_legend_param = list(title="Z-score"), cluster_columns = TRUE, top_annotation = ha,
              right_annotation = rightanno, show_column_dend = FALSE)
draw(ht, annotation_legend_list = list(lgd_anno, lgd_anno2))
dev.off()

######### miloR ##########

library(miloR)

cd8_milo <- as.SingleCellExperiment(iso_cd8)
cd8_milo <- Milo(cd8_milo)
cd8_milo <- buildGraph(cd8_milo, d=10, k=20, reduced.dim = "HARMONY")
cd8_milo <- makeNhoods(cd8_milo, prop = 0.05, d=10, k=20, refined = TRUE, reduced_dims = "HARMONY")
plotNhoodSizeHist(cd8_milo)
cd8_milo <- countCells(cd8_milo, meta.data = data.frame(colData(cd8_milo)), sample="PatientID")
cd8_milo <- calcNhoodDistance(cd8_milo, d=10, reduced.dim = "HARMONY")

cd8_design <- data.frame(colData(cd8_milo))[,c("PatientID", "Time")]
cd8_design$Time <- as.factor(cd8_design$Time)
cd8_design <- distinct(cd8_design)
rownames(cd8_design) <- c(cd8_design$PatientID)

da_results <- testNhoods(cd8_milo, design = ~ Time, design.df = cd8_design)
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)

cd8_milo <- buildNhoodGraph(cd8_milo)
milo <- plotNhoodGraphDA(cd8_milo, da_results, alpha=1,size_range=c(1,5)) + guides(edge_width="none")


da_results <- annotateNhoods(cd8_milo, da_results, coldata_col = "Sample")
da_results <- annotateNhoods(cd8_milo, da_results, coldata_col = "Resp")
da_results <- annotateNhoods(cd8_milo, da_results, coldata_col = "PSI_snn_res.0.4")
#ggplot(da_results, aes(type_fraction)) + geom_histogram(bins=50)
da_results$Sample <- ifelse(da_results$Sample_fraction < 0.7, "Mixed", da_results$Sample)
da_results$Resp <- ifelse(da_results$Resp_fraction < 0.7, "Mixed", da_results$Resp)
#da_results$PSI_snn_res.0.4 <- ifelse(da_results$PSI_snn_res.0.4_fraction < 0.7, "Mixed", da_results$PSI_snn_res.0.4)
table(da_results$PSI_snn_res.0.4, da_results$logFC > 0)

da_results <- annotateNhoods(cd8_milo, da_results, coldata_col = "Time")
da_results$Time <- ifelse(da_results$Time_fraction < 0.7, "Mixed", da_results$Time)

#milo <- milo + guides(size=guide_legend(title="Neighborhood Size"))

pdf("Documents/Isoplexis/Milo/cd8_milo_umap_new.pdf", height=5, width=7)
milo
dev.off()

pdf("Documents/Isoplexis/Milo/cd8_sample_beeswarm.pdf", height=2, width=3.5)
plotDAbeeswarm(da_results, group.by = "Sample") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

pdf("Documents/Isoplexis/Milo/cd8_resp_beeswarm_nomixed.pdf", height=3, width=4)
plotDAbeeswarm(da_results, group.by = "Resp") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

da_results$Resp <- ifelse(da_results$Resp_fraction < 0.7, "Mixed", da_results$Resp)
pdf("Documents/Isoplexis/Milo/cd8_resp_beeswarm_mixed.pdf", height=3, width=4)
plotDAbeeswarm(da_results, group.by = "Resp") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

pdf("Documents/Isoplexis/Milo/cd8_time_beeswarm.pdf", height=2, width=3.5)
plotDAbeeswarm(da_results, group.by = "Time") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

pdf("Documents/Isoplexis/Milo/cd8_res0.4_beeswarm.pdf", height=3, width=7)
plotDAbeeswarm(da_results, group.by = "PSI_snn_res.0.4") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

cd4_milo <- as.SingleCellExperiment(iso_cd4)
cd4_milo <- Milo(cd4_milo)
cd4_milo <- buildGraph(cd4_milo, d=10, k=20, reduced.dim = "HARMONY")
cd4_milo <- makeNhoods(cd4_milo, prop = 0.05, d=10, k=20, refined = TRUE, reduced_dims = "HARMONY")
cd4_milo <- countCells(cd4_milo, meta.data = data.frame(colData(cd4_milo)), sample="PatientID")
cd4_milo <- calcNhoodDistance(cd4_milo, d=10, reduced.dim = "HARMONY")
plotNhoodSizeHist(cd4_milo)

cd4_design <- data.frame(colData(cd4_milo))[,c("PatientID", "Time")]
cd4_design$Time <- as.factor(cd4_design$Time)
cd4_design <- distinct(cd4_design)
rownames(cd4_design) <- c(cd4_design$PatientID)

da_results_cd4 <- testNhoods(cd4_milo, design = ~ Time, design.df = cd4_design)
da_results_cd4 <- annotateNhoods(cd4_milo, da_results_cd4, coldata_col = "Sample")
da_results_cd4 <- annotateNhoods(cd4_milo, da_results_cd4, coldata_col = "Resp")
da_results_cd4 <- annotateNhoods(cd4_milo, da_results_cd4, coldata_col = "Time")
da_results_cd4 <- annotateNhoods(cd4_milo, da_results_cd4, coldata_col = "PSI_snn_res.0.3")
#ggplot(da_results, aes(type_fraction)) + geom_histogram(bins=50)
da_results_cd4$Sample <- ifelse(da_results_cd4$Sample_fraction < 0.7, "Mixed", da_results_cd4$Sample)
#da_results_cd4$Resp <- ifelse(da_results_cd4$Resp_fraction < 0.7, "Mixed", da_results_cd4$Resp)
#da_results_cd4$PSI_snn_res.0.3 <- ifelse(da_results_cd4$PSI_snn_res.0.3_fraction < 0.7, "Mixed", da_results_cd4$PSI_snn_res.0.3)
da_results_cd4$Time <- ifelse(da_results_cd4$Time_fraction < 0.7, "Mixed", da_results_cd4$Time)

cd4_milo <- buildNhoodGraph(cd4_milo)
milo_cd4 <- plotNhoodGraphDA(cd4_milo, da_results_cd4, alpha=1,size_range=c(1,5)) + guides(edge_width="none")
#milo <- milo + guides(size=guide_legend(title="Neighborhood Size"))

pdf("Documents/Isoplexis/Milo/cd4_milo_umap_new.pdf", height=5, width=7)
milo_cd4
dev.off()

pdf("Documents/Isoplexis/Milo/cd4_sample_beeswarm.pdf", height=2, width=3.5)
plotDAbeeswarm(da_results_cd4, group.by = "Sample") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

pdf("Documents/Isoplexis/Milo/cd4_resp_beeswarm_nomixed.pdf", height=3, width=4)
plotDAbeeswarm(da_results_cd4, group.by = "Resp") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

da_results_cd4$Resp <- ifelse(da_results_cd4$Resp_fraction < 0.7, "Mixed", da_results_cd4$Resp)
pdf("Documents/Isoplexis/Milo/cd4_resp_beeswarm_mixed.pdf", height=3, width=4)
plotDAbeeswarm(da_results_cd4, group.by = "Resp") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

pdf("Documents/Isoplexis/Milo/cd4_time_beeswarm.pdf", height=2, width=3.5)
plotDAbeeswarm(da_results_cd4, group.by = "Time") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()
pdf("Documents/Isoplexis/Milo/cd4_res0.3_beeswarm.pdf", height=3, width=7)
plotDAbeeswarm(da_results_cd4, group.by = "PSI_snn_res.0.3") + xlab(NULL) + ylab(NULL) +  theme(axis.text = element_text(size=10), axis.title = element_blank())
dev.off()

########## Testing groups by Resp ##########

a <- c("IL1b", "IL13", "IL17F", "RANTES", "sCD40L", "CCL11", "IP10", "TNFb", "IL21", "MIP1a", "GranzymeB", "sCD137")
b <- c("TGFb1", "MCP1", "MCP4", "IL7")
c <- c("TNFa", "IFNg", "IL8", "IL10", "MIP1b")
d <- c("IL4", "IL22", "IL9", "IL2", "IL15", "GMCSF",  "IL12","Perforin", 
       "IL6", "IL17A")
e <- c("IL5")

iso_cd4$a <- rowSums(FetchData(iso_cd4, a), na.rm = TRUE)
iso_cd4$b <- rowSums(FetchData(iso_cd4, b), na.rm = TRUE)
iso_cd4$c <- rowSums(FetchData(iso_cd4, c), na.rm = TRUE)
iso_cd4$d <- rowSums(FetchData(iso_cd4, d), na.rm = TRUE)
iso_cd4$e <- rowSums(FetchData(iso_cd4, e), na.rm = TRUE)

cd4_base <- subset(iso_cd4, subset=Time %in% "Baseline")

cd4_new <- FetchData(cd4_base, c("a","b", "c", "d", "e", "Resp", "Patient", "Sample"))
cd4_new_sum <- cd4_new %>% 
  group_by(Patient,Resp, Sample) %>% 
  summarise(a = mean(a), b = mean(b), c=mean(c), d=mean(d), e=mean(e))
cd4_sum_sub <- cd4_new_sum  %>%
  pivot_longer(cols=a:e, names_to="Group", values_to = "PSI")
cd4_sum_sub$Resp <- as.factor(cd4_sum_sub$Resp)
stat.test.cd4 <- cd4_sum_sub %>%
  group_by(Group,Sample) %>%
  pairwise_wilcox_test(PSI~ Resp, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x="Resp", fun="max", step.increase = 0.01)

pdf("Documents/Isoplexis/BoxPlots/baseline_resp_cd4.pdf", width=8, height=8)
ggboxplot(cd4_sum_sub, "Resp", "PSI", fill="Resp", 
          width=0.75, lwd=0.25, palette = c("brown4", "steelblue3"), bxp.errorbar = TRUE) + 
  facet_grid2(Sample~Group, scales="free_y", independent = "y") +
  theme_test() +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle=90), legend.position = "none") +
  stat_pvalue_manual(stat.test.cd4, label="p.adj", tip.length = 0.0)
dev.off()

#For CD8 cells
a <- c("IL12", "IL15",  "GMCSF", "TGFb1")
b <- c("IL1b", "IL7", "IL13", "IL17F", "IL21", "RANTES", "sCD40L", "CCL11",
       "IP10", "MCP1", "MCP4")
c <- c("IL17A", "IL2", "IL6",  "IL10", "Perforin")
d <- c("TNFb", "IL4", "IL5")
e <- c("TNFa", "IFNg", "MIP1b", "IL9", "IL8", "IL22")
f <- c("MIP1a", "sCD137", "GranzymeB")

iso_cd8$a <- rowSums(FetchData(iso_cd8, a), na.rm = TRUE)
iso_cd8$b <- rowSums(FetchData(iso_cd8, b), na.rm = TRUE)
iso_cd8$c <- rowSums(FetchData(iso_cd8, c), na.rm = TRUE)
iso_cd8$d <- rowSums(FetchData(iso_cd8, d), na.rm = TRUE)
iso_cd8$e <- rowSums(FetchData(iso_cd8, e), na.rm = TRUE)
iso_cd8$f <- rowSums(FetchData(iso_cd8, f), na.rm = TRUE)

cd8_base <- subset(iso_cd8, subset=Time %in% "Baseline")

cd8_new <- FetchData(cd8_base, c("a","b", "c", "d", "e", "f", "Resp", "Patient", "Sample"))
cd8_new_sum <- cd8_new %>% 
  group_by(Patient,Resp, Sample) %>% 
  summarise(a = mean(a), b = mean(b), c=mean(c), d=mean(d), e=mean(e), f = mean(f))
cd8_sum_sub <- cd8_new_sum  %>%
  pivot_longer(cols=a:f, names_to="Group", values_to = "PSI")
cd8_sum_sub$Resp <- as.factor(cd8_sum_sub$Resp)
stat.test.cd8 <- cd8_sum_sub %>%
  group_by(Group,Sample) %>%
  pairwise_wilcox_test(PSI~ Resp, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x="Resp", fun="max", step.increase = 0.01)


pdf("Documents/Isoplexis/BoxPlots/baseline_resp_cd8.pdf", width=8, height=8)
ggboxplot(cd8_sum_sub, "Resp", "PSI", fill="Resp", 
          width=0.75, lwd=0.25, palette = c("brown4", "steelblue3"), bxp.errorbar = TRUE) + 
  facet_grid2(Sample~Group, scales="free_y", independent = "y") +
  theme_test() +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle=90), legend.position = "none") +
  stat_pvalue_manual(stat.test.cd8, label="p.adj", tip.length = 0.0)
dev.off()