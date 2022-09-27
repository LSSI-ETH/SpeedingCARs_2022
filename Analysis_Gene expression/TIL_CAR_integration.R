library(colorspace)
library(Seurat)
library(mclust)
library(dplyr)
library(sctransform)
library(DESeq2)
library(cowplot)
library(patchwork)
library(harmony)
library(Rcpp)
library(tuple)
library(tidyr)
library(tidyverse)
#library(ggtree)
library(stringdist)
library(ape)
library(data.table)
library(fastcluster)
library(ggplot2)
library(circlize)
library(pheatmap)
library(slingshot)
library(princurve)
library(SingleCellExperiment)
library(TSCAN)
library(tuple)
library(tidyr)
library(Matrix)
set.seed(1234)

#Read in CAR data:
car_all_cd8 <- subset(car_all, cells=WhichCells(car_all, expression = CD4==0))
car_all_cd8 <- subset(car_all_cd8, cells=cells_to_select)
car_all_cd8 <- AddMetaData(car_all_cd8, rownames(car_all_cd8@meta.data) %in% CIC_cells, "CIC")

car_all_cd8[["percent.mt"]] <- PercentageFeatureSet(object = car_all_cd8, pattern = "^MT-")
head(car_all_cd8@meta.data, 5)
VlnPlot(car_all_cd8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(object = car,feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(object = car,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))
#nFeature_RNA: genes per cell; nCount_RNA: molecules detected per cell
car_all_cd8 <- subset(x = car_all_cd8, subset = nFeature_RNA > 150 & nFeature_RNA < 6500 & percent.mt < 10)
saveRDS(car_all_cd8, "/Users/fbieberich/Downloads/car_final_20220204.rds")
seurat_object.list <- SplitObject(car_cd8, split.by = "orig.ident")
for (i in 1:length(seurat_object.list)) {
    seurat_object.list[[i]] <- SCTransform(seurat_object.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}
# Save variable features:
car.features <- SelectIntegrationFeatures(object.list = seurat_object.list, nfeatures = 3000)

seurat_object.merged <- merge(seurat_object.list[[1]], y = seurat_object.list[2:length(seurat_object.list)], project = "car", merge.data = TRUE)

VariableFeatures(seurat_object.merged) <- car.features
seurat_object.merged <- RunPCA(object = seurat_object.merged, assay = "SCT", npcs = 50)
DimPlot(seurat_object.merged, reduction = "pca")
car.h <- RunHarmony(object = seurat_object.merged,
                                    assay.use = "SCT",
                                   reduction = "pca",
                                    dims.use = 1:50,
                                    group.by.vars = "orig.ident",
                                    plot_convergence = TRUE)

# Read in TIL data from Liu et al. Temporal tracing of TILs
pd1.data <- readRDS("GSE179994_all.Tcell.rawCounts(1).rds.gz")

pt_selection <- c("P29.tr.1",
"P29.ut",
"P33.tr.1",
"P33.ut",
"P35.tr.1",
"P35.ut")

pd1.data_pt1 <- pd1.data[,str_extract(colnames(pd1.data), "^.{6}") %in% pt_selection[c(2,4,6)]]
pd1.data_pt2 <- pd1.data[,str_extract(colnames(pd1.data), "^.{8}") %in% pt_selection[c(1,3,5)]]
pd1.data_filt <- cbind(pd1.data_pt1,pd1.data_pt2)

pd1.seurat <- CreateSeuratObject(counts = pd1.data_filt, project = "pd1_liu", min.cells = 3, min.features = 200)
pd1.seurat <- AddMetaData(pd1.seurat, str_extract(rownames(pd1.seurat@meta.data), "^.{6}"), "orig.ident")
table(pd1.seurat$orig.ident)

# QC
pd1.seurat[["percent.mt"]] <- PercentageFeatureSet(object = pd1.seurat, pattern = "^MT-")
head(pd1.seurat@meta.data, 5)

VlnPlot(pd1.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#nFeature_RNA: genes per cell; nCount_RNA: molecules detected per cell
pd1.seurat <- subset(x = pd1.seurat, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 10)

patient_select <- c("P1.tr.1", "P1.tr.2","P1.tr.3", "P1.ut", "P13.tr.1", "P13.tr.2", "P13.ut",
                             "P19.tr", "P19.ut", "P30.tr", "P30.ut","P36.tr.1", "P38.tr.1")
#pd1.seurat <- AddMetaData(pd1.seurat, str_extract(rownames(pd1.seurat@meta.data), "^.{6}"), "orig.ident")
#table(str_extract(rownames(pd1.seurat@meta.data), "^.{8}"))
#pd1_sub.seurat <- subset(pd1.seurat, cells = rownames(pd1.seurat@meta.data)[
#  pd1.seurat$orig.ident %in% patient_select])

origin <- rep(NA, length(rownames(pd1.seurat@meta.data)))
for(i in 1:length(patient_select)){
origin[grep(patient_select[i], rownames(pd1.seurat@meta.data))] <- patient_select[i]
}
table(origin)

pd1.seurat <- AddMetaData(pd1.seurat, origin, "orig.ident")
pd1_sub.seurat <- subset(pd1.seurat, cells = rownames(pd1.seurat@meta.data)[pd1.seurat$orig.ident %in% patient_select])

# Combine TIL Liu et al data with CAR data:
pd1_subcd8.seurat <- subset(pd1_sub.seurat, cells=WhichCells(pd1_sub.seurat, expression = CD4==0))
saveRDS(pd1_subcd8.seurat, "/Users/fbieberich/Downloads/pd1_subcd8.seurat20220204.rds")
pd1_subcd8.seurat <- readRDS("/Users/fbieberich/Downloads/pd1_subcd8.seurat20220204.rds")

seurat_object.list <- SplitObject(pd1_subcd8.seurat, split.by = "orig.ident")
seurat_object.list2 <- SplitObject(car_all_cd8, split.by = "orig.ident")
seurat_object.list <- c(seurat_object.list, seurat_object.list2)

#CellCycle addition:
for (i in 1:length(seurat_object.list)) {
    seurat_object.list[[i]] <- NormalizeData(seurat_object.list[[i]])
}
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
for (i in 1:length(seurat_object.list)) {
    seurat_object.list[[i]] <- CellCycleScoring(seurat_object.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
}

for (i in 1:length(seurat_object.list)) {
    seurat_object.list[[i]] <- SCTransform(seurat_object.list[[i]], vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = FALSE)
}
# Save variable features:
lung.features <- SelectIntegrationFeatures(object.list = seurat_object.list, nfeatures = 3000)
# Merge data and run PCA:
reference.seurat <- merge(seurat_object.list[[1]], y = seurat_object.list[2:length(seurat_object.list)], project = "reference", merge.data = TRUE)
# Made on cluster
VariableFeatures(reference.seurat) <- lung.features
reference.seurat <- RunPCA(object = reference.seurat, assay = "SCT", npcs = 50)
DimPlot(reference.seurat, reduction = "pca")
reference.seurat <- RunHarmony(object = reference.seurat,
                                    assay.use = "SCT",
                                   reduction = "pca",
                                    dims.use = 1:50,
                                    group.by.vars = "orig.ident",
                                    plot_convergence = TRUE)
reference.seurat <- RunUMAP(object = reference.seurat, assay = "SCT", reduction = "harmony", return.model = TRUE, dims = 1:50)
reference.seurat <- FindNeighbors(object = reference.seurat, assay = "SCT", reduction = "harmony", dims = 1:50)
reference.seurat <- FindClusters(object = reference.seurat, resolution = .5)
DimPlot(reference.seurat, label = T)
table(reference.seurat$orig.ident)

ref_liupd1_pd1 <- reference.seurat
pd1_sub.seurat <- readRDS("/Users/fbieberich/Downloads/ref_liuPD1.rds")
table(pd1_sub.seurat$orig.ident)

ref_car_liupd1 <- reference.seurat

reference.seurat <- readRDS("/Users/fbieberich/Downloads/car_ref_liuPD1_CD8_20220127.rds")

reference.seurat <- ref_liupd1
rownames(ref_cd8@meta.data)
ref_cd4 <- subset(reference.seurat, cells=WhichCells(reference.seurat, expression = CD4>0))
ref_cd8 <- subset(reference.seurat, cells=WhichCells(reference.seurat, expression = CD4==0))

ref_cd8 <- RunUMAP(object = ref_cd8, assay = "SCT", reduction = "harmony", return.model = TRUE, dims = 1:50)
ref_cd8 <- FindNeighbors(object = ref_cd8, assay = "SCT", reduction = "harmony", dims = 1:50)
ref_cd8 <- FindClusters(object = ref_cd8, resolution = .4)
DimPlot(reference.seurat)
table(reference.seurat$orig.ident)
FeaturePlot(reference.seurat, features = c("SELL","IL7R", "MKI67",  "GZMK", "ENTPD1", "NR4A2"), order = T, label=T)
FeaturePlot(reference.seurat, features = c("PDCD1","GZMK", "NR4A2"), order = T, label=F, ncol = 3)
FeaturePlot(reference.seurat, features = c("ENTPD1","HAVCR2", "LAYN"), order = T, label=T, ncol = 3)

DimPlot(reference.seurat, cells.highlight = list(rownames(reference.seurat@meta.data)[reference.seurat$orig.ident == 'D3_hT'],
                                                 rownames(reference.seurat@meta.data)[reference.seurat$orig.ident == "D2_WT"],
                                                 rownames(reference.seurat@meta.data)[reference.seurat$CIC == 'TRUE']), cols.highlight = c('red', 'lightgreen', 'lightblue'))
