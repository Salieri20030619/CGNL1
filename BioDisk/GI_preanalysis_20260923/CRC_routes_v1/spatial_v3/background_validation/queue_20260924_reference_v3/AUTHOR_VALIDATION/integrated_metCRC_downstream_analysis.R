library(nicheDE)
library(enrichR)
#read in object
Merged_NDE = readRDS('merged_NDE_22261_fisher.rds')

#Figure 4 Pathway enrichment analysis for fibroblasts and tumor
#get fibroblast tumor niche genes
fibro_tum_pos = get_niche_DE_genes(Merged_NDE,'I',index='stromal',niche = 'tumor_epithelial',positive = T,alpha = 0.05)

#Figure 4B) perform pathway enrichment analysis
fibro_tum_processes = enrichr(fibro_tum_pos[,1],databases = 'Reactome_2016')
View(fibro_tum_processes$Reactome_2016)


#Figure 4E) Niche LR analysis for fibroblasts and tumor
data("niche_net_ligand_target_matrix")
data("ramilowski_ligand_receptor_list")
fibro_tumor_LR = niche_LR_spot(Merged_NDE,'tumor_epithelial','stromal',niche_net_ligand_target_matrix,
                                   ramilowski_ligand_receptor_list,25,50,0.05)

#Figure 5A) Pathway enrichment analysis for macrophages and tumor
#get macro,tumor+ niche-DE genes
macro_tum_pos = get_niche_DE_genes(Merged_NDE,'I',index='myeloid',niche = 'tumor_epithelial',positive = T,alpha = 0.05)
#perform pathway enrichment analysis
macro_tum_processes = enrichr(macro_tum_pos[,1],databases = 'Reactome_2016')
View(macro_tum_processes$Reactome_2016)


#Figure 5B) LR analysis for macrophages and tumor
data("niche_net_ligand_target_matrix")
data("ramilowski_ligand_receptor_list")
macro_tumor_LR = niche_LR_spot(Merged_NDE,'tumor_epithelial','myeloid',niche_net_ligand_target_matrix,
                               ramilowski_ligand_receptor_list,25,50,0.05,truncation_value = 5)

#Figure 5C) TAM markers and normal macrophage markers
tum_marker = niche_DE_markers(Merged_NDE,'myeloid','tumor_epithelial','hepatocytes+cholangiocytes',0.05)
norm_marker = niche_DE_markers(Merged_NDE,'myeloid','hepatocytes+cholangiocytes','tumor_epithelial',0.05)





#Figure 5G) pathway analysis of TAM markers
TAM_processes = enrichr(tum_marker[,1],databases = 'Reactome_2016')
View(TAM_processes$Reactome_2016)







# Figure 5F) Umap module score for TAM marker genes

#Read in data
sobj = readRDS('scrna_seq_macro_fine.rds')
#read in SAVER results
A = readRDS('saver.rds')

sobj[["Saver"]] <- Seurat::CreateAssayObject(counts = A$estimate)
Seurat::DefaultAssay(sobj) <- "Saver"

#get kuppfer and TAM markers
kup_markers = norm_marker[,1]
tum_markers = tum_marker[,1]

#type seurat pipeline
sobj = FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sobj)
sobj <- ScaleData(sobj, features = all.genes,do.scale = F, do.center = F)
sobj = RunPCA(sobj)

#clustering
sobj <- FindNeighbors(sobj, dims = 1:50)
#sobj <- FindClusters(sobj, resolution = 0.5)
sobj <- RunUMAP(sobj, dims = 1:50)
DimPlot(sobj, reduction = "umap")

#get module score for kupffer markers
sobj <- AddModuleScore(sobj,
                       features = list(kup_markers),
                       name="kup_markers")
#get module score for TAM markers
sobj <- AddModuleScore(sobj,
                       features = list(tum_markers),
                       name="tum_markers")



#kupffer marker heatmap
gene = 'kup_markers1'
FeaturePlot(sobj,features = gene,cols = c("white", "blue"))+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
                                                                  axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  xlab('')+ylab('')+ggtitle('')



#TAM marker heatmap
gene = 'tum_markers1'
FeaturePlot(sobj,features = gene,cols = c("white", "blue"))+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
                                                                  axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  xlab('')+ylab('')+ggtitle('')





