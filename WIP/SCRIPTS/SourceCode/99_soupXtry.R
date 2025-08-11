# Make the soup
toc = Seurat::Read10X(file.path( '/igm/projects/250718_GSL-EM-4483/FF061/outs/per_sample_outs/SOM301/count/sample_filtered_feature_bc_matrix'))

tod = Seurat::Read10X(file.path("/igm/projects/250718_GSL-EM-4483/FF061/outs/per_sample_outs/SOM301/count/sample_raw_feature_bc_matrix"))

# Cellranger automatically filtering out deprecated probes from their FRP protocol. Based on the description of Cellranger multi's outputs, it seems like the raw matrix includes these probes while the filtered matrix does not. Perhaps these are responsible for the discrepancy? (from github issue 148 reply). So filter out those excluded and only get ones in toc.

# Get toc genes
filt_genes <- rownames(toc)

# Subset tod to keep only the genes in toc
tod_sub <- tod[rownames(tod) %in% filt_genes, ]

obj  <- CreateSeuratObject(counts = toc)

sc = SoupChannel(tod_sub, toc)

obj <- NormalizeData(obj, verbose = F)
obj <- FindVariableFeatures(obj, verbose = F)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, verbose = F)
obj <- RunUMAP(obj, dims = 1:30, verbose = F)
obj <- FindNeighbors(obj, dims = 1:30, verbose = F)
obj <- FindClusters(obj, resolution = 0.1, verbose = T)
obj <- SCTransform(obj)

Seurat::DefaultAssay(obj) <- "SCT"
Seurat::DefaultAssay(ref_obj) <- "SCT"

# run mouse annotation
obj <- ref_mousebrain_annotation(ref_obj,
                                 obj,
                                 save = F)

# Ref cell type prediction will be stored in seurat_obj$predicted.id.
Idents(obj) <- obj$predicted.id

# Visualize
DimPlot(obj, reduction = "umap", label = T, repel = T) + NoLegend()

# get names and shorten them
ln <- unique(obj$predicted.id)

sn <- c("Ex_Cortex","Cajal_Retzius",  "Interneuron", "Ex_CA3","Ex_CA1",  "Astrocytes",  "Interneuron", "Granule_DG", "Vascular",  "Astrocytes",
        "NPC", "Mature_OL", "MFOL", "Interneuron", "Microglia",  "Interneuron","MGE","Microglia","OPC", "Trilaminar","Pericytes",
        "Ependymal","Granule_Blast_DG","NFOL","Vascular", "COP", "Interneuron", "CGE", "CGE",  "Vascular", "Interneuron", "Macrophages","RGL_DG", "Vascular")

# Mapping sn-ln
name_map <- setNames(sn, ln)

# Function to replace long names to short names
shortened <- function(name) {ifelse(name %in% names(name_map), name_map[name], name)}

# Apply function
obj$ref_short <- shortened(obj$predicted.id)

# Check UMAP
Idents(obj) <- "ref_short"
DimPlot(obj, reduction = "umap", label = T, repel = T) + NoLegend()
oligsub <- subset(obj, integrated$ref_short %in% c("OPC", "COP", "NFOL", "MFOL", "Mature_OL"))

meta <- obj@meta.data
umap <- obj@reductions$umap@cell.embeddings
soup.channel  <- setClusters(sc, setNames(meta$ref_short, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
head(meta)

soup.channel  <- autoEstCont(soup.channel)
head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 20)
plotMarkerDistribution(soup.channel)

# Check the results
summary(soup.channel$metaData$rho)
hist(soup.channel$metaData$rho, breaks = 20, main = "Contamination Distribution")

# Get the corrected count matrix
out <- adjustCounts(soup.channel)

# Create Seurat object with corrected counts
sobj_corrected <- CreateSeuratObject(counts = out, project = "SOM301_SoupX")
sobj_corrected    <- NormalizeData(sobj_corrected, verbose = F)
sobj_corrected    <- FindVariableFeatures(sobj_corrected, verbose = F)
sobj_corrected    <- ScaleData(sobj_corrected, verbose = F)
sobj_corrected    <- RunPCA(sobj_corrected, verbose = F)
sobj_corrected    <- RunUMAP(sobj_corrected, dims = 1:30, verbose = F)
sobj_corrected    <- FindNeighbors(sobj_corrected, dims = 1:30, verbose = F)
sobj_corrected    <- FindClusters(sobj_corrected, resolution = 0.1, verbose = T)
sobj_corrected <- SCTransform(sobj_corrected)

# run mouse annotation
sobj_corrected <- ref_mousebrain_annotation(ref_obj,
                                            sobj_corrected,
                                            save = F)

# Ref cell type prediction will be stored in seurat_obj$predicted.id.
Idents(sobj_corrected) <- sobj_corrected$predicted.id

# Visualize
DimPlot(sobj_corrected, reduction = "umap", label = T, repel = T) + NoLegend()

# get names and shorten them
ln <- unique(sobj_corrected$predicted.id)

sn <- c("Ex_Cortex","Cajal_Retzius",  "Interneuron", "Ex_CA3","Ex_CA1",  "Astrocytes", "Granule_DG", "Vascular",  "Astrocytes",
        "NPC", "Mature_OL", "MFOL", "Interneuron", "OPC", "Interneuron", "Microglia",  "Interneuron","MGE","Microglia", "Trilaminar","Pericytes",
        "Ependymal","Granule_Blast_DG","NFOL","Vascular", "COP", "Interneuron", "CGE", "CGE",  "Vascular", "Interneuron", "Macrophages","RGL_DG", "Vascular")

# Mapping sn-ln
name_map <- setNames(sn, ln)

# Function to replace long names to short names
shortened <- function(name) {ifelse(name %in% names(name_map), name_map[name], name)}

# Apply function
sobj_corrected$ref_short <- shortened(sobj_corrected$predicted.id)

# Check the highly contaminating genes:
Idents(obj) <- "ref_short"
Idents(sobj_corrected) <- "ref_short"

VlnPlot2(obj, features = c("Slc1a2", "Plp1"), show.mean = F, pt = F, box = F) |
    VlnPlot2(sobj_corrected, features = c("Slc1a2", "Plp1"), show.mean = F, pt = F, box = F)

VlnPlot2(obj, features = "Slc35a2",show.mean = F, pt = F, box = F) |
    VlnPlot2(sobj_corrected, features = "Slc35a2",show.mean = F, pt = F, box = F)