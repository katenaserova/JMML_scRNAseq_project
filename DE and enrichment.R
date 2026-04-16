


# Load JMML reference data


library(Seurat)
library(openxlsx)
library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)


# Load annotated JMML-reference datasets
sce_kras <- readH5AD("KRAS_annotated.h5ad")
sce_nras <- readH5AD("NRAS_annotated.h5ad")


# Convert JMML-references to Seurat
kras <- as.Seurat(
  sce_kras,
  counts = "counts",
  data = "X"
)

nras <- as.Seurat(
  sce_nras,
  counts = "counts",
  data = "X"
)


# Check of JMML-reference annotations
table(kras$celltype_l2_predicted)
table(kras$celltype_l1_predicted)
table(nras$cell_type)


# QC metrics of JMML-reference
VlnPlot(
  kras,
  features = c("n_genes_by_counts", "total_counts", "pct_counts_mt"),
  ncol = 3
)

VlnPlot(
  nras,
  features = c("n_genes_by_counts", "total_counts", "pct_counts_mt"),
  ncol = 3
)


# Define marker genes for comparison
markers_list <- list(
  HSC = c("HLF", "CD34"),
  MPP = c("KIT", "GATA2"),
  Cycling = c("MKI67", "TOP2A"),
  Myeloid = c("MPO", "LYZ", "AZU1")
)


# Select KRAS JMML populations
cells_of_interest <- c(
  "HSC",
  "GMP",
  "EMP",
  "CLP",
  "Early Eryth"
)

kras_sub <- subset(
  kras,
  subset = celltype_l2_predicted %in% cells_of_interest
)

# Check selected populations
table(kras_sub$celltype_l2_predicted)


# Plot marker expression
DotPlot(
  kras_sub,
  features = unique(unlist(markers_list)),
  group.by = "celltype_l2_predicted"
) + RotatedAxis()

# ===================================================
# HSC DE analysis (JMML HSC vs aneuploid/diploid HSC)




# Subset HSC populations
kras_hsc <- subset(
  kras,
  subset = celltype_l2_predicted == "HSC"
)

JMML_hsc <- subset(
  combined_JMML_AML_harmony,
  subset = sample == "2" &
    celltype_refined == "HSC" &
    copykat.pred.2 %in% c("diploid", "aneuploid")
)

table(JMML_hsc$copykat.pred.2)



# Match genes
common_genes_hsc <- intersect(
  rownames(kras_hsc),
  rownames(JMML_hsc)
)


# Extract counts
kras_counts <- LayerData(
  kras_hsc,
  assay = "originalexp",
  layer = "counts"
)[common_genes_hsc, ]

patient_2_counts <- LayerData(
  JMML_hsc,
  assay = "RNA",
  layer = "counts.2"
)[common_genes_hsc, ]


# Creating of merged object
merged_counts_hsc <- cbind(kras_counts, patient_2_counts)
hsc_clean <- CreateSeuratObject(counts = merged_counts_hsc)


# Metadata
kras_meta <- data.frame(
  comparison_group = rep("HSC_JMML", ncol(kras_counts)),
  source = rep("KRAS", ncol(kras_counts)),
  population = rep("HSC", ncol(kras_counts))
)
rownames(kras_meta) <- colnames(kras_counts)

patient_2_meta <- data.frame(
  comparison_group = JMML_hsc$copykat.pred.2,
  source = rep("patient_2", ncol(patient_2_counts)),
  population = rep("HSC", ncol(patient_2_counts))
)
rownames(patient_2_meta) <- colnames(patient_2_counts)

hsc_clean@meta.data <- rbind(kras_meta, patient_2_meta)


# 6. Normalize + identities
hsc_clean <- NormalizeData(hsc_clean)

DefaultAssay(hsc_clean) <- "RNA"
Idents(hsc_clean) <- "comparison_group"

table(Idents(hsc_clean))





# DE
deg_hsc_jmml_vs_aneu <- FindMarkers(
  hsc_clean,
  ident.1 = "HSC_JMML",
  ident.2 = "aneuploid",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

deg_hsc_jmml_vs_dip <- FindMarkers(
  hsc_clean,
  ident.1 = "HSC_JMML",
  ident.2 = "diploid",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

deg_hsc_aneu_vs_dip <- FindMarkers(
  hsc_clean,
  ident.1 = "aneuploid",
  ident.2 = "diploid",
  logfc.threshold = 0.25,
  min.pct = 0.1
)


# Filter significant genes
split_deg_sig <- function(deg_table) {
  deg_table$gene <- rownames(deg_table)
  
  up <- deg_table[
    deg_table$avg_log2FC > 0.25 & deg_table$p_val_adj < 0.05,
  ]
  up <- up[order(up$p_val_adj), ]
  
  down <- deg_table[
    deg_table$avg_log2FC < -0.25 & deg_table$p_val_adj < 0.05,
  ]
  down <- down[order(down$p_val_adj), ]
  
  list(up = up, down = down)
}

res_jmml_vs_aneu <- split_deg_sig(deg_hsc_jmml_vs_aneu)
res_jmml_vs_dip  <- split_deg_sig(deg_hsc_jmml_vs_dip)
res_aneu_vs_dip  <- split_deg_sig(deg_hsc_aneu_vs_dip)


# 9. Save results
wb <- createWorkbook()

addWorksheet(wb, "JMML_vs_aneu_up")
writeData(wb, "JMML_vs_aneu_up", res_jmml_vs_aneu$up)

addWorksheet(wb, "JMML_vs_aneu_down")
writeData(wb, "JMML_vs_aneu_down", res_jmml_vs_aneu$down)

addWorksheet(wb, "JMML_vs_dip_up")
writeData(wb, "JMML_vs_dip_up", res_jmml_vs_dip$up)

addWorksheet(wb, "JMML_vs_dip_down")
writeData(wb, "JMML_vs_dip_down", res_jmml_vs_dip$down)

addWorksheet(wb, "aneu_vs_dip_up")
writeData(wb, "aneu_vs_dip_up", res_aneu_vs_dip$up)

addWorksheet(wb, "aneu_vs_dip_down")
writeData(wb, "aneu_vs_dip_down", res_aneu_vs_dip$down)

saveWorkbook(
  wb,
  "DEG_HSC.xlsx",
  overwrite = TRUE
)


# Check the results
head(res_jmml_vs_aneu$up, 20)
head(res_jmml_vs_aneu$down, 20)

head(res_jmml_vs_dip$up, 20)
head(res_jmml_vs_dip$down, 20)

head(res_aneu_vs_dip$up, 20)
head(res_aneu_vs_dip$down, 20)


# ================================================
# GP DE analysis (JMML GP vs aneuploid/diploid GP)



# DE
deg_gp_jmml_vs_aneu <- FindMarkers(
  gp_clean,
  ident.1 = "GP_JMML",
  ident.2 = "aneuploid",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

deg_gp_jmml_vs_dip <- FindMarkers(
  gp_clean,
  ident.1 = "GP_JMML",
  ident.2 = "diploid",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

deg_gp_aneu_vs_dip <- FindMarkers(
  gp_clean,
  ident.1 = "aneuploid",
  ident.2 = "diploid",
  logfc.threshold = 0.25,
  min.pct = 0.1
)


# Filter significant genes
res_gp_jmml_vs_aneu <- split_deg_sig(deg_gp_jmml_vs_aneu)
res_gp_jmml_vs_dip  <- split_deg_sig(deg_gp_jmml_vs_dip)
res_gp_aneu_vs_dip  <- split_deg_sig(deg_gp_aneu_vs_dip)




# Save results
wb_gp <- createWorkbook()

addWorksheet(wb_gp, "JMML_vs_aneu_up")
writeData(wb_gp, "JMML_vs_aneu_up", res_gp_jmml_vs_aneu$up)

addWorksheet(wb_gp, "JMML_vs_aneu_down")
writeData(wb_gp, "JMML_vs_aneu_down", res_gp_jmml_vs_aneu$down)

addWorksheet(wb_gp, "JMML_vs_dip_up")
writeData(wb_gp, "JMML_vs_dip_up", res_gp_jmml_vs_dip$up)

addWorksheet(wb_gp, "JMML_vs_dip_down")
writeData(wb_gp, "JMML_vs_dip_down", res_gp_jmml_vs_dip$down)

addWorksheet(wb_gp, "aneu_vs_dip_up")
writeData(wb_gp, "aneu_vs_dip_up", res_gp_aneu_vs_dip$up)

addWorksheet(wb_gp, "aneu_vs_dip_down")
writeData(wb_gp, "aneu_vs_dip_down", res_gp_aneu_vs_dip$down)

saveWorkbook(
  wb_gp,
  "DEG_GP.xlsx",
  overwrite = TRUE
)


# Check the results
head(res_gp_jmml_vs_aneu$up, 20)
head(res_gp_jmml_vs_aneu$down, 20)

head(res_gp_jmml_vs_dip$up, 20)
head(res_gp_jmml_vs_dip$down, 20)

head(res_gp_aneu_vs_dip$up, 20)
head(res_gp_aneu_vs_dip$down, 20)





# =====================================================
# Monocyte DE analysis (JMML vs all tJMML monocytes)



patient_1_mono <- subset(
  patient_1_obj,
  subset = celltype_refined %in% c("Monocytes", "Monocyte-like", "IF Monocytes")
)

patient_3_mono <- subset(
  patient_3_obj,
  subset = celltype_refined %in% c("Monocytes", "Monocyte/DC-like", "IF Monocytes")
)

DefaultAssay(patient_1_mono) <- "RNA"
DefaultAssay(patient_3_mono) <- "RNA"

all_patient_mono <- merge(
  patient_1_mono,
  y = patient_3_mono,
  add.cell.ids = c("patient_1", "patient_3")
)

all_patient_mono$group <- "tJMML"


# Merge JMML and JMML-transformed monocytes
mono_de_all <- merge(
  jmml_mono,
  y = all_patient_mono,
  add.cell.ids = c("JMML", "tJMML")
)

mono_de_all <- NormalizeData(mono_de_all)
mono_de_all <- JoinLayers(mono_de_all, assay = "RNA")

DefaultAssay(mono_de_all) <- "RNA"
Idents(mono_de_all) <- "group"

table(Idents(mono_de_all))


# DE
deg_mono_jmml_vs_tJMML <- FindMarkers(
  mono_de_all,
  ident.1 = "JMML",
  ident.2 = "tJMML",
  logfc.threshold = 0.25,
  min.pct = 0.1
)


# Filter significant genes
res_mono_jmml_vs_tJMML <- split_deg_sig(deg_mono_jmml_vs_tJMML)


# saving
wb_mono <- createWorkbook()

addWorksheet(wb_mono, "JMML_vs_tJMML_up")
writeData(
  wb_mono,
  "JMML_vs_tJMML_up",
  res_mono_jmml_vs_tJMML$up
)

addWorksheet(wb_mono, "JMML_vs_tJMML_down")
writeData(
  wb_mono,
  "JMML_vs_tJMML_down",
  res_mono_jmml_vs_tJMML$down
)

saveWorkbook(
  wb_mono,
  "DEG_monocytes_JMML_vs_tJMML.xlsx",
  overwrite = TRUE
)


# check the results

head(res_mono_jmml_vs_tJMML$up, 20)
head(res_mono_jmml_vs_tJMML$down, 20)



# ===========================================
# Enrichment analysis for monocyte populations

library(readxl)
library(dplyr)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Input DEG file from monocyte DE analysis
input_file <- "DEG_monocytes_JMML_vs_tJMML.xlsx"

# Remove low-informative genes
clean_genes_mono <- function(genes) {
  genes <- unique(genes)
  genes[!(
    grepl("^MT-", genes)   |
      grepl("^RPL", genes)   |
      grepl("^RPS", genes)   |
      grepl("^AC", genes)    |
      grepl("^AL", genes)    |
      grepl("^LINC", genes)  |
      grepl("^MIR", genes)   |
      grepl("^SNOR", genes)  |
      grepl("^SCARNA", genes)
  )]
}

# Read significant DEGs
deg_up <- read_excel(input_file, sheet = "JMML_vs_tJMML_up")
deg_down <- read_excel(input_file, sheet = "JMML_vs_tJMML_down")

# Extract gene symbols
genes_up <- deg_up %>%
  dplyr::pull(gene) %>%
  na.omit() %>%
  as.character() %>%
  unique()

genes_down <- deg_down %>%
  dplyr::pull(gene) %>%
  na.omit() %>%
  as.character() %>%
  unique()

# Filter genes
genes_up_clean <- clean_genes_mono(genes_up)
genes_down_clean <- clean_genes_mono(genes_down)

# Convert gene symbols to ENTREZ IDs
up_df <- bitr(
  genes_up_clean,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_df <- bitr(
  genes_down_clean,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# GO enrichment
ego_up <- enrichGO(
  gene = up_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ego_down <- enrichGO(
  gene = down_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# KEGG enrichment
ekegg_up <- enrichKEGG(
  gene = up_df$ENTREZID,
  organism = "hsa"
)

ekegg_down <- enrichKEGG(
  gene = down_df$ENTREZID,
  organism = "hsa"
)

# Preview results
head(as.data.frame(ego_up), 10)
head(as.data.frame(ego_down), 10)

head(as.data.frame(ekegg_up), 10)
head(as.data.frame(ekegg_down), 10)

# Save enrichment results
write.xlsx(
  list(
    JMML_up_genes_clean = data.frame(gene = genes_up_clean),
    tJMML_up_genes_clean = data.frame(gene = genes_down_clean),
    GO_JMML_up = as.data.frame(ego_up),
    GO_tJMML_up = as.data.frame(ego_down),
    KEGG_JMML_up = as.data.frame(ekegg_up),
    KEGG_tJMML_up = as.data.frame(ekegg_down)
  ),
  file = "Monocytes_JMML_vs_tJMML_enrichment.xlsx"
)

# GO and KEGG dotplots
dotplot(
  ego_up,
  showCategory = 10,
  title = "GO BP: JMML monocytes up"
)

dotplot(
  ego_down,
  showCategory = 10,
  title = "GO BP: patient monocytes up"
)

dotplot(
  ekegg_up,
  showCategory = 10,
  title = "KEGG: JMML monocytes up"
)

dotplot(
  ekegg_down,
  showCategory = 10,
  title = "KEGG: patient monocytes up"
  
  
)


# ===========================================================
# Enrichment analysis for Granulocytic Progenitors population



# Input DEG file from GP DE analysis
input_file <- "DEG_GP.xlsx"

# Remove low-informative genes
clean_genes <- function(genes) {
  genes <- unique(genes)
  genes[!(
    grepl("^MT-", genes)   |
      grepl("^RPL", genes)   |
      grepl("^RPS", genes)   |
      grepl("^AC", genes)    |
      grepl("^AL", genes)    |
      grepl("^LINC", genes)  |
      grepl("^MIR", genes)   |
      grepl("^SNOR", genes)  |
      grepl("^SCARNA", genes)
  )]
}

# Run enrichment for one pair of DEG sheets
run_enrichment_pair <- function(sheet_up, sheet_down) {
  deg_up <- read_excel(input_file, sheet = sheet_up)
  deg_down <- read_excel(input_file, sheet = sheet_down)
  
  genes_up <- deg_up %>%
    dplyr::pull(gene) %>%
    na.omit() %>%
    as.character() %>%
    unique()
  
  genes_down <- deg_down %>%
    dplyr::pull(gene) %>%
    na.omit() %>%
    as.character() %>%
    unique()
  
  genes_up_clean <- clean_genes(genes_up)
  genes_down_clean <- clean_genes(genes_down)
  
  up_df <- bitr(
    genes_up_clean,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  down_df <- bitr(
    genes_down_clean,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  ego_up <- enrichGO(
    gene = up_df$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ego_down <- enrichGO(
    gene = down_df$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ekegg_up <- enrichKEGG(
    gene = up_df$ENTREZID,
    organism = "hsa"
  )
  
  ekegg_down <- enrichKEGG(
    gene = down_df$ENTREZID,
    organism = "hsa"
  )
  
  list(
    up_genes_clean = data.frame(gene = genes_up_clean),
    down_genes_clean = data.frame(gene = genes_down_clean),
    GO_up_obj = ego_up,
    GO_down_obj = ego_down,
    KEGG_up_obj = ekegg_up,
    KEGG_down_obj = ekegg_down,
    GO_up = as.data.frame(ego_up),
    GO_down = as.data.frame(ego_down),
    KEGG_up = as.data.frame(ekegg_up),
    KEGG_down = as.data.frame(ekegg_down)
  )
}

# Run enrichment for all GP comparisons
res_jmml_vs_aneu <- run_enrichment_pair(
  "JMML_vs_aneu_up",
  "JMML_vs_aneu_down"
)

res_jmml_vs_dip <- run_enrichment_pair(
  "JMML_vs_dip_up",
  "JMML_vs_dip_down"
)

res_aneu_vs_dip <- run_enrichment_pair(
  "aneu_vs_dip_up",
  "aneu_vs_dip_down"
)

# Save enrichment results
write.xlsx(
  list(
    JMML_vs_aneu_up_genes_clean = res_jmml_vs_aneu$up_genes_clean,
    JMML_vs_aneu_down_genes_clean = res_jmml_vs_aneu$down_genes_clean,
    JMML_vs_aneu_GO_up = res_jmml_vs_aneu$GO_up,
    JMML_vs_aneu_GO_down = res_jmml_vs_aneu$GO_down,
    JMML_vs_aneu_KEGG_up = res_jmml_vs_aneu$KEGG_up,
    JMML_vs_aneu_KEGG_down = res_jmml_vs_aneu$KEGG_down,
    
    JMML_vs_dip_up_genes_clean = res_jmml_vs_dip$up_genes_clean,
    JMML_vs_dip_down_genes_clean = res_jmml_vs_dip$down_genes_clean,
    JMML_vs_dip_GO_up = res_jmml_vs_dip$GO_up,
    JMML_vs_dip_GO_down = res_jmml_vs_dip$GO_down,
    JMML_vs_dip_KEGG_up = res_jmml_vs_dip$KEGG_up,
    JMML_vs_dip_KEGG_down = res_jmml_vs_dip$KEGG_down,
    
    aneu_vs_dip_up_genes_clean = res_aneu_vs_dip$up_genes_clean,
    aneu_vs_dip_down_genes_clean = res_aneu_vs_dip$down_genes_clean,
    aneu_vs_dip_GO_up = res_aneu_vs_dip$GO_up,
    aneu_vs_dip_GO_down = res_aneu_vs_dip$GO_down,
    aneu_vs_dip_KEGG_up = res_aneu_vs_dip$KEGG_up,
    aneu_vs_dip_KEGG_down = res_aneu_vs_dip$KEGG_down
  ),
  file = "GP_enrichment.xlsx"
)

# Preview results
head(res_jmml_vs_aneu$GO_up, 10)
head(res_jmml_vs_aneu$GO_down, 10)
head(res_jmml_vs_aneu$KEGG_up, 10)
head(res_jmml_vs_aneu$KEGG_down, 10)

head(res_jmml_vs_dip$GO_up, 10)
head(res_jmml_vs_dip$GO_down, 10)
head(res_jmml_vs_dip$KEGG_up, 10)
head(res_jmml_vs_dip$KEGG_down, 10)

head(res_aneu_vs_dip$GO_up, 10)
head(res_aneu_vs_dip$GO_down, 10)
head(res_aneu_vs_dip$KEGG_up, 10)
head(res_aneu_vs_dip$KEGG_down, 10)

# GO dotplots
dotplot(
  res_jmml_vs_aneu$GO_up_obj,
  showCategory = 10,
  title = "GO BP: GP JMML up vs aneuploid"
)

dotplot(
  res_jmml_vs_aneu$GO_down_obj,
  showCategory = 10,
  title = "GO BP: GP aneuploid up vs JMML"
)

dotplot(
  res_jmml_vs_dip$GO_up_obj,
  showCategory = 10,
  title = "GO BP: GP JMML up vs diploid"
)

dotplot(
  res_jmml_vs_dip$GO_down_obj,
  showCategory = 10,
  title = "GO BP: GP diploid up vs JMML"
)

dotplot(
  res_aneu_vs_dip$GO_up_obj,
  showCategory = 10,
  title = "GO BP: GP aneuploid up vs diploid"
)

dotplot(
  res_aneu_vs_dip$GO_down_obj,
  showCategory = 10,
  title = "GO BP: GP diploid up vs aneuploid"
)

# KEGG dotplots
dotplot(
  res_jmml_vs_aneu$KEGG_up_obj,
  showCategory = 10,
  title = "KEGG: GP JMML up vs aneuploid"
)

dotplot(
  res_jmml_vs_aneu$KEGG_down_obj,
  showCategory = 10,
  title = "KEGG: GP aneuploid up vs JMML"
)

dotplot(
  res_jmml_vs_dip$KEGG_up_obj,
  showCategory = 10,
  title = "KEGG: GP JMML up vs diploid"
)

dotplot(
  res_jmml_vs_dip$KEGG_down_obj,
  showCategory = 10,
  title = "KEGG: GP diploid up vs JMML"
)

dotplot(
  res_aneu_vs_dip$KEGG_up_obj,
  showCategory = 10,
  title = "KEGG: GP aneuploid up vs diploid"
)

dotplot(
  res_aneu_vs_dip$KEGG_down_obj,
  showCategory = 10,
  title = "KEGG: GP diploid up vs aneuploid"
)


# ===========================================================
# Enrichment analysis for Hematopoietic stem cells population



# Input DEG file from HSC DE analysis
input_file <- "DEG_HSC.xlsx"

# Remove low-informative genes
clean_genes <- function(genes) {
  genes <- unique(genes)
  genes[!(
    grepl("^MT-", genes)   |
      grepl("^RPL", genes)   |
      grepl("^RPS", genes)   |
      grepl("^AC", genes)    |
      grepl("^AL", genes)    |
      grepl("^LINC", genes)  |
      grepl("^MIR", genes)   |
      grepl("^SNOR", genes)  |
      grepl("^SCARNA", genes)
  )]
}

# Run enrichment for one pair of DEG sheets
run_enrichment_pair <- function(sheet_up, sheet_down) {
  deg_up <- read_excel(input_file, sheet = sheet_up)
  deg_down <- read_excel(input_file, sheet = sheet_down)
  
  genes_up <- deg_up %>%
    dplyr::pull(gene) %>%
    na.omit() %>%
    as.character() %>%
    unique()
  
  genes_down <- deg_down %>%
    dplyr::pull(gene) %>%
    na.omit() %>%
    as.character() %>%
    unique()
  
  genes_up_clean <- clean_genes(genes_up)
  genes_down_clean <- clean_genes(genes_down)
  
  up_df <- bitr(
    genes_up_clean,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  down_df <- bitr(
    genes_down_clean,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  ego_up <- NULL
  if (nrow(up_df) > 1) {
    ego_up <- enrichGO(
      gene = up_df$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      readable = TRUE
    )
  }
  
  ego_down <- NULL
  if (nrow(down_df) > 1) {
    ego_down <- enrichGO(
      gene = down_df$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      readable = TRUE
    )
  }
  
  ekegg_up <- NULL
  if (nrow(up_df) > 1) {
    ekegg_up <- tryCatch(
      enrichKEGG(
        gene = up_df$ENTREZID,
        organism = "hsa"
      ),
      error = function(e) NULL
    )
  }
  
  ekegg_down <- NULL
  if (nrow(down_df) > 1) {
    ekegg_down <- tryCatch(
      enrichKEGG(
        gene = down_df$ENTREZID,
        organism = "hsa"
      ),
      error = function(e) NULL
    )
  }
  
  list(
    up_genes_raw = data.frame(gene = genes_up),
    down_genes_raw = data.frame(gene = genes_down),
    up_genes_clean = data.frame(gene = genes_up_clean),
    down_genes_clean = data.frame(gene = genes_down_clean),
    up_mapped = up_df,
    down_mapped = down_df,
    GO_up_obj = ego_up,
    GO_down_obj = ego_down,
    KEGG_up_obj = ekegg_up,
    KEGG_down_obj = ekegg_down,
    GO_up = if (!is.null(ego_up)) as.data.frame(ego_up) else data.frame(),
    GO_down = if (!is.null(ego_down)) as.data.frame(ego_down) else data.frame(),
    KEGG_up = if (!is.null(ekegg_up)) as.data.frame(ekegg_up) else data.frame(),
    KEGG_down = if (!is.null(ekegg_down)) as.data.frame(ekegg_down) else data.frame()
  )
}

# Run enrichment for all HSC comparisons
res_jmml_vs_aneu <- run_enrichment_pair(
  "JMML_vs_aneu_up",
  "JMML_vs_aneu_down"
)

res_jmml_vs_dip <- run_enrichment_pair(
  "JMML_vs_dip_up",
  "JMML_vs_dip_down"
)

res_aneu_vs_dip <- run_enrichment_pair(
  "aneu_vs_dip_up",
  "aneu_vs_dip_down"
)

# Save enrichment results
write.xlsx(
  list(
    JMML_vs_aneu_up_genes_raw = res_jmml_vs_aneu$up_genes_raw,
    JMML_vs_aneu_down_genes_raw = res_jmml_vs_aneu$down_genes_raw,
    JMML_vs_aneu_up_genes_clean = res_jmml_vs_aneu$up_genes_clean,
    JMML_vs_aneu_down_genes_clean = res_jmml_vs_aneu$down_genes_clean,
    JMML_vs_aneu_up_mapped = res_jmml_vs_aneu$up_mapped,
    JMML_vs_aneu_down_mapped = res_jmml_vs_aneu$down_mapped,
    JMML_vs_aneu_GO_up = res_jmml_vs_aneu$GO_up,
    JMML_vs_aneu_GO_down = res_jmml_vs_aneu$GO_down,
    JMML_vs_aneu_KEGG_up = res_jmml_vs_aneu$KEGG_up,
    JMML_vs_aneu_KEGG_down = res_jmml_vs_aneu$KEGG_down,
    
    JMML_vs_dip_up_genes_raw = res_jmml_vs_dip$up_genes_raw,
    JMML_vs_dip_down_genes_raw = res_jmml_vs_dip$down_genes_raw,
    JMML_vs_dip_up_genes_clean = res_jmml_vs_dip$up_genes_clean,
    JMML_vs_dip_down_genes_clean = res_jmml_vs_dip$down_genes_clean,
    JMML_vs_dip_up_mapped = res_jmml_vs_dip$up_mapped,
    JMML_vs_dip_down_mapped = res_jmml_vs_dip$down_mapped,
    JMML_vs_dip_GO_up = res_jmml_vs_dip$GO_up,
    JMML_vs_dip_GO_down = res_jmml_vs_dip$GO_down,
    JMML_vs_dip_KEGG_up = res_jmml_vs_dip$KEGG_up,
    JMML_vs_dip_KEGG_down = res_jmml_vs_dip$KEGG_down,
    
    aneu_vs_dip_up_genes_raw = res_aneu_vs_dip$up_genes_raw,
    aneu_vs_dip_down_genes_raw = res_aneu_vs_dip$down_genes_raw,
    aneu_vs_dip_up_genes_clean = res_aneu_vs_dip$up_genes_clean,
    aneu_vs_dip_down_genes_clean = res_aneu_vs_dip$down_genes_clean,
    aneu_vs_dip_up_mapped = res_aneu_vs_dip$up_mapped,
    aneu_vs_dip_down_mapped = res_aneu_vs_dip$down_mapped,
    aneu_vs_dip_GO_up = res_aneu_vs_dip$GO_up,
    aneu_vs_dip_GO_down = res_aneu_vs_dip$GO_down,
    aneu_vs_dip_KEGG_up = res_aneu_vs_dip$KEGG_up,
    aneu_vs_dip_KEGG_down = res_aneu_vs_dip$KEGG_down
  ),
  file = "HSC_enrichment.xlsx"
)

# Preview results
head(res_jmml_vs_aneu$GO_up, 10)
head(res_jmml_vs_aneu$GO_down, 10)
head(res_jmml_vs_aneu$KEGG_up, 10)
head(res_jmml_vs_aneu$KEGG_down, 10)

head(res_jmml_vs_dip$GO_up, 10)
head(res_jmml_vs_dip$GO_down, 10)
head(res_jmml_vs_dip$KEGG_up, 10)
head(res_jmml_vs_dip$KEGG_down, 10)

head(res_aneu_vs_dip$GO_up, 10)
head(res_aneu_vs_dip$GO_down, 10)
head(res_aneu_vs_dip$KEGG_up, 10)
head(res_aneu_vs_dip$KEGG_down, 10)

# GO dotplots
if (!is.null(res_jmml_vs_aneu$GO_up_obj)) {
  print(dotplot(
    res_jmml_vs_aneu$GO_up_obj,
    showCategory = 10,
    title = "GO BP: HSC JMML up vs aneuploid"
  ))
}

if (!is.null(res_jmml_vs_aneu$GO_down_obj)) {
  print(dotplot(
    res_jmml_vs_aneu$GO_down_obj,
    showCategory = 10,
    title = "GO BP: HSC aneuploid up vs JMML"
  ))
}

if (!is.null(res_jmml_vs_dip$GO_up_obj)) {
  print(dotplot(
    res_jmml_vs_dip$GO_up_obj,
    showCategory = 10,
    title = "GO BP: HSC JMML up vs diploid"
  ))
}

if (!is.null(res_jmml_vs_dip$GO_down_obj)) {
  print(dotplot(
    res_jmml_vs_dip$GO_down_obj,
    showCategory = 10,
    title = "GO BP: HSC diploid up vs JMML"
  ))
}

if (!is.null(res_aneu_vs_dip$GO_up_obj)) {
  print(dotplot(
    res_aneu_vs_dip$GO_up_obj,
    showCategory = 10,
    title = "GO BP: HSC aneuploid up vs diploid"
  ))
}

if (!is.null(res_aneu_vs_dip$GO_down_obj)) {
  print(dotplot(
    res_aneu_vs_dip$GO_down_obj,
    showCategory = 10,
    title = "GO BP: HSC diploid up vs aneuploid"
  ))
}

# KEGG dotplots
if (!is.null(res_jmml_vs_aneu$KEGG_up_obj)) {
  print(dotplot(
    res_jmml_vs_aneu$KEGG_up_obj,
    showCategory = 10,
    title = "KEGG: HSC JMML up vs aneuploid"
  ))
}

if (!is.null(res_jmml_vs_aneu$KEGG_down_obj)) {
  print(dotplot(
    res_jmml_vs_aneu$KEGG_down_obj,
    showCategory = 10,
    title = "KEGG: HSC aneuploid up vs JMML"
  ))
}

if (!is.null(res_jmml_vs_dip$KEGG_up_obj)) {
  print(dotplot(
    res_jmml_vs_dip$KEGG_up_obj,
    showCategory = 10,
    title = "KEGG: HSC JMML up vs diploid"
  ))
}

if (!is.null(res_jmml_vs_dip$KEGG_down_obj)) {
  print(dotplot(
    res_jmml_vs_dip$KEGG_down_obj,
    showCategory = 10,
    title = "KEGG: HSC diploid up vs JMML"
  ))
}

if (!is.null(res_aneu_vs_dip$KEGG_up_obj)) {
  print(dotplot(
    res_aneu_vs_dip$KEGG_up_obj,
    showCategory = 10,
    title = "KEGG: HSC aneuploid up vs diploid"
  ))
}

if (!is.null(res_aneu_vs_dip$KEGG_down_obj)) {
  print(dotplot(
    res_aneu_vs_dip$KEGG_down_obj,
    showCategory = 10,
    title = "KEGG: HSC diploid up vs aneuploid"
  ))
}
