

# CopyKAT analysis for all patients


library(pheatmap)

# ---------------------------
# Patient 2
# ---------------------------

# Load CopyKAT
copykat_patient_2 <- read.table(
  "patient_2_CNV/Patient_2_FULL_copykat_prediction.txt",
  header = TRUE
)

# Add predictions
combined_JMML_AML_harmony$copykat.pred.patient_2 <- NA

m_patient_2 <- match(
  rownames(combined_JMML_AML_harmony@meta.data),
  copykat_patient_2$cell.names
)

combined_JMML_AML_harmony$copykat.pred.patient_2 <- copykat_patient_2$copykat.pred[m_patient_2]

# Subset
patient_2_obj <- subset(
  combined_JMML_AML_harmony,
  subset = sample == "patient_2"
)

# Proportions
prop_table_patient_2 <- round(
  prop.table(
    table(
      patient_2_obj$copykat.pred.patient_2,
      patient_2_obj$celltype_refined,
      useNA = "no"
    ),
    margin = 2
  ),
  3
)

# Selected populations
cells_of_interest_patient_2 <- c(
  "HSC", "MPP", "Granulocytic Progenitors"
)

prop_small_patient_2 <- prop_table_patient_2[, cells_of_interest_patient_2, drop = FALSE]

# Heatmap
pheatmap(
  prop_small_patient_2,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "CopyKAT proportions in Patient 2"
)

# UMAP
DimPlot(
  patient_2_obj,
  group.by = "copykat.pred.patient_2",
  reduction = "umap"
)


# ---------------------------
# Patient 1
# ---------------------------

# Load CopyKAT
copykat_patient_1 <- read.table(
  "copykat_patient_3_1/Patient_1_FULL_copykat_prediction.txt",
  header = TRUE
)

# Add predictions
combined_JMML_AML_harmony$copykat.pred.patient_1 <- NA

m_patient_1 <- match(
  rownames(combined_JMML_AML_harmony@meta.data),
  copykat_patient_1$cell.names
)

combined_JMML_AML_harmony$copykat.pred.patient_1 <- copykat_patient_1$copykat.pred[m_patient_1]

# Subset
patient_1_obj <- subset(
  combined_JMML_AML_harmony,
  subset = sample == "patient_1"
)

# Proportions
prop_table_patient_1 <- round(
  prop.table(
    table(
      patient_1_obj$copykat.pred.patient_1,
      patient_1_obj$celltype_refined,
      useNA = "no"
    ),
    margin = 2
  ),
  3
)

# Selected populations
cells_of_interest_patient_1 <- c(
  "Activated macrophages",
  "Macrophages",
  "cDC",
  "IF Monocytes",
  "Monocyte-like",
  "Monocytes",
  "Granulocytic Progenitors",
  "Myeloid Blasts",
  "Erythroid Progenitors",
  "MPP"
)

cells_present_patient_1 <- intersect(
  cells_of_interest_patient_1,
  colnames(prop_table_patient_1)
)

prop_small_patient_1 <- prop_table_patient_1[, cells_present_patient_1, drop = FALSE]

# Heatmap
pheatmap(
  prop_small_patient_1,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "CopyKAT proportions in Patient 1"
)

# UMAP
DimPlot(
  patient_1_obj,
  group.by = "copykat.pred.patient_1",
  reduction = "umap"
)


# ---------------------------
# Patient 3
# ---------------------------

# Load CopyKAT
copykat_patient_3 <- read.table(
  "copykat_patient_3_1/Patient_3_FULL_copykat_prediction.txt",
  header = TRUE
)

# Add predictions
combined_JMML_AML_harmony$copykat.pred.patient_3 <- NA

m_patient_3 <- match(
  rownames(combined_JMML_AML_harmony@meta.data),
  copykat_patient_3$cell.names
)

combined_JMML_AML_harmony$copykat.pred.patient_3 <- copykat_patient_3$copykat.pred[m_patient_3]

# Subset
patient_3_obj <- subset(
  combined_JMML_AML_harmony,
  subset = sample == "patient_3"
)

# Proportions
prop_table_patient_3 <- round(
  prop.table(
    table(
      patient_3_obj$copykat.pred.patient_3,
      patient_3_obj$celltype_refined,
      useNA = "no"
    ),
    margin = 2
  ),
  3
)

# Selected populations
cells_of_interest_patient_3 <- c(
  "cDC", "IF Monocytes",
  "Macrophages", "Monocyte/DC-like", "Monocytes",
  "Cycling Myeloid", "MPP", "HSC", "pDC"
)

cells_present_patient_3 <- intersect(
  cells_of_interest_patient_3,
  colnames(prop_table_patient_3)
)

prop_small_patient_3 <- prop_table_patient_3[, cells_present_patient_3, drop = FALSE]

# Heatmap
pheatmap(
  prop_small_patient_3,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "CopyKAT proportions in Patient 3"
)

# UMAP
DimPlot(
  patient_3_obj,
  group.by = "copykat.pred.patient_3",
  reduction = "umap"
)
