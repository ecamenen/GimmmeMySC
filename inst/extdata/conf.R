set.seed(1234)
options(future.globals.maxSize = 2000 * 1024^2)
options(ggrepel.max.overlaps = 1000)
max_dim <- 100
n_dim <- 50
dims <- seq(30)
dims <- dims[-c(1)]
species <- "mouse"
assay <- "SCT"
SCT <- TRUE
reduction <- "pca"
if (assay == "ATAC") {
    assay2 <- "RNA"
    SCT <- FALSE
    reduction <- "lsi"
    path <- "ATAC"
} else {
    assay2 <- assay
    path <- "RNA"
}

if (species == "human") {
    func_format <- identity
} else {
    func_format <- str_to_sentence
}

target_type <- "KO"
split.by <- "Type"
pal_discrete_sc <- palette_discrete()[-7] -> pal_discrete_sc2
# pal_discrete_sc <- palette_continuous()(15)
# pal_discrete_sc2 <-  palette_continuous()(29)

best_resolution <- 0.35
snn_cluster <- paste0(
    assay,
    "_snn_res.",
    best_resolution
)

path_root <- file.path("C:", "Users", "etien")
path_root_data <- file.path(path_root, "DATA")
path_data <- file.path(path_root_data, "dobino", path)
path_fig <- file.path(path_root, "bin", "GimmmeMySC", "vignettes",  knitr::opts_chunk$get("fig.path"))
# seurat_dataset <- paste0("rna_allcell", "_", target_type)
seurat_dataset <- "rna_bcell_integrated"

