set.seed(1234)
options(future.globals.maxSize = 2000 * 1024^2)
options(ggrepel.max.overlaps = 1000)
max_dim <- 100
n_dim <- 30
dims <- seq(n_dim)
species <- "mouse"
assay <- "RNA"
SCT <- FALSE
if (assay == "ATAC") {
    assay2 <- "RNA"
    SCT <- FALSE
    reduction <- "lsi"
    path <- "ATAC"
    dims <- dims[-c(1)]
} else {
    if (SCT) {
        assay <- "SCT"
    }
    assay2 <- assay
    reduction <- "pca"
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

best_resolution <- 0.35
snn_cluster <- paste0(
    assay,
    "_snn_res.",
    best_resolution
)

path_root <- file.path("C:", "Users", "etien")
path_root_data <- file.path(path_root, "DATA")
path_data <- file.path(path_root_data, "dobino", path)
path_project <- file.path(path_data, "cellranger_pasteur")
path_fig <- file.path(path_root, "bin", "GimmmeMySC", "vignettes",  knitr::opts_chunk$get("fig.path"))
seurat_dataset <- "rna_allcell4"
# seurat_dataset <- "rna_bcell_integrated2abc"
# seurat_dataset <- paste0(seurat_dataset, "_", target_type)
