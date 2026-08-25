set.seed(1234)
options(future.globals.maxSize = 2000 * 1024^2)
options(ggrepel.max.overlaps = 1000)
multiome <- TRUE
if (multiome) {
    reduction2 <- "wnn.umap"
} else {
    reduction2 <- "umap"
}
max_dim <- 100
n_dim <- 30
dims <- seq(n_dim)
species <- "mouse"
assay <- "RNA"
assay2 <- assay
SCT <- TRUE
integrated <- FALSE
if (assay == "ATAC") {
    SCT <- FALSE
    reduction <- "lsi"
    path <- "ATAC"
} else {
    if (SCT) {
        assay <- "SCT"
    }
    reduction <- "pca"
    path <- "RNA"
}

if (species == "human") {
    func_format <- identity
} else {
    func_format <- format_gene
}

target_type <- "KO"
split.by <- "Type"
pal_discrete_sc <- palette_discrete()[-7] -> pal_discrete_sc2

best_resolution <- 0.35
snn_cluster <- paste0(
    ifelse(multiome, "w", paste0(assay, "_")),
    "snn_res.",
    best_resolution
)

path_root <- file.path("C:", "Users", "etien")
path_root_data <- file.path(path_root, "DATA")
path_data0 <- file.path(path_root_data, "dobino")
path_data <- file.path(path_data0, path)
path_project <- file.path(path_data, "cellranger")
path_fig <- file.path(path_root, "bin", "GimmmeMySC", "vignettes",  knitr::opts_chunk$get("fig.path"))
seurat_dataset <- "rna_allcell4"
# seurat_dataset <- "rna_bcell_integrated2abc"
# seurat_dataset <- paste0(seurat_dataset, "_", target_type)
reorder_ident <- c("2", "5", "0", "1", "3", "4", "6", "7")

# l_samples <- list.dirs(path_project, recursive = FALSE) %>% basename() %>% .[!str_detect(., "__")]
l_samples <- c("WT1", "WT11")
l_samples <- str_remove_all(l_samples, "cellranger_onlyATAC_count_")
types  <- str_remove_all(l_samples, "_?\\d+")

col_inds <- list(
    # c("#FB9A99", "#E31A1C" ),
    c("#B2DF8A", "#33A02C")
)
color_types <- sapply(col_inds, last)

library("rlist")
col_inds <- list.mapv(unique(types), f(i, j) ~colorRampPalette(col_inds[[j]])(length(types[types == i]))) %>% as.character()

cell_id <- paste0(
    types,
    "_",
    str_remove_all(l_samples, types) %>%
        str_remove_all("_")
)


if (multiome) {
    techno <- "atac_"
} else {
    techno <- ""
}

features_supp0 <- c("unmapped", ifelse(techno == "atac_", "dup", "duplicate"), paste0("mitochondrial", ifelse(multiome, "_atac", "")))
features_supp <- c("percent_low_map_qc", paste0("percent_", features_supp0))

features_atac <- c(
    paste0("nCount", "_", assay),
    paste0("nFeature", "_", assay),
    "nucleosome_signal",
    "TSS.enrichment",
    "percent_reads_in_peaks",
    "blacklist_ratio"
)
