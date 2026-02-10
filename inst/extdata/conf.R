set.seed(1234)
options(future.globals.maxSize = 2000 * 1024^2)
options(ggrepel.max.overlaps = 1000)
max_dim <- 100
n_dim <- 50
dims <- seq(30)
dims <- dims[-c(1)]
species <- "mouse"
assay <- "ATAC"
if (assay == "ATAC") {
    assay2 <- "ATAC_V5"
} else {
    assay2 <- "RNA"
}
SCT <- FALSE
reduction <- "lsi"

if (species == "human") {
    func_format <- identity
} else {
    func_format <- str_to_sentence
}

target_type <- "KO"
pal_discrete_sc <- palette_discrete()[-7]

best_resolution <- 0.35
snn_cluster <- paste0(
    assay,
    "_snn_res.",
    best_resolution
)

