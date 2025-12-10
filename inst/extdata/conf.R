set.seed(1234)
options(future.globals.maxSize = 2000 * 1024^2)
options(ggrepel.max.overlaps = 1000)
max_dim <- 100
n_dim <- 30
species <- "mouse"
assay <- "RNA"
SCT <- FALSE

if (species == "human") {
    func_format <- identity
} else {
    func_format <- str_to_sentence
}
