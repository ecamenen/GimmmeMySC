#### B cell markers
l_bcell_markers <- list(
    General = c(
        "Pax5", "Ebf1", "Spi1", "Ikzf1", "Ikzf3", # facteurs de transcription clés
        "Pou2f2", "Ets1", "Bach2",                # maintien du programme B
        "Ptprc", "Fnip1", "Hdac5", "Hdac9"
    ),

    FrA = c( # B220+CD43+ très précoces
        "Enpep", "Foxp1", "Il7r", "Kit", "Spn"
    ),

    FrB = c( # Pro-B précoces
        "Cdkn3", "Dntt", "Igll1", "Irf8", "Lef1", "Rag1", "Rag2", "Vpreb1a", "Vpreb1b", "Il7r"
    ),

    FrC = c( # Pro-B tardives (pré-BCR en recombinaison)
        "Polm", "Rag1", "Rag2", "Spib"
    ),

    FrD = c( # Pré-B prolifératives
        "Ccna2", "Ccnb1", "Ccnd2", "Ccne1", "Fancg", "Gfi1b",
        "Il2ra", "Mki67", "Pcna", "Polr2a", "Reln", "Top2a"
    ),

    FrE = c( # B immatures IgM+
        "Cd19", "Cd79a", "Cd79b", "Ighm", "Igkc", "Iglc2", "Iglc3", "Bcl6"
    ),

    FrF = c( # B matures (transit moelle → périphérie)
        "Il4i1", "Ms4a1",
        "Nfkbiz", "Notch2", "Ptk2b", "Ptprj", "Tlr9"
    ),

    Plasmablast = c(
        "Dnajb9", "Mmp14", "Prdm1", "Sdc1", "Xbp1"
    ),

    MemoryB = c(
        "Ahr", "Tnfaip3", "Zfp36l1", "Cd27", "Cd80", "Cd86", "Ms4a1"
    )

    # Regulatory<- c(
    #   "Clcf1", "Irf2bp2"
    # )
)

usethis::use_data(l_bcell_markers, overwrite = TRUE)


##### Microarray markers

fileIn <- file.path(path, "B_cell_subet_Marker_Genes-JCB.xlsx")
name <- "Gene name"

jcb_markers <- list.map(
    excel_sheets(fileIn),
    f(i) ~ read_excel(fileIn, sheet = i) %>%
        pull(!!sym(name)) %>%
        na.omit() %>%
        func_format() %>%
        .[. %in% Features(seurat)]
)
usethis::use_data(jcb_markers, overwrite = TRUE)

#### Data integration

path <- file.path(
    golem::get_golem_wd(),
    "inst",
    "extdata"
)

wt <- readRDS(file.path(path, "all_cell_wt.rds"))
ko <- readRDS(file.path(path, "all_cell_ko.rds"))

wt@assays$SCT@scale.data <- matrix()
ko@assays$SCT@scale.data <- matrix()

seurat <- merge(wt, ko)
# seurat[["RNA"]] <- JoinLayers(seurat[["RNA"]])
# usethis::use_data(bcell_integrated, overwrite = TRUE)
seurat <- subset(seurat, cell_type_formatted == "B cell")
save(seurat, file = file.path(path, "integrated_bcells2.rda"))

################################################

seurat <- merge(wt, ko)
seurat$cell_type_formatted <- remove_parenthesis(seurat$cell_subtype)
seurat$cell_subtype_formatted <- keep_parenthesis(seurat$cell_subtype) %>% factor(levels = reorder_celltype(.))
seurat$cell_subtype_formatted2 <- seurat$cell_subtype_formatted
seurat[["percent_ribosomal"]] <- PercentageFeatureSet(seurat, pattern = "^((RP[LS])|(Rp[ls]))")
seurat$Type <- factor(seurat$Type, levels = c("WT", "KO"))

seurat <- subset(seurat, cell_type_formatted %in% "B cell")
seurat <- subset(seurat, !str_detect(cell_subtype_formatted, "FRA|CLP"))

save(seurat, file = file.path(path, "integrated_bcells5.rda"))

################################################
library("readxl")

excel_file <- "Signatures-Articles.xlsx"
excel_file <- "B_cell_subet_Marker_Genes-JCB.xlsx"
file_path <- file.path(golem::get_golem_wd(), "inst", "extdata", excel_file)

term2gene_microarray <- excel_sheets(file_path) %>%
    list.map(~read_excel(file_path, sheet = ., col_names = TRUE)[[2]] %>% str_to_sentence()) %>% set_names(names(.) %>% str_remove_all(" vs_all")) %>%
enframe(name = "term", value = "name") %>%
    unnest(cols = c(name))

usethis::use_data(term2gene_microarray, overwrite = TRUE)
