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
        "Cd19", "Cd79a", "Cd79b", "Ighm", "Igkc", "Il4i1", "Ms4a1",
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
