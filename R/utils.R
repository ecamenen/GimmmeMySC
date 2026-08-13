#' @export
plot_dim_cycle <- function(x, ncol = 2, ...) {
    Idents(x) <- "Phase"
    p1 <- plot_feature(
        x,
        features = "CC.Difference",
        pt.size = .5,
        ...
    ) %>% list.map(
        f(x) ~x + labs(title = NULL, color = "S to G2/M")
    )
    print(p1)

    plot_dim_cycle0(x, ncol = ncol, ...)
}

#' @export
plot_dim_cycle0 <- function(x, ncol = 2, ...) {
    Idents(x) <- "Phase"
    plot_dim(
        object = x,
        label = FALSE,
        title = TRUE,
        combine = FALSE,
        ncol = ncol,
        ...
    ) %>% list.map(
        f(x) ~x +
            labs(x = NULL, y = NULL) +
            scale_color_manual(
                values = palette_discrete(),
                limits = c("G1", "G2M", "S")
            )
    )
}

#' @export
plot_pca_qc <- function(x) {
    p1 <- plot_eig(x) +
        geom_hline(yintercept = mean(x), col = "gray50", lwd = 1, lty = 2) +
        ylab("% Explained\n variance")

    p2 <- cumsum(x) %>%
        plot_eig()+
        ylim(0, NA) +
        ylab("% Cumulated\n explained variance") +
        geom_hline(yintercept = 80, color = "gray50", lwd = 1, lty = 2)

    p1 + p2
}

#' @export
kable_cluster <-  function(
        x,
        clusters,
        file = NULL,
        assay = "RNA"
    ) {
    cl_size <- map_dfr(
        clusters,
        ~ cluster_size(x@meta.data[, .x]) %>%
            pull(1) %>%
            tibble(
                Clusters = .x,
                cluster = seq_along(.),
                inertia = .
            )
    ) %>% pivot_wider(
        names_from = cluster,
        values_from = inertia,
        names_prefix = ""
    ) %>%
        as.data.frame() %>%
        mutate(Clusters = str_remove_all(Clusters, paste0("(", assay, "_)?w?snn_res.")))

    cl_size[is.na(cl_size)] <- ""
    p <- kable0(cl_size)

    if (!is.null(file))
        save_kable(
            p,
            file = file,
            zoom = 2
        )
    print(p)
}

#' @export
kable_sc <- function(x, func = str_clean, digits = 3, file = NULL) {
    p <- x %>%
        set_colnames(colnames(.) %>% func()) %>%
        kable0(digits = digits)
    if (!is.null(file))
        save_kable(
            p,
            file = file,
            zoom = 2
        )
    p
}

#' @export
kable_ncell <- function(x, label_type = "pruned.labels", digits = 2, ...) {
    pull(x, label_type) %>%
        cluster_size(digits = digits) %>%
        kable_sc(digits = digits, ...)
}

#' @export
plot_res <- function(x, clusters, cols = palette_continuous(), ncol = 4, ...) {
    p_cl <- list.map(
        clusters,
        f(i) ~
            plot_dim(
                x,
                group.by = i,
                cols = cols(x@meta.data[, i] %>% unique() %>% length()),
                ...
            ) +
            NoLegend() +
            labs(title = str_remove(i, ".*(_w)?snn_res."))
    )
    p1 <- plot_grid(plotlist = p_cl[seq(8)], align = "hv", ncol = ncol, nrow = 2)
    print(p1)
    if (length(p_cl) > 8)
    plot_grid(plotlist = p_cl[8:15], align = "hv", ncol = ncol, nrow = 2)
}


#' @export
plot_module <- function(x, markers, ncol = 3, cols = pal_discrete_sc, ...) {

    res <- AddModuleScore(
        object = x,
        features = markers,
        name = names(markers)
    )

    score_cols <- paste0(names(markers), seq_along(names(markers)))
    df_medians <- res[[]] %>%
        summarise(across(all_of(score_cols), \(x) median(x, na.rm = TRUE)))

    list.map(
        score_cols,
        f(i, j) ~ {
            plot_violin_sc(
                object = res,
                features = i,
                cols = cols,
                ...
            ) %>%
                pluck(1) +
                labs(title = names(markers)[j]) +
                geom_hline(
                    yintercept = df_medians[[j]],
                    linetype = "dashed",
                    linewidth = 1
                )
        }
    ) %>%
    plot_grid(
        plotlist = .,
        ncol  = ncol,
        align = "hv"
    )
}

#' @export
plot_markers <- function(x, markers, cols = pal_discrete_sc) {
    p1 <- plot_mviolin(x, markers, cols = cols)
    print(p1)
    list.map(
        markers,
        f(i, j, k) ~ plot_dot(
            x,
            features = func_format(i) %>% head(50) %>% unique()
        ) +
            labs(title = k)
    )
}

#' @export
cluster_size <- function(x, order = TRUE, digits = 1) {
    res <- fct_count(x) %>%
        mutate(freq = round(n/sum(n) * 100, digits)) %>%
        filter(!is.na(f)) %>%
        column_to_rownames("f")

    if (order)
        arrange(res, desc(n))
    else
        res
}


#' @export
kable_stats <- function(x, digits = 3, file = NULL, ...) {
    p <- print_sc_stats(x, ...) %>%
        kable0(digits = digits)
    if (!is.null(file))
        save_kable(
            p,
            file = file,
            zoom = 2
        )
    print(p)
}

#' @export
print_sc_stats <- function(
        x,
        assay = "RNA",
        features = c(paste0("nFeature", "_", assay), paste0("nCount", "_", assay), "percent_mitochondrial", "percent_ribosomal"),
        probs = c(0, .025, .05, .1, seq(.25, .75, .25), .9, .95, .975, 1)
) {
    sapply(
        features,
        function(i) quantile(
            x[[i]],
            probs = probs,
            na.rm = TRUE)
    ) %>%
        set_colnames(colnames(.) %>% str_clean(assay))
}

#' @export
print_sc_stats2 <- function(
        x,
        y = x@meta.data,
        features = "percent_ribosomal",
        idents = levels(Idents(x)),
        probs = c(0, .025, .05, .1, seq(.25, .75, .25), .9, .95, .975, 1)
) {
    sapply(
        idents,
        function(i) quantile(
            unlist(y[WhichCells(x, idents = i), features]),
            probs = probs,
            na.rm = TRUE)
    )
}

#' @export
print_filtered_cells <- function(before, after) {
    before <- ncol(before)
    after <- ncol(after)
    pct <- ((after / before) * 100 - 100) %>% `*`(-1) %>% round(1) %>% paste0("%")
    paste0("Before: ", before, " cells / After: ", after, " cells (", pct, " removed).")
}

#' @export
integrate_multisamples <- function(x, labels = names(x)) {
    res <- merge(
        x = x[[1]],
        y = x[-1],
        add.cell.ids = labels,
        merge.dr = TRUE
    )

    res$sample <- rownames(res[[]])
    res[[]] <- separate(
        res[[]],
        col = "sample",
        into = c("Type", "Patient", "Barcode"),
        sep = "_"
    )
    res$sample <- NULL
    res$Patient <- paste0(res$Type, res$Patient)

    Idents(res) <- "Patient"
    return(res)
}

#' @export
calc_other_cluster_avg <- function(avg_matrix) {
    all_clusters <- colnames(select(avg_matrix, -c("gene", starts_with("except"))))

    for (target_cluster in all_clusters) {
        other_clusters <- setdiff(all_clusters, target_cluster)
        col_name <- paste0("except_", target_cluster)
        avg_matrix[[col_name]] <- rowMeans(avg_matrix[, other_clusters, drop = FALSE], na.rm = TRUE)
    }

    return(avg_matrix)
}

#' @export
pct_by_ident_type <- function(
        x,
        clusters = Idents(x),
        group.by = "Type"
) {
    expr_bin <- GetAssayData(x, layer = "data") > 0
    group <- interaction(clusters, x[[]][, group.by], drop = TRUE)
    sapply(
        levels(group),
        function(g) {
            cells <- colnames(x)[group == g]
            expr_bin[, cells, drop = FALSE] %>%
                as.matrix() %>%
                rowMeans()
        }) %>%
        as.data.frame()
}

#' @export
dea_sc <- function(
        x,
        ids = NULL,
        group.by = "ident",
        grouping.var = "Type",
        assay = "RNA"
    ) {
    cts <- AverageExpression(
        x,
        group.by = group.by,
        slot = ifelse(assay == "SCT", "data", "counts"),
        return.seurat = FALSE
    ) %>%
        pluck(assay) %>%
        as.data.frame() %>%
        rename_with(~ str_remove(., "^g")) %>%
        mutate(gene = rownames(.)) %>%
        calc_other_cluster_avg()

    if (!is.null(ids)) {
        it <- seq(nrow(ids))
    } else {
        it <- levels(Idents(x))
    }
    cluster_marker_genes0 <- list.map(
        it,
        f(i) ~ {
            if (is.null(grouping.var)) {
                grouping.var <- "orig.ident"
            }
            if (!is.null(ids)) {
                ident.1 <- ids[i, 1]
                ident.2 <- ids[i, 2]
                cluster <- paste(ident.1, ident.2, sep = " vs ")
            } else {
                ident.1 <- i
                ident.2 <- NULL
                cluster <- i
            }
            res <- FindConservedMarkers(
                ident.1 = ident.1,
                ident.2 = ident.2,
                grouping.var = grouping.var,
                object = x,
                logfc.threshold = .Machine$double.xmin,
                min.pct = .Machine$double.xmin,
                assay = assay
            ) %>%
            mutate(
                cluster = cluster,
                gene = rownames(.),
                log2FoldChange = rowMeans(select(., ends_with("avg_log2FC")), na.rm = TRUE),
                padj = rowMeans(select(., ends_with("p_val_adj")), na.rm = TRUE),
                pval = rowMeans(select(., ends_with("p_val")), na.rm = TRUE),
                `pct.1` = rowMeans(select(., ends_with("pct.1")), na.rm = TRUE),
                `pct.2` = rowMeans(select(., ends_with("pct.2")), na.rm = TRUE)
            ) %>%
                select(-contains(unique(x$orig.ident)[1]))

            if (!is.null(ids)) {
                res <- left_join(
                    res,
                    select(cts, all_of(c("gene", ids[i, 1], ids[i, 2]))) %>%
                        set_colnames(c("gene", "exp.1", "exp.2")),
                    by = "gene")
            } else {
                res <- left_join(
                    res,
                    select(cts, all_of(c("gene", i, paste0("except_", i)))) %>%
                        set_colnames(c("gene", "exp.1", "exp.2")),
                    by = "gene")
            }
        }
    )
    cluster_marker_genes <- list.rbind(cluster_marker_genes0)

    min_p <- cluster_marker_genes %>%
        filter(padj > 0) %>%
        summarise(min_padj = min(padj)) %>%
        pull(min_padj)

    cluster_marker_genes$padj[which(cluster_marker_genes$padj == 0)] <- min_p

    return(cluster_marker_genes)
}

#' @export
kable_dea <- function(x, file = NULL) {
    p <- list.map(
        x,
        f(i, ii, iii) ~
            fct_count(i$Expression) %>%
            filter(f!="ns") %>%
            set_colnames(c("Nb DEG", iii))) %>%
    Reduce(function(i, j) left_join(i, j, by = "Nb DEG"), .) %>%
    bind_rows(
        summarise(
              .,
              `Nb DEG` = "Total",
              across(-`Nb DEG`, \(x) sum(x, na.rm = TRUE))
            )
        ) %>%
    kable0()
    if (!is.null(file))
        save_kable(
            p,
            file = file,
            zoom = 2
        )
    p
}

#' @export
dea2deseq <- function(
        x,
        pct_threshold = .25,
        exp_threshold = .1
    ) {
    cls <- pull(x, "cluster") %>% unique()
    is_id2 <- str_detect(cls, " vs ") %>% any()
    if (is_id2) {
        func <- identity
    } else {
        func <- as.numeric
    }
    x %>%
    mutate(cluster = factor(cluster, levels = as.character(func(cls)))
        ) %>%
    group_split(cluster, .keep = TRUE) %>%
    set_names(map_chr(., ~as.character(first(.x$cluster)))) %>%
    list.map(
        select(., -cluster) %>%
            mutate(
                Expression = ifelse(
                    !((pct.1 < pct_threshold | exp.1 < exp_threshold) & (pct.2 < pct_threshold | exp.2 < exp_threshold)),
                    Expression,
                    "ns"
                ),
                pfc = -log10(padj) * log2FoldChange,
                rank_p = dense_rank(log10(padj)),
                rank_fc = dense_rank(desc(abs(log2FoldChange))),
                rank_pfc = dense_rank((rank_p + rank_fc) / 3)
            )
        )
}

#' @export
volcano_sc <- function(x, top_genes, fc_threshold = log2(1.5), p_threshold = .05, ...) {
    list.map(
        discard(x, function(i) nrow(i) == 0),
        f(i, j, k) ~ {
            volcano_plot(
                i,
                title = k,
                top_genes = top_genes[[j]],
                fc_threshold = fc_threshold,
                cex = .75,
                cex_genes = 4,
                p_threshold = p_threshold,
                max.overlaps = 1e3,
                force = 1e2,
                ...
            ) + labs(x = NULL, y = NULL, title = k) +
                theme(legend.position = "none")
        })
}

#' @export
plot_mqc <- function(
        x,
        features = c("nFeature_RNA", "nCount_RNA", "percent_mitochondrial", "percent_ribosomal", "percent_hemoglobin"),
        file = NULL,
        nrow = 2,
        digits = 3,
        ...
    ) {
    kable_stats(x, features = features, probs = c(0, 0.1, 0.5, 0.9, 1), digits = digits, file = file)

    p <- list.map(
        features[c(1, 3)],
        f(i, j) ~ {
            density_scatter(
                x,
                x = i,
                y = features[2],
                log_x = TRUE,
                log_y = TRUE
            )
        }
    ) %>% plot_grid(plotlist = ., align = "hv")
    print(p)

    Idents(x) <- "Patient"
    plot_sc_violin(
        x,
        features = features,
        normalize = rep(10, length(features)),
        nrow = 2,
        breaks = list(
            10^seq(-1, 10, by = 1),
            10^seq(-1, 10, by = 1),
            10^seq(-1, 3, by = 1),
            10^seq(-1, 3, by = 1),
            10^seq(-1, 3, by = 1)
        ),
        ...
    )
}

#' @export
plot_mqc_atac <- function(x, colour = palette_discrete(), nrow = 2, digits = 3, file = NULL, ...) {
    kable_stats(x, features = c(features_atac, features_supp), probs = c(0, 0.1, 0.5, 0.9, 1), digits = digits, file = file)

    p <- list.map(
        c("nCount_ATAC", "percent_mitochondrial"),
        f(i) ~ {
            density_scatter(
                x,
                x = i,
                y = "nFeature_ATAC",
                assay = "ATAC",
                log_x = TRUE,
                log_y = TRUE
            )
        }
    ) %>% plot_grid(plotlist = ., align = "hv")
    print(p)

    Idents(x) <- "Patient"
    p <- plot_sc_violin(
        x,
        features = features_atac,
        assay = "ATAC",
        cols = colour,
        nrow = nrow,
        normalize = c(10, 10, 2, 10, 0, 0),
        breaks = list(
            10^seq(-1, 10, by = 1),
            10^seq(-1, 10, by = 1),
            c(.1, .5, 1, 2, 4),
            c(1, 5, 10, 25),
            c(0, 20, 40, 60, 80, 100),
            NULL
        ),
        ...
    )
    print(p)

    plot_sc_violin(
        x,
        features = features_supp,
        assay = "ATAC",
        cols = colour,
        nrow = nrow,
        normalize = c(2, 10, 2, 10),
        ...
    )
}

#' @export
subsampling_sc <- function(x, group.by = "Type") {
    meta <- x@meta.data
    meta$cell_id <- rownames(meta)

    cells_to_keep <- meta %>%
        group_by(!!sym(group.by)) %>%
        sample_n(size = min(table(meta[, group.by]), n()), replace = FALSE) %>%
        pull(cell_id)

    seurat_subset <- subset(x, cells = cells_to_keep)
}


#' @export
reorder_idents <- function(x, y = "predicted.labels") {
    x@meta.data[, y] <- factor(
        x@meta.data[, y],
        levels = 0:(length(unique(x@meta.data[, y])) - 1)
    )
    Idents(x) <- x@meta.data[, y]
    return(x)
}

#' @export
organise_plots <- function(x, i = 2) {
    n <- length(x)

    if (n == 4*i) return(x)

    bloc <- n / i
    missing <- 4 - bloc

    c(
        x[seq_len(bloc)],
        rep(list(NULL), missing),
        x[seq_len(bloc) + bloc],
        rep(list(NULL), missing)
    )
}



#' @export
plot_marker_split <- function(x, markers, n_line = nlevels(x), split.by = "Type", cols = palette_discrete(), ...) {
    p1 <- list.map(
        markers,
        f(i, j, k) ~ {
            plot_dot(x, head(i, 50) %>% unique(), split.by = split.by, cols = cols) +
            geom_hline(yintercept = seq(2, (n_line - 1) * 2, by = 2) + 0.5) +
            ggtitle(k)
        }
    )
    print(p1)

    plot_mviolin(x, markers, split.by = split.by, cols = cols, ...)
}

#' @export
plot_feature_split <- function(x, markers, split.by = "Type", cols = palette_discrete(), ...) {
    list.map(
        markers,
        f(i, j, k) ~ {
            head(i, 12) %>% unique() %>%
                split(ceiling(seq_along(.)/4)) %>%
                list.map(
                    f(it) ~ {
                        plot_feature(
                            x,
                            features = it,
                            split.by = split.by,
                            cols = brewer.pal(9, "Reds"),
                            pt.size = .5,
                            ...
                        ) %>%
                            organise_plots() %>%
                            theme_multiple(
                                ncol = NULL,
                                nrow = 2,
                                title = k
                            )
                    }
                )
        }
    )
}

#' @export
reorder_celltype <- function(x) {
    b_cell <- c("B.FRF", "B.FRE", "preB.FRD", "preB.FRC", "proB.FRBC", "proB.FRA", "proB.CLP")

    all_levels <- unique(x)
    b_cell <- intersect(b_cell, all_levels)
    other_levels <- setdiff(all_levels, b_cell)

    sc <- grep("^SC\\.", other_levels, value = TRUE)
    dc <- grep("^DC\\.", other_levels, value = TRUE)

    other <- setdiff(other_levels, c(sc, dc))

    c(
        b_cell,
        sort(sc),
        sort(dc),
        sort(other)
    )
}

#' @export
format_celltype <- function(
        x,
        annotation,
        type = "cell_type",
        label_type = "pruned.labels",
        lim = 50,
        path_fig = NULL
    ) {
    table_annotation <- annotation[[type]] %>%
        as.data.frame()
    x@meta.data[, paste0(type, "_raw")] <- pull(table_annotation, label_type)

    if (!is.null(path_fig))
    kable_ncell(
        table_annotation,
        label_type,
        digits = 2,
        file =  file.path(path_fig, paste("table_", type, "_before.png"))
    )
    p2 <- pull(table_annotation, label_type) %>%
        plot_bar_sc()
    print(p2)

    if (!is.null(path_fig))
        kable_ncell(
            table_annotation,
            label_type,
            digits = 2,
            file =  file.path(path_fig, paste0("table_", type, "_after.png"))
        )

    format_celltype0(x, table_annotation, type, label_type, lim)
}

#' @export
format_celltype0 <- function(
        x,
        table_annotation = x@meta.data,
        type = "cell_type",
        label_type = paste0(type, "_raw"),
        lim = 50
) {
    to_remove <- table(pull(table_annotation, label_type)) %>%
        .[. < lim] %>%
        names()
    table_annotation[pull(table_annotation, label_type) %in% to_remove, label_type] <- NA

    x@meta.data[, type] <-
        table_annotation[match(rownames(x[[]]), rownames(table_annotation)), label_type] %>%
        factor(levels = reorder_celltype(.))

    return(x)
}

#' @export
set_percent <- function(x) {
    x[["percent_mitochondrial"]] <- PercentageFeatureSet(x, pattern = "^((MT)|(mt))-")
    x[["percent_ribosomal"]] <- PercentageFeatureSet(x, pattern = "^((RP[LS])|(Rp[ls]))")
    x[["percent_hemoglobin"]] <- PercentageFeatureSet(x, pattern = "^((^HB[^(P)])|(Hb[^(p)]))")
    return(x)
}

#' @export
calculate_doublet_by_condition <- function(x, split_by = "Patient", total = TRUE, label = "scDblFinder") {
    res <- list.map(
        unique(x[[]][, split_by]) %>% as.character() %>% sort(),
        f(i) ~ filter(x[[]], !!sym(split_by) == i) %>% pull(all_of(label)) %>% table()
    ) %>%
        list.rbind() %>%
        as.data.frame() %>%
        set_colnames(colnames(.) %>% str_to_sentence())

    if (isTRUE(total)) {
        res <- res %>%
            rbind(".Total" = colSums(.))
    }
    return(res)
}

#' @export
kable_doublet <- function(x, file = NULL, ...) {
    p <- calculate_doublet_by_condition(x, ...) %>%
        kable0()

    if (!is.null(file))
        save_kable(
            p,
            file = file,
            zoom = 2
        )
    print(p)
}

#' @export
plot_bar_doublet <- function(x, normalize = FALSE, pct = TRUE, ...) {
    res <- calculate_doublet_by_condition(x, ...) %>%
        t() %>%
        as.data.frame()

    if (isTRUE(normalize)) {
        res <- res %>%
            mutate(across(everything(), ~ .x / sum(.x) * 100))
    }
    plot_bar_2cat(
        res,
        count = TRUE,
        stats = FALSE,
        pct = pct,
        colour = palette_discrete()[c(1, 3)],
        colour_text = "black"
    )
}

#' @export
plot_bar_kept <- function(x, y, l_samples, ...) {
    list.map(
        l_samples,
        f(i) ~ list.map(
            c(x, y),
            f(j) ~ filter(j[[]], Patient == i) %>% nrow()
        ) %>%
            list.rbind()
    ) %>%
    list.cbind() %>%
    as.data.frame() %>%
    set_colnames(l_samples) %>%
    t() %>%
    as.data.frame() %>%
    rbind(`.Total` = colSums(.)) %>%
    mutate(V3 = V1 - V2) %>%
    select(-V1) %>%
    t() %>%
    as.data.frame() %>%
    set_rownames(c("Kept", "Removed"))  %>%
    plot_bar_2cat(
        count = TRUE,
        stats = FALSE,
        ...
    )
}

#' @export
peak_per_chr <- function(x, chr) {
    chrY_ranges <- GRanges(chr, IRanges(start = 1, end = 5e8))

    chrY_counts <- FeatureMatrix(
        fragments = Fragments(x),
        features = chrY_ranges,
        cells = colnames(x)
    )

    Matrix::colSums(chrY_counts)
}

#' @export
pct_peak_per_chr <- function(x, chr) {
    peaks <- granges(x)

    y_peaks <- peaks[seqnames(peaks) == chr]
    y_peak_names <- paste0(seqnames(y_peaks), "-", start(y_peaks), "-", end(y_peaks))

    atac_counts <- GetAssayData(x, assay = "ATAC", layer = "counts")
    atac_y <- atac_counts[y_peak_names, ]

    total_fragments <- Matrix::colSums(atac_counts)
    chrY_fragments <- Matrix::colSums(atac_y)

    (chrY_fragments / total_fragments) * 100
}

#' @export
downsample_by_cluster <- function(x, group_by = "Type") {

    df <- data.frame(
        cell = colnames(x),
        cluster = Idents(x),
        group_by = x[[group_by]][,1],
        stringsAsFactors = FALSE
    )

    selected_cells <- unlist(
        lapply(split(df, df$cluster), function(i) {

            counts <- table(i$group_by)
            n_min <- min(counts)

            unlist(
                lapply(names(counts), function(j) {
                    sample(i$cell[i$group_by == j], n_min)
                })
            )
        })
    )

    subset(x, cells = selected_cells)
}

#' @export
remove_low_cells <- function(x, label = "cell_subtype_formatted", n = 25) {
    cell_counts <- table(x@meta.data[, label])
    cell_types_low <- names(cell_counts[cell_counts < n])

    x@meta.data[, label][x@meta.data[, label] %in% cell_types_low] <- NA
    x@meta.data[, label] <- fct_drop(x@meta.data[, label])

    return(x)
}

#' @export

load_dataset <- function(path_data, seurat_dataset, split.by, target_type, pal_discrete_sc, pal_discrete_sc2, integrated = TRUE, snn_cluster) {
    load(file.path(path_data, paste0(seurat_dataset, ".rda")))
    Idents(seurat) <- snn_cluster
    if (integrated) {
        seurat_subset <- subsampling_sc(seurat) %>%
            subset(subset = Type == target_type)
    } else {
        seurat_subset <- seurat
    }

    if (split.by != "Type") {
        if (integrated) {
            load(file.path(path_data, paste0(seurat_dataset, "_", target_type, ".rda")))
            Idents(seurat) <- snn_cluster
        } else {
            load(file.path(path_data, paste0(seurat_dataset, ".rda")))
            Idents(seurat) <- snn_cluster
        }
        seurat_subset <- subset(seurat, cells = Cells(seurat_subset))
        n_idents <- nlevels(Idents(seurat))
        if (n_idents > length(pal_discrete_sc))
            pal_discrete_sc <- palette_continuous()(n_idents)
    }
    seurat_subset$cell_subtype_formatted <- fct_drop(seurat_subset$cell_subtype_formatted)

    if (split.by == "Type") {
        seurat <- subset(seurat, subset = Type == target_type)
    }
    seurat$cell_subtype_formatted <- fct_drop(seurat$cell_subtype_formatted)

    seurat <- remove_low_cells(seurat) %>%
        remove_low_cells("cell_subtype_formatted2")

    if (split.by != "Type") {
        pal_discrete_sc2 <- palette_continuous()(length(unique(seurat$cell_subtype_formatted)))
    }
    list(seurat = seurat, seurat_subset = seurat_subset, pal_discrete_sc = pal_discrete_sc, pal_discrete_sc2 = pal_discrete_sc2)
}

#' @export
plot_feature_soupex <- function(x, markers, tsne_matrix, assay = "RNA") {

    DefaultAssay(x) <- assay
    x[["tsne"]] <-  CreateDimReducObject(
        embeddings = tsne_matrix,
        key = "tSNE_",
        assay = assay
    )

    list.map(
        unique(x$Patient),
        f(i) ~ {
            WhichCells(x, expression = Patient == i) %>%
                plot_feature(
                    x,
                    markers,
                    cells = .,
                    pt.size = .25,
                    cols = brewer.pal(9, "Reds"),
                    slot = "counts"
                ) %>%
                plot_grid(plotlist = ., ncol = 3) %>%
                plot_grid(
                    ggdraw() +
                        draw_label(
                            i,
                            color = "#2F5496",
                            size = 20,
                            fontface = "bold.italic"
                        ),
                    .,
                    ncol = 1,
                    rel_heights = c(0.1, 1)
                )
        }
    )
}

#' @export
format_gene <- function(x) {
    str_to_lower(x) %>%
        str_replace_all("(^|\\-)([a-z])", function(x) toupper(x))
}

#' @export
inverser_vs <- function(texte) {
    str_replace(texte, "(.*) vs (.*)", "\\2 vs \\1")
}

#' @export
filter_depth_cor <- function(object, cor.cutoff = 0.5) {
    embeddings <- Embeddings(object, reduction = "lsi")
    cor.values <- apply(embeddings, 2, function(x) cor(x, object$nCount_ATAC, method = "spearman"))
    keep.mask <- abs(cor.values) <= cor.cutoff
    as.numeric(gsub("[^0-9]", "", names(cor.values)[keep.mask]))
}

#' @export
print_anchors <- function(x, n_anchors, rna, atac, ...) {
    data.frame(
        n_anchors = nrow(x),
        pct_anchors = nrow(anchor_df) / n_anchors * 100,
        pct_rna = 100 * length(unique(x[, 1])) / ncol(rna),
        pct_atac = 100 * length(unique(x[, 2])) / ncol(atac)
    ) %>%
        kable_sc(...)
}

#' @export
add_dr <- function(x, sample, assay = "atac", dr = "umap") {
    cls <- file.path(path_project, sample, "outs", "analysis", "clustering", assay) %>%
        list.files(full.names = TRUE) %>%
        .[str_detect(., "graphclust$")] %>%
        list.files(full.names = TRUE) %>%
        .[str_detect(., "clusters.csv$")] %>%
        read_csv() %>%
        rename_with(~ paste0("cluster_", assay), .cols = "Cluster")
    meta_data <- file.path(path_project, sample, "outs", "analysis", "dimensionality_reduction", assay, paste0(dr, "_projection.csv")) %>%
        read_csv() %>%
        left_join(cls, by = "Barcode") %>%
        filter(Barcode %in% Cells(x)) %>%
        column_to_rownames("Barcode")
    x <- AddMetaData(x, metadata = meta_data)
    x[[paste0(dr, "_before_", assay)]] <- meta_data %>%
        .[, paste0(toupper(dr), "-", seq(2))] %>%
        as.matrix() %>%
        CreateDimReducObject(
            embeddings = .,
            key = paste0(toupper(dr), "_before_", assay, "_"),
            assay = DefaultAssay(x)
        )
    return(x)
}

#' @export
find_doublets <- function(seurat, assay = assay) {
    if (.Platform$OS.type == "windows") {
        bp_param <- SnowParam(workers = n_cores, type = "SOCK",  exportglobals = TRUE, progressbar = TRUE)
    } else {
        bp_param <- MulticoreParam(workers = n_cores, progressbar = TRUE)
    }

    sce <- scDblFinder(
        GetAssayData(seurat[[assay]], layer = "counts"),
        samples = seurat$Patient,
        aggregateFeatures = assay == "ATAC",
        nfeatures = ifelse(assay == "ATAC", 30, 1352),
        k = if (assay == "ATAC") 20 else NULL,
        processing = "normFeatures",
        dims = n_dim,
        cluster = as.factor(seurat[[]][, paste0("cluster_", str_to_lower(ifelse(assay == "ATAC", assay, "gex")))]),
        BPPARAM = bp_param,
        artificialDoublets = if (assay == "ATAC") 1 else NULL
    )

    seurat$scDblFinder <- sce[, row.names(seurat[[]])]$scDblFinder.class
    seurat$doublet_score <- sce[, row.names(seurat[[]])]$scDblFinder.score

    return(seurat)
}

#' @export
plot_qc_doublet <- function(seurat, cols = col_inds) {
    kable_doublet(seurat, file = file.path(path_fig, "table_doublets_scDblFinder.png"))

    plot_bar_doublet(seurat) %>% plot()

    plot_violin_doublet <- function(x, ...) plot_sc_violin(x, features = "doublet_score", ncol = 1, nrow = 1, ...)
    Idents(seurat) <- "scDblFinder"
    p1 <- plot_violin_doublet(seurat)
    Idents(seurat) <- "Patient"
    p1 | plot_violin_doublet(seurat, cols = cols)
}

#' @export
plot_dim_doublet <- function(seurat) {
    selected <- WhichCells(seurat, expression = scDblFinder == "doublet")
    list.map(
        unique(seurat$Patient),
        f(x) ~ {
            cells <- WhichCells(seurat, expression = Patient == x)
            intersect(selected, cells) %>%
                plot_dim_selected(seurat, selected = ., cells = cells) +
                ggtitle(x) +
                scale_color_manual(
                    labels = c("Singlets", "Doublets"),
                    values = c("grey80", "red")
                )
        }) %>%
        plot_grid(plotlist = ., ncol = 2)
}

#' @export
plot_feature_qc2 <- function(x, features, assay = "RNA") {
    list.map(
        unique(seurat$Patient),
        f(x) ~ {
            plot_feature_qc(
                seurat,
                features,
                assay,
                cells = WhichCells(seurat, expression = Patient == x)
                ) %>%
                plot_grid(plotlist = p, nrow = 2, align = "hv")
        }
    )
}

#' @export
plot_feature_qc <- function(x, features, assay = "RNA", ...) {
    map(
        features,
        ~FeaturePlot(
            object = seurat,
            features = .,
            alpha = 0.5,
            cols = c("yellow", "blue"),
            ...
        ) %>%
            theme_sc_dim() +
            ggtitle(str_clean(., assay))
    )
}
