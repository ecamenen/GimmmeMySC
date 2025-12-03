#' @export
plot_dim_cycle <- function(x, ...) {
    Idents(x) <- "Phase"
    p1 <- plot_feature(
        x,
        features = "CC.Difference",
        pt.size = .5,
        ...
    ) %>% pluck(1) +
        labs(title = NULL, color = "S to G2/M")

    p2 <- plot_dim(
        object = x,
        label = FALSE,
        title = TRUE,
        ...
    ) +  labs(x = NULL, y = NULL)
    p1 + p2
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
        seurat,
        clusters,
        file = NULL,
        assay = "RNA"
    ) {
    cl_size <- map_dfr(
        clusters,
        ~ cluster_size(seurat@meta.data[, .x]) %>%
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
        mutate(Clusters = str_remove_all(Clusters, paste0(assay, "_snn_res.")))

    cl_size[is.na(cl_size)] <- ""
    p <- kable0(cl_size)

    if (!is.null(file))
        save_kable(
            p,
            file = file,
            zoom = 2
        )
    p
}

#' @export
kable_sc <- function(x, func = str_clean, digits = getOption("digits"), file = NULL) {
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
plot_res <- function(seurat, clusters, cols = palette_continuous(), ncol = 4, ...) {
    p_cl <- list.map(
        clusters,
        f(x, y, z) ~
            plot_dim(
                seurat,
                group.by = x,
                cols = cols(seurat@meta.data[, x] %>% unique() %>% length()),
                ...
            ) +
            NoLegend() +
            labs(title = str_remove(x, ".*_snn_res."))
    )
    p1 <- plot_grid(plotlist = p_cl[seq(8)], align = "hv", ncol = ncol, nrow = 2)
    print(p1)
    if (length(p_cl) > 8)
    plot_grid(plotlist = p_cl[8:15], align = "hv", ncol = ncol, nrow = 2)
}


#' @export
plot_module <- function(x, markers, ncol = 3, cols = pal_discrete_sc, ...) {
    list.map(
        markers,
        f(i, j, k) ~ {
            tmp <- list(i) %>%
                AddModuleScore(x, features = .)
            stats <- quantile(tmp$Cluster1, probs = c(0.01, 0.99))
            tmp <- subset(tmp, subset = Cluster1 >= stats[1] & Cluster1 <= stats[2])
            plot_violin_sc(
                object = tmp,
                features = "Cluster1",
                cols = cols,
                ...
            ) %>%
                pluck(1) +
                labs(title = str_remove_all(k, " ?vs_all"))
        }
    ) %>%
    plot_grid(
        plotlist = .,
        ncol  = ncol,
        align = "hv"
    )
}

#' @export
plot_markers <- function(x, markers) {
    p1 <- plot_mfeature(x, markers, nrow = 3)
    print(p1)
    p2 <- plot_mviolin(x, markers, cols = pal_discrete_sc)
    print(p2)
}

#' @export
plot_markers2 <- function(x, markers) {
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
cluster_size <- function(x, order = TRUE) {
    res <- fct_count(x) %>%
        filter(!is.na(f)) %>%
        column_to_rownames("f") %>%
        mutate(freq = round(n/sum(n) * 100, 1))

    if(order)
        arrange(res, desc(n))
    else
        res
}

#' @export
print_sc_stats <- function(
        x,
        assay = "RNA",
        features = c(paste0("nFeature", "_", assay), paste0("nCount", "_", assay), "percent_mitochondrial"),
        probs = c(0, .025, .05, .1, seq(.25, .75, .25), .9, .95, .975, 1)
) {
    sapply(
        features,
        function(i) quantile(
            x[[i]],
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
integrate_multisamples <- function(
        x,
        labels = paste0(
            types,
            "_",
            # str_remove_all(l_samples, types) %>% str_remove_all("_")
            l_samples
        )
) {
    seurat <- merge(
        x = x[[1]],
        y = x[-1],
        add.cell.ids = labels
    )

    seurat$sample <- rownames(seurat[[]])
    seurat[[]] <- separate(
        seurat[[]],
        col = "sample",
        into = c("Type", "Patient", "Barcode"),
        sep = "_"
    )
    return(seurat)
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
        seurat,
        clusters = Idents(seurat),
        group.by = seurat$Type
) {

    expr_bin <- GetAssayData(seurat, layer = "data") > 0
    group <- interaction(clusters, group.by, drop = TRUE)
    percent_matrix <- sapply(
        levels(group),
        function(g) {
            cells <- colnames(seurat)[group == g]
            rowMeans(expr_bin[, cells, drop = FALSE])
        })
    as.matrix(percent_matrix) %>% as.data.frame()
}

#' @export
dea_sc <- function(
        seurat,
        ids = NULL,
        group.by = "ident",
        grouping.var = "Patient",
        assay = "RNA"
    ) {
    cts <- AverageExpression(
        seurat,
        group.by = group.by,
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
        it <- levels(Idents(seurat))
    }
    cluster_marker_genes0 <- list.map(
        it,
        f(x) ~ {
            if (!is.null(ids)) {
                ident.1 = ids[x, 1]
                ident.2 = ids[x, 2]
                cluster = paste(ident.1, ident.2, sep = " vs ")
            } else {
                ident.1 = x
                ident.2 = NULL
                cluster = x
            }
            res <- FindConservedMarkers(
                ident.1 = ident.1,
                ident.2 = ident.2,
                grouping.var = grouping.var,
                object = seurat,
                logfc.threshold = .Machine$double.xmin
            ) %>%
            mutate(
                cluster = cluster,
                gene = rownames(.),
                log2FoldChange = rowMeans(select(., ends_with("avg_log2FC")), na.rm = TRUE),
                padj = rowMeans(select(., ends_with("p_val_adj")), na.rm = TRUE),
                `pct.1` = rowMeans(select(., ends_with("pct.1")), na.rm = TRUE),
                `pct.2` = rowMeans(select(., ends_with("pct.2")), na.rm = TRUE)
            )

            if (!is.null(ids)) {
                res <- left_join(
                    res,
                    select(cts, all_of(c("gene", ids[x, 1], ids[x, 2]))) %>%
                        set_colnames(c("gene", "exp.1", "exp.2")),
                    by = "gene")
            } else {
                res <- left_join(
                    res,
                    select(cts, all_of(c("gene", x, paste0("except_", x)))) %>%
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
    is_id2 <- str_detect(cls, " vs ")
    if (is_id2) {
        func <- identity
    } else {
        func <- as.numeric
    }
    res <- x %>%
    mutate(cluster = factor(cluster, levels = as.character(sort(func(cls))))
        ) %>%
    group_split(cluster, .keep = TRUE) %>%
    set_names(map_chr(., ~as.character(first(.x$cluster)))) %>%
    list.map(
        select(., -cluster) %>%
            mutate(
                Expression = ifelse(
                    pct.1 > pct_threshold | exp.1 > exp_threshold,
                    Expression,
                    "ns"
                )
            )
    )

    if (is_id2) {
        list.map(
            res,
                mutate(
                    .,
                    Expression = ifelse(
                        pct.2 > pct_threshold | exp.2 > exp_threshold,
                        Expression,
                        "ns"
                    )
                )
        )
    } else {
        res
    }
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
            ) + labs(x = NULL, y = NULL) +
                theme(legend.position = "none")
        })
}


#' @export
format_dea_sc <- function(
        x,
        pct_threshold = .25,
        exp_threshold = .1,
        fc_threshold = log2(1.5),
        p_threshold = .05
    ) {
    is_id2 <- pull(x, "cluster") %>%
        unique() %>%
        str_detect(" vs ")

    list.map(
        unique(x$cluster),
        f(i) ~ {
            res <- filter(x, cluster == i) %>%
                filter(pct.1 > pct_threshold |
                           exp.1 > exp_threshold)

            if (is_id2) {
                res <- filter(
                    res,
                    pct.2 > pct_threshold |
                        exp.2 > exp_threshold
                )
            }

            top_genes(
                res,
                return_rank = TRUE,
                fc_threshold = fc_threshold,
                p_threshold = p_threshold
            )
        }
    ) %>%
        discard(~ nrow(.x) == 0)
}

#' @export
plot_mqc <- function(
        x,
        features = c("nFeature_RNA", "nCount_RNA", "percent_mitochondrial", "percent_ribosomal"),
        file = NULL
    ) {
    p <- print_sc_stats(x, probs = c(0, 0.1, 0.5, 0.9, 1)) %>%
        kable_sc(digits = 3)
    if (!is.null(file))
        save_kable(
            p,
            file = file,
            zoom = 2
        )
    print(p)

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
        cols = col_inds,
        nrow = 1,
        normalize = c(10, 10, 10, 10)
    )
}

#' @export
subsampling_sc <- function(x, group.by = "Patient") {
    meta <- seurat@meta.data
    meta$cell_id <- rownames(meta)

    cells_to_keep <- meta %>%
        group_by(!!sym(group.by)) %>%
        sample_n(size = min(table(meta[, group.by]), n()), replace = FALSE) %>%
        pull(cell_id)

    seurat_subset <- subset(seurat, cells = cells_to_keep)
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
plot_marker_split <- function(x, markers, n_line = nlevels(x), split.by = "Type") {
    p1 <- list.map(
        markers,
        f(i, j, k) ~ {
            plot_dot(x, head(i, 50) %>% unique(), split.by = split.by) +
            geom_hline(yintercept = seq(2, (n_line - 1) * 2, by = 2) + 0.5) +
            ggtitle(k)
        }
    )
    print(p1)

    p2 <- list.map(
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
                            pt.size = .5
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
    print(p2)

    plot_mviolin(x, markers, split.by = split.by)
}
