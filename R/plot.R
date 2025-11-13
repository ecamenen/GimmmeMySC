#' @export
plot_sc_violin <- function(
        x,
        features = c(paste0("nFeature", "_", assay), paste0("nCount", "_", assay), "percent_mitochondrial"),
        normalize = NULL,
        nrow = 3,
        ncol = 3,
        cols = brewer.pal(9, "Set1"),
        width = 15,
        breaks = NULL,
        metadata = TRUE,
        ...
) {
    if(is.null(normalize)) {
        normalize <- rep(TRUE, length(features))
    }
    list.map(
        features,
        f(i, j) ~{
            p <- VlnPlot(
                x,
                features = i,
                pt.size = 0,
                # group.by = "Patient",
                cols = cols,
                # cols = cols[j],
                ncol = ncol,
                ...
            ) +
                theme_custom() +
                theme(
                    axis.title.x = element_blank(),
                    axis.line = element_line(linewidth = 1),
                    panel.grid.major.y = element_line(colour = "grey", linetype = 2)
                )

            if (normalize[j] > 0) {
                if (!metadata) {x
                    df <- t(x[[x@active.assay]]@data[i, ])
                } else {
                    df <- x@meta.data[, i]
                }
                p <- p + GimmeMyPlot:::axis_log(df, "y", normalize[j], breaks[[j]])
            } else {
                if (!is.null(breaks[[j]]))
                    p <- p + scale_y_continuous(labels = label_number_auto(), breaks = breaks[[j]])
                else
                    p <- p + scale_y_continuous(labels = label_number_auto())
            }
            p +
                NoLegend()  +
                ggtitle(
                    str_clean(i) %>%
                        str_wrap(width = width)
                )  +
                ylab(NULL) +
                # theme(axis.text.x = element_blank()) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
    ) %>% wrap_plots(nrow = nrow)
}

#' @export
print_sc_stats <- function(
        x,
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
plot_sc_cor <- function(
        x,
        col = palette_discrete()[seq_along(Idents(x))],
        features = c("percent_mitochondrial", paste0("nFeature", "_", assay))
) {
    lapply(
        features,
        function(i) {
            FeatureScatter(
                x,
                feature1 = paste0("nCount", "_", assay),
                feature2 = i,
                cols = col
            ) +
                geom_smooth(method  = "lm", se= TRUE, colour = "gray30") +
                NoLegend() +
                GimmeMyPlot:::theme_custom(cex = 1.5) +
                theme(
                    axis.line = element_line(linewidth = 1),
                    panel.grid.major = element_line(colour = "grey", linetype = 2),
                    axis.text = element_text(size = 13, color = "gray50")
                ) +
                labs(x = "# Count", y = str_clean(i))
        }
    ) %>%
        plot_grid(plotlist = ., ncol = 2)
}

#' @export
plot_eig <- function(x) {
    x %>%
        data.frame(
            dim = seq_along(.),
            var = .
        ) %>%
        ggplot(aes(x = dim)) +
        geom_line(aes(y = var), color = "red", lwd = 1) +
        geom_point(aes(y = var), color = "red", size = 2, pch = 3) +
        labs(
            title = NULL,
            x = "Dimension",
            y = "Standard deviation"
        ) +
        xlim(1, NA) +
        theme_custom()
}

#' @export
density_scatter <- function(object, x, y, ...) {
    DensityScatter(
        object,
        x,
        y,
        quantiles = TRUE,
        ...
    ) +
        theme_custom() +
        xlab(str_clean(x)) +
        ylab(str_clean(y))
}

#' @export
variable_plot <- function(x) {
    VariableFeaturePlot(x) %>%
        LabelPoints(
            points = head(VariableFeatures(x), 10),
            repel = TRUE
        ) +
        scale_x_continuous(
            trans = log10_trans(),
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
        ) +
        scale_y_continuous(
            trans = log2_trans(),
            breaks = trans_breaks("log2", function(x) 2^x),
            labels = label_number_auto()
        ) +
        theme_custom() +
        ylab("Variance") +
        NoLegend()
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
density_plot <- function(x, feature = nFeature_RNA) {
    as.data.frame(x[[]]) %>%
        ggplot(aes(color=orig.ident, x=feature, fill= orig.ident)) +
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_custom() +
        # geom_vline(xintercept = 300) +
        ylab(str_clean(i))
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
        add.cell.ids =labels
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
plot_feature <- function(
     object = object,
     features =  features,
     cols = pal_sc,
     reduction = "umap",
     pt.size = 1,
     alpha = .5,
     ...
    ) {
    FeaturePlot(
        object = object,
        features =  features,
        reduction = reduction,
        pt.size = pt.size,
        alpha = alpha,
        combine = FALSE,
        # min.cutoff = "q10",
        # max.cutoff = "q90",
        ...
    )  %>%
        map(
            ~theme_sc_dim(.) +
                labs(x = NULL, y = NULL) +
                scale_color_gradientn(colors = cols)
        ) %>%
        set_names(features)
}

#' @export
plot_mfeature <- function(
        object = object,
        features =  features,
        pt.size = .5,
        ncol = 4,
        nrow = NULL,
        func = function(x) func_format(x) %>% .[. %in% Features(object)] %>% head(12),
        cols = brewer.pal(9, "Reds"),
        ...
    ) {
    list.map(
        features,
        f(x, y, z) ~ {
            x <- unique(func(x))
            names(x) <- NULL
            plot_feature(
                object,
                features = x,
                pt.size = pt.size,
                cols = cols,
                ...
            ) %>%
            theme_multiple(ncol = ncol, nrow = nrow, title = z)
        }
    )
}

#' @export
plot_dot <- function(
        object = object,
        features = features,
        cols = NULL,
        split.by = NULL,
        ...
    ) {
    if (is.null(cols)) {
        if (!is.null(split.by)) {
            cols <- palette_discrete()
        } else {
            cols <- pal_sc
        }
    }
    p <- DotPlot(
        object = object,
        features = features,
        dot.scale = 8,
        assay = assay,
        cols = cols,
        split.by = split.by,
        # scale = FALSE,
        ...
    ) +
        RotatedAxis() +
        theme_custom() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = NULL, y = NULL, size = "Percentage Expressed") +
        guides(size = guide_legend(order = 1))
    if (is.null(split.by)) {
        p + scale_color_gradientn(colors = cols)
    } else {
        p
    }
}

#' @export
plot_violin_sc <- function(
        object = object,
        features =  features,
        pt.size = 0,
        cols = palette_discrete(),
        ...
    ) {
    VlnPlot(
        object = object,
        features =  features,
        pt.size = pt.size,
        cols = cols,
        combine = FALSE,
        ...
    ) %>%
        map(~. +theme_violin_sc() + labs(y = NULL) + NoLegend())
}

#' @export
plot_mviolin <- function(
        object = object,
        features =  features,
        func = function(x) func_format(x) %>% .[. %in% Features(object)] %>% head(12),
        ncol = 4,
        nrow = 3,
        ...
) {
    list.map(
        features,
    f(x, y, z) ~ {
        plot_violin_sc(
            object = object,
            features = func(x),
            ...
        ) %>%
        theme_multiple(ncol = ncol, nrow = nrow, title = z)
    }
    )
}

#' @export
plot_dim <- function(
        object = object,
        reduction = "umap",
        cols = palette_discrete(),
        label = TRUE,
        title = FALSE,
        axis = FALSE,
        ...
    ) {
        p <- DimPlot(
            object = object,
            reduction = reduction,
            cols = cols,
            pt.size = .5,
            alpha = .25,
            repel = TRUE,
            label = label,
            na.value = "white",
            ...
        ) %>%
            theme_sc_dim()
        if (isFALSE(title)) {
            p <- p + labs(title = NULL)
        }
        if (isFALSE(axis)) {
            p <- p + labs(x = NULL, y = NULL)
        }
        return(p)
}
