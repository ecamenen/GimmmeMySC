#' @export
plot_sc_violin <- function(
        x,
        assay = "RNA",
        features = c(paste0("nFeature", "_", assay), paste0("nCount", "_", assay), "percent_mitochondrial", "percent_ribosomal"),
        normalize = NULL,
        nrow = 3,
        ncol = 3,
        cols = palette_discrete(),
        width = 15,
        breaks = NULL,
        metadata = TRUE,
        ...
) {
    if(is.null(normalize)) {
        normalize <- rep(TRUE, length(features))
    }
    plot_list <- list.map(
        features,
        f(i, j) ~{
            p <- VlnPlot(
                x,
                features = i,
                pt.size = 0,
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
    )

    if (ncol * nrow < length(plot_list)) {
        ncol <- round(length(plot_list) / nrow)
    }
    grid.arrange(
        grobs = plot_list,
        nrow = nrow,
        ncol = ncol,
        padding = unit(0, "cm")
    )
}

#' @export
plot_sc_cor <- function(
        x,
        assay = "RNA",
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
                theme_custom(cex = 1.5) +
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
        theme_custom() +
        theme(
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
            panel.grid.minor = element_line(color = "grey80", linewidth = 0.15)
        )
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
plot_feature <- function(
     object = object,
     features =  features,
     cols = pal_sc,
     pt.size = 1,
     alpha = .5,
     limits = NULL,
     normalize = FALSE,
     split.by = NULL,
     ...
    ) {
    features <- features %>% .[features %in% Features(object) | features %in% colnames(object[[]])]
    if (isTRUE(normalize)) {
        cts <- GetAssayData(object = object)
        limits <- map(features, ~cts[., ] %>% quantile(c(0, 1)))
    } else {
        limits <- replicate(length(features), NULL, simplify = FALSE)
    }

    if(!is.null(split.by)) {
        n_groups <- length(unique(object[[]][, split.by]))
    } else {
        n_groups <- 1
    }

    limits <- set_names(limits, features)

    FeaturePlot(
        object = object,
        features =  features,
        pt.size = pt.size,
        alpha = alpha,
        combine = FALSE,
        split.by = split.by,
        ...
    ) %>%
        set_names(rep(features, n_groups)) %>%
        list.map(
            f(i, j, k) ~theme_sc_dim(i) +
                labs(x = NULL, y = NULL) +
                scale_color_gradientn(colors = cols, limits = limits[[k]])
        )
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
        # assay = assay,
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


plot_dot2 <- function(x, markers, ...) {
    last_row_numbers <- markers %>%
        filter(!duplicated(gene)) %>%
        mutate(row_num = row_number()) %>%
        group_by(cluster) %>%
        slice_tail(n = 1) %>%
        pull(row_num) %>%
        sort()

    plot_dot(
        x,
        features = unique(pull(markers, "gene")),
        ...
    ) +
        geom_vline(xintercept = last_row_numbers + 0.5)
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
        cols = palette_discrete(),
        label = TRUE,
        title = FALSE,
        axis = FALSE,
        pt.size = .5,
        alpha = .25,
        combine = TRUE,
        na.value = "white",
        ...
    ) {
        p <- DimPlot(
            object = object,
            cols = cols,
            pt.size = pt.size,
            alpha = alpha,
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

#' @export
plot_bar_sc <- function(x, title = NULL, ratio = 3, label_y = "percent", ...) {
    if (is(x, "Seurat")) {
        res <- Idents(x)
        sample_size <- ncol(x)
    } else {
        res <- x
        sample_size <- length(x)
    }
    res %>%
        cluster_size() %>%
        select(n) %>%
        plot_bar_sc0(
            label_y = label_y,
            sample_size = sample_size,
            ratio = ratio,
            title = title,
            ...
        )
}

#' @export
plot_bar_sc0 <- function(x, title = NULL, ratio = 3, label_y = "percent", ...) {
    plot_bar(
        x,
        label_x = "none",
        label_y = label_y,
        digits = 1,
        n_max = 15,
        ratio = ratio,
        ...
    ) + labs(title = title)
}

#' @export
clustree_sc <- function(x, prefix, clusters) {
    clustree(x, prefix = prefix) +
    scale_color_discrete(type = palette_continuous(gray = FALSE)(length(clusters))) +
    theme(
        legend.title = element_text(face = "italic", size = 15),
        legend.text = element_text(size = 10)
    )
}


#' @export
plot_ident_celltype <- function(x, cols = palette_discrete(), cell_label = "cell_subtype_formatted") {

    if ("0" %in% levels(Idents(x))) {
        legend <- reorder_celltype(x[[]][, cell_label])
    } else if ("0" %in% levels(x[[]][, cell_label])) {
        legend <- levels(x[[]][, cell_label])
    } else {
        legend <- NULL
    }

    table(Idents(x), x[[]][, cell_label]) %>%
        as.data.frame.matrix() %>%
        t() %>%
        as.data.frame() %>%
        mutate(across(everything(), ~ .x / sum(.x) * 100)) %>%
        plot_bar_2cat(
            count = TRUE,
            stats = FALSE,
            colour =  cols,
            pct  = FALSE,
            threshold = 1,
            digits = 1,
            legend = legend
        ) +
        labs(y = "% Cells")
}

table_ident_group <- function(
        x,
        group.by = "Patient",
        idents = "RNA__snn_res.0.35"
    ) {
    list.map(
        unique(x[[]][, group.by]),
        f(i) ~ filter(x[[]], !!sym(group.by) == i) %>%
            pull(idents) %>%
            fct_drop() %>%
            table()
    ) %>%
        list.rbind() %>%
        as.data.frame()
}


bar_ident_group0 <- function(
        x,
        group.by = "Patient",
        idents = "RNA__snn_res.0.35",
        ...
    ) {
    x %>%
        mutate(across(everything(), ~ .x / sum(.x) * 100)) %>%
        plot_bar_2cat(
            count = TRUE,
            stats = FALSE,
            pct  = FALSE,
            digits = 1,
            ...
        ) +
        labs(y = "% Cells")
}

#' @export
bar_ident_group <- function(
        x,
        group.by = "Patient",
        idents = "RNA__snn_res.0.35",
        ...
    ) {
    table_ident_group(x, group.by, idents) %>%
        bar_ident_group0(...) +
        geom_hline(yintercept = 50, color = "gray")
}

#' @export
bar_group_ident <- function(x, group.by = "Patient", idents = "RNA__snn_res.0.35", colour = palette_discrete()) {

    if ("0" %in% levels(x[[]][, idents])) {
        legend <- levels(x[[]][, idents])
    } else {
        legend <- reorder_celltype(x[[]][, idents])
    }
    table_ident_group(x, group.by, idents) %>%
        t() %>%
        as.data.frame() %>%
        bar_ident_group0(
            colour = colour,
            legend = legend
            )
}

#' @export
confusing_table <- function(
        x,
        new_ident = "predicted.labels",
        old_ident = "RNA_snn_res.0.35",
        as_percent = TRUE,
        margin = 1,
        cluster_rows = FALSE,
        cluster_cols = FALSE
) {
    tab_raw <- table(
        x[[]][, old_ident],
        x[[]][, new_ident]
    )

    # tab_raw <- tab_raw[
    #     as.character(0:(nrow(tab_raw) - 1)),
    #     as.character(0:(ncol(tab_raw) - 1))
    # ]

    if (as_percent) {
        breaks <- seq(0, 100, by = 10)
        display_matrix <- prop.table(tab_raw, margin = margin) * 100
        tab <- round(display_matrix, 1)
        legend_labels <- breaks
    } else {
        tab <- tab_raw
        tab_log <- log10(tab_raw + 1)
        breaks <- seq(min(tab_log), max(tab_log), length.out = 12)
        display_matrix <- tab_log
        legend_labels <- 10^breaks - 1
        legend_labels <- round(legend_labels)
        legend_labels[legend_labels < 0] <- 0
    }


    legend_breaks <- breaks
    palette_custom <- colorRampPalette(c("white", brewer.pal(9, "Reds")))(length(breaks) - 1)

    pheatmap(
        display_matrix,
        color = palette_custom,
        breaks = breaks,
        legend_breaks = legend_breaks,
        legend_labels = legend_labels,
        display_numbers = tab,
        na_col = "white",
        number_color = "white",
        fontsize_number = 10,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols
    )
}

#' @export
heatmap_celltype <- function(
        x,
        new_ident = "predicted.labels",
        old_ident = "cell_subtype_formatted",
        is_cluster = FALSE,
        ...
    ) {
    confusing_table(
        x,
        new_ident = new_ident,
        old_ident = old_ident,
        cluster_rows = is_cluster,
        cluster_cols = is_cluster,
        ...
    )
}

#' @export
plot_dim_cell <- function(x, cols = palette_discrete(), ...) {
    Idents(x) <- factor(Idents(x), levels =  reorder_celltype(Idents(x)))
    plot_dim(x, ...) +
        scale_color_manual(na.translate = FALSE, values = cols)
}
