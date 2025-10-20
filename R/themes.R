#' @export
pal_sc <- c(
    rev(brewer.pal(9, "Blues")[-c(1, 3, 5, 7, 9)]),
    "white",
    # "grey90",
    brewer.pal(9, "Reds")[-c(1, 3, 5, 7, 9)]
)

#' @export
kable_sc <- function(x, func = str_clean, digits = getOption("digits")) {
    x %>%
        set_colnames(colnames(.) %>% func()) %>%
        kable0(digits = digits)
}

#' @export
str_clean <- function(x) {
    str_replace_all(x, "[_\\.]", " ") %>%
        str_replace_all("^n(Count|Feature)", "# \\1") %>%
        str_replace_all("# Feature RNA", "# Gene RNA") %>%
        str_replace_all("# Count ATAC", "# Fragment ATAC") %>%
        str_replace_all("# Feature ATAC", "# Peak ATAC") %>%
        str_remove_all(assay) %>%
        str_replace_all("percent", "% ") %>%
        str_to_sentence() %>%
        str_replace_all("# Count", "# UMI") %>%
        str_replace_all(" Dup$"," Duplicated") %>%
        str_trim()
}

#' @export
theme_sc_dim <- function(p, dims = seq(2)) {
    p +
        theme_custom() +
        labs(
            x = paste("Dim.", dims[1]),
            y = paste("Dim.", dims[2])
        ) +
        geom_hline(yintercept = 0, linewidth = 1, color = "gray30") +
        geom_vline(xintercept = 0, linewidth = 1, color = "gray30")
}

#' @export
format_annot <- function(x) {
    str_replace_all(x, "(cell)s", "\\1") %>%
        str_replace_all("_", " ") %>%
        str_replace_all("\\:", ", ") %>%
        str_replace_all("([Mm]onocyte)s", "\\1") %>%
        str_replace_all("([Nn]eutrophil)s", "\\1") %>%
        str_replace_all("([Bb]asophil)s", "\\1") %>%
        str_replace_all("([Pp]rogenitor)s", "\\1")  %>%
        str_replace_all("\\.Fr", "\\.FR") %>%
        str_replace_all("(B cell), pro", "\\1")
}

remove_parenthesis <- function(x) {
    str_remove(x, "\\s*\\(.*\\)")
}

keep_parenthesis <- function(x) {
    str_extract(x, "(?<=\\()[^\\)]+")
}
