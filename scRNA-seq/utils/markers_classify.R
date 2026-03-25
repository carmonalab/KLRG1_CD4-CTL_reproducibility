# markers_classify.R
#
# Generic helpers for annotating any gene-symbol set as cell-surface receptors
# using GO molecular-function and cellular-component terms.

# classify_gene_symbols_surface_receptors --------------------------------------
#
# Purpose:
#   Classify a vector of gene symbols into receptor-related sets.
#
# Inputs:
#   gene_symbols    : character vector of HGNC gene symbols.
#   receptor_terms  : regex used to detect receptor molecular function terms.
#   membrane_terms  : regex used to detect membrane/cell-surface localization.
#
# Output:
#   A list containing mapping tables and three key gene sets:
#   - receptor_genes
#   - membrane_genes
#   - membrane_receptor_genes (intersection of the two)

classify_gene_symbols_surface_receptors <- function(
    gene_symbols,
    receptor_terms = "receptor activity",
    membrane_terms = "plasma membrane|cell surface|external side of plasma membrane"
) {
  stopifnot(is.character(gene_symbols))

  gene_symbols <- unique(gene_symbols)
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & nzchar(gene_symbols)]

  if (length(gene_symbols) == 0) {
    stop("gene_symbols is empty after removing missing/blank values.")
  }

  # 1) Map gene symbols to Entrez IDs.
  symbol2entrez <- AnnotationDbi::mapIds(
    x = org.Hs.eg.db::org.Hs.eg.db,
    keys = gene_symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )

  gene_map <- tibble::tibble(
    gene = names(symbol2entrez),
    ENTREZID = as.character(symbol2entrez)
  ) %>%
    dplyr::filter(!is.na(ENTREZID))

  # 2) Retrieve GO terms in MF and CC ontologies.
  go_annotations <- AnnotationDbi::select(
    x = org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(gene_map$ENTREZID),
    columns = c("GO", "ONTOLOGY"),
    keytype = "ENTREZID"
  ) %>%
    dplyr::filter(!is.na(GO), ONTOLOGY %in% c("CC", "MF")) %>%
    dplyr::mutate(go_term = as.character(AnnotationDbi::Term(GO.db::GOTERM[GO])))

  # 3) Detect genes with receptor activity (MF).
  receptor_genes <- gene_map %>%
    dplyr::inner_join(go_annotations %>% dplyr::filter(ONTOLOGY == "MF"), by = "ENTREZID") %>%
    dplyr::filter(grepl(receptor_terms, go_term, ignore.case = TRUE)) %>%
    dplyr::pull(gene) %>%
    unique()

  # 4) Detect genes localized at the membrane/cell surface (CC).
  membrane_genes <- gene_map %>%
    dplyr::inner_join(go_annotations %>% dplyr::filter(ONTOLOGY == "CC"), by = "ENTREZID") %>%
    dplyr::filter(grepl(membrane_terms, go_term, ignore.case = TRUE)) %>%
    dplyr::pull(gene) %>%
    unique()

  # 5) Final cell-surface receptor set.
  membrane_receptor_genes <- intersect(receptor_genes, membrane_genes)

  list(
    gene_symbols = gene_symbols,
    gene_map = gene_map,
    go_annotations = go_annotations,
    receptor_genes = receptor_genes,
    membrane_genes = membrane_genes,
    membrane_receptor_genes = membrane_receptor_genes
  )
}

# annotate_marker_table_with_receptor_class ------------------------------------
#
# Purpose:
#   Add a binary receptor annotation column to any marker table.
#
# Inputs:
#   marker_table             : data frame with at least one gene-symbol column.
#   membrane_receptor_genes  : character vector from
#                              classify_gene_symbols_surface_receptors().
#   gene_col                 : marker_table column containing gene symbols.
#   class_col                : output column name to add.
#   positive_label           : label for receptor genes.
#   negative_label           : label for non-receptor genes.
#
# Output:
#   marker_table with an added class_col.

annotate_marker_table_with_receptor_class <- function(
    marker_table,
    membrane_receptor_genes,
    gene_col = "gene",
    class_col = "membrane_receptor",
    positive_label = "yes",
    negative_label = "no"
) {
  stopifnot(is.data.frame(marker_table))
  stopifnot(gene_col %in% colnames(marker_table))

  marker_table %>%
    dplyr::mutate(
      !!class_col := ifelse(
        .data[[gene_col]] %in% membrane_receptor_genes,
        positive_label,
        negative_label
      )
    )
}
