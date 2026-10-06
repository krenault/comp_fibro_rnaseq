strip <- function(x) sub("\\.[0-9]+$", "", as.character(x))

geo_suppl_url <- function(gse, filename) {
  n <- as.integer(sub("^GSE", "", gse))
  prefix <- sprintf("GSE%dnnn", n %/% 1000L)
  sprintf(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s",
    prefix, gse, filename
  )
}

download_geo_suppl <- function(gse, filename, dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(dest_dir, filename)
  if (file.exists(dest) && file.info(dest)$size > 100) return(dest)
  url <- geo_suppl_url(gse, filename)
  message("Downloading ", url)
  utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  dest
}

map_mouse_ids_to_ensg <- function(ids) {
  ids <- unique(strip(ids))
  ids <- ids[nzchar(ids) & !is.na(ids)]
  if (!length(ids)) return(tibble::tibble(from_id = character(), gene = character()))
  is_ens <- grepl("^ENSMUSG", ids)
  out <- list()
  chunks <- function(x, n = 1500L) split(x, ceiling(seq_along(x) / n))
  map_chunk <- function(s, id_type) {
    tryCatch({
      o <- babelgene::orthologs(genes = s, species = "mouse", human = FALSE)
      nm <- names(o)
      from_col <- if (id_type == "ensmusg" && "ensembl" %in% nm) "ensembl" else "symbol"
      tibble::tibble(
        from_id = strip(as.character(o[[from_col]])),
        gene = strip(as.character(o$human_ensembl))
      ) %>%
        dplyr::filter(nzchar(from_id), nzchar(gene), !is.na(gene)) %>%
        dplyr::distinct(from_id, .keep_all = TRUE)
    }, error = function(e) {
      message("map fail (", id_type, "): ", conditionMessage(e))
      tibble::tibble(from_id = character(), gene = character())
    })
  }
  if (any(is_ens)) {
    out <- c(out, lapply(chunks(ids[is_ens]), map_chunk, id_type = "ensmusg"))
  }
  if (any(!is_ens)) {
    out <- c(out, lapply(chunks(ids[!is_ens]), map_chunk, id_type = "symbol"))
  }
  dplyr::bind_rows(out) %>% dplyr::distinct(from_id, .keep_all = TRUE)
}

collapse_to_human <- function(df, map, gene_raw_col = "gene_raw") {
  df %>%
    dplyr::mutate(gene_raw = strip(.data[[gene_raw_col]])) %>%
    dplyr::inner_join(map, by = c("gene_raw" = "from_id")) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      log2FoldChange = mean(log2FoldChange, na.rm = TRUE),
      padj = suppressWarnings(min(padj, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(padj = ifelse(is.infinite(padj), NA_real_, padj)) %>%
    dplyr::filter(is.finite(log2FoldChange), grepl("^ENSG", gene))
}

run_deseq_pair <- function(mat, group, treat, ctrl = "control") {
  keep <- group %in% c(treat, ctrl)
  m <- mat[, keep, drop = FALSE]
  cd <- data.frame(
    row.names = colnames(m),
    group = factor(group[keep], levels = c(ctrl, treat))
  )
  keep_g <- rowSums(m >= 5) >= 2L
  dds <- DESeq2::DESeqDataSetFromMatrix(m[keep_g, , drop = FALSE], cd, ~ group)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  as.data.frame(DESeq2::results(dds, contrast = c("group", treat, ctrl))) %>%
    tibble::rownames_to_column("gene_raw") %>%
    dplyr::filter(is.finite(log2FoldChange)) %>%
    dplyr::transmute(
      gene_raw = strip(gene_raw),
      log2FoldChange = as.numeric(log2FoldChange),
      padj = as.numeric(padj)
    )
}

lfc_from_abundance <- function(mat, treat_cols, ctrl_cols, pseudocount = 1) {
  tt <- as.matrix(mat[, treat_cols, drop = FALSE])
  cc <- as.matrix(mat[, ctrl_cols, drop = FALSE])
  tt[!is.finite(tt)] <- 0
  cc[!is.finite(cc)] <- 0
  tibble::tibble(
    gene_raw = strip(rownames(mat)),
    log2FoldChange = rowMeans(log2(tt + pseudocount)) - rowMeans(log2(cc + pseudocount)),
    padj = NA_real_
  ) %>% dplyr::filter(is.finite(log2FoldChange), nzchar(gene_raw))
}

write_geo_deseq <- function(hum, path, dataset, stress_panel, comparison, method) {
  out <- hum %>%
    dplyr::mutate(
      dataset = dataset,
      stress_panel = stress_panel,
      comparison = comparison,
      method = method
    )
  readr::write_csv(out, path)
  message("Wrote ", basename(path), "  genes=", nrow(out))
  invisible(out)
}

ensg_to_symbol <- function(genes) {
  genes <- unique(strip(genes))
  tryCatch({
    AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db, keys = genes,
      columns = "SYMBOL", keytype = "ENSEMBL"
    ) %>%
      tibble::as_tibble() %>%
      dplyr::transmute(gene = strip(ENSEMBL), symbol = SYMBOL) %>%
      dplyr::filter(!is.na(symbol)) %>%
      dplyr::distinct(gene, .keep_all = TRUE)
  }, error = function(e) tibble::tibble(gene = character(), symbol = character()))
}

safe_spearman <- function(x, y, min_n = 80L) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < min_n) return(list(n = n, rho = NA_real_, p = NA_real_))
  if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(list(n = n, rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(n = n, rho = unname(ct$estimate), p = ct$p.value)
}

frac_agree_fun <- function(z) {
  s <- sign(stats::median(z, na.rm = TRUE))
  if (!is.finite(s) || s == 0) mean(abs(z) < 1e-8, na.rm = TRUE) else mean(sign(z) == s, na.rm = TRUE)
}
