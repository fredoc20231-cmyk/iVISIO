# app.R
# Complete 10x Visium spatial-omics analysis application in R Shiny.
#
# Supported input:
#   * A Space Ranger output directory containing spatial/, filtered_feature_bc_matrix.h5
#     (or filtered_feature_bc_matrix/), and tissue_positions_list.csv or tissue_positions.csv.
#   * Optional metadata CSV with a barcode column.
#
# Main workflow:
#   Load -> QC -> normalization -> dimensional reduction -> clustering -> spatial maps
#   -> gene/feature exploration -> marker genes -> differential expression
#   -> spatially variable features -> spot-level signature scoring -> export/report.
#
# Install packages outside this app. Do not install packages at application startup.
# Bioconductor packages can be installed with BiocManager::install().
# Remove Shiny's built-in request-size cap. Reverse proxies, hosting platforms,
# browser infrastructure, and cloud load balancers may still impose limits.
# Use a finite value only when your deployment requires one.
options(shiny.maxRequestSize = Inf)

required_cran <- c("shiny", "bslib", "shinycssloaders", "DT", "plotly", "ggplot2",
                   "dplyr", "readr", "tidyr", "zip", "rmarkdown", "png", "jsonlite")
required_bioc <- c("Seurat", "SeuratObject", "sp", "sctransform", "SummarizedExperiment")

missing_packages <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
missing <- missing_packages(c(required_cran, required_bioc))
if (length(missing)) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "),
    "\nInstall dependencies before starting the app; package installation is intentionally not performed by app.R.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(shinycssloaders)
  library(DT)
  library(plotly)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(Seurat)
  library(SeuratObject)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

safe_filename <- function(x) {
  x <- trimws(as.character(x %||% "ivisio_analysis"))[1]
  if (is.na(x) || !nzchar(x)) x <- "ivisio_analysis"
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  if (!nzchar(x)) x <- "ivisio_analysis"
  substr(x, 1, 100)
}

log_message <- function(rv, text) {
  rv$log <- c(rv$log, sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text))
}

require_object <- function(rv) {
  if (is.null(rv$obj)) stop("Load a Visium dataset first.", call. = FALSE)
  rv$obj
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[[1]] else NULL
}

find_visium_matrix <- function(root) {
  first_existing(c(
    file.path(root, "filtered_feature_bc_matrix.h5"),
    file.path(root, "outs", "filtered_feature_bc_matrix.h5"),
    file.path(root, "filtered_feature_bc_matrix", "matrix.mtx.gz"),
    file.path(root, "outs", "filtered_feature_bc_matrix", "matrix.mtx.gz")
  ))
}

find_spatial_dir <- function(root) {
  candidates <- c(file.path(root, "spatial"), file.path(root, "outs", "spatial"))
  hit <- candidates[dir.exists(candidates)]
  if (length(hit)) hit[[1]] else NULL
}

build_spatial_overlay <- function(root) {
  if (!requireNamespace("png", quietly = TRUE)) stop("The png package is required for direct spatial overlays.", call. = FALSE)
  spatial_dir <- find_spatial_dir(root)
  image_path <- first_existing(c(file.path(spatial_dir, "tissue_hires_image.png"), file.path(spatial_dir, "tissue_lowres_image.png")))
  scale_path <- first_existing(c(file.path(spatial_dir, "scalefactors_json.json")))
  position_path <- first_existing(c(file.path(spatial_dir, "tissue_positions_list.csv"), file.path(spatial_dir, "tissue_positions.csv")))
  if (is.null(image_path) || is.null(scale_path) || is.null(position_path)) stop("The image, scale-factor JSON, and tissue-position CSV are all required for the spatial overlay.", call. = FALSE)
  img <- png::readPNG(image_path)
  img_info <- dim(img)
  height <- img_info[1]; width <- img_info[2]
  scale <- jsonlite::fromJSON(scale_path)
  scale_factor <- if (grepl("hires", basename(image_path), ignore.case = TRUE)) scale$tissue_hires_scalef else scale$tissue_lowres_scalef
  pos <- read.csv(position_path, header = FALSE, stringsAsFactors = FALSE)
  if (tolower(as.character(pos[1, 1])) == "barcode") pos <- read.csv(position_path, header = TRUE, stringsAsFactors = FALSE)
  if (ncol(pos) < 6) stop("Tissue-position CSV must contain six columns.", call. = FALSE)
  pos <- pos[, 1:6, drop = FALSE]
  names(pos) <- c("barcode", "tissue", "row", "col", "imagerow", "imagecol")
  pos$barcode <- as.character(pos$barcode)
  pos$tissue <- as.integer(pos$tissue)
  pos$imagerow <- as.numeric(pos$imagerow) * scale_factor
  pos$imagecol <- as.numeric(pos$imagecol) * scale_factor
  pos <- pos[pos$tissue == 1 & is.finite(pos$imagerow) & is.finite(pos$imagecol), , drop = FALSE]
  pos$x <- pos$imagecol
  pos$y <- height - pos$imagerow
  list(path = image_path, raster = as.raster(img), width = width, height = height, scale_factor = scale_factor, positions = pos)
}

normalize_visium_barcode <- function(x) sub("-[0-9]+$", "", trimws(as.character(x)))

# Discrete palette generator: uses a colorblind-friendly, high-contrast set for
# small numbers of groups (clusters/cell types) and falls back to an
# interpolated hue wheel for larger numbers so labels stay distinguishable.
ivisio_discrete_palette <- function(n) {
  base_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFD92F",
                "#A65628", "#F781BF", "#66C2A5", "#8DA0CB", "#E78AC3", "#A6D854",
                "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D")
  if (n <= length(base_pal)) return(base_pal[seq_len(n)])
  grDevices::colorRampPalette(base_pal)(n)
}

spatial_label_overlay_plot <- function(obj, overlay, group_col, point_size = 2.2) {
  if (is.null(group_col) || !group_col %in% colnames(obj@meta.data)) stop("Choose a metadata column containing clusters or labels.", call. = FALSE)
  pos <- overlay$positions
  pos$barcode <- as.character(pos$barcode)
  meta_barcodes <- rownames(obj@meta.data)
  matched <- !is.na(pos$barcode) & pos$barcode %in% meta_barcodes
  if (!any(matched)) {
    # Some pipelines append/remove a suffix such as '-1'; try a normalized match once,
    # the same fallback already used for expression overlays.
    meta_norm <- normalize_visium_barcode(meta_barcodes)
    pos_norm <- normalize_visium_barcode(pos$barcode)
    matched <- !is.na(pos_norm) & pos_norm %in% meta_norm
    if (!any(matched)) stop("No tissue-position barcodes matched the Seurat object. Check that the tissue-position CSV and loaded object come from the same Visium library.", call. = FALSE)
    meta_index <- match(pos_norm[matched], meta_norm)
  } else {
    meta_index <- match(pos$barcode[matched], meta_barcodes)
  }
  pos <- pos[matched, , drop = FALSE]
  meta <- as.character(obj@meta.data[meta_index, group_col])
  keep <- !is.na(meta) & nzchar(meta)
  if (!any(keep)) stop("No tissue-position barcodes have labels in the selected metadata column.", call. = FALSE)
  df <- pos[keep, , drop = FALSE]
  df$label <- factor(meta[keep])
  df$hover <- sprintf("Spot: %s<br>%s: %s", df$barcode, group_col, df$label)
  n_groups <- nlevels(df$label)
  ggplot(df, aes(x = x, y = y, color = label, text = hover)) +
    annotation_raster(overlay$raster, xmin = 0, xmax = overlay$width, ymin = 0, ymax = overlay$height, interpolate = TRUE) +
    geom_point(size = point_size, alpha = 0.95, shape = 16) +
    scale_color_manual(values = ivisio_discrete_palette(n_groups), name = group_col) +
    coord_fixed(xlim = c(0, overlay$width), ylim = c(0, overlay$height), expand = FALSE) +
    labs(title = sprintf("iVisio spatial labels: %s", group_col), subtitle = sprintf("Labeled tissue spots: %d | Groups: %d", nrow(df), n_groups), x = NULL, y = NULL) +
    theme_void(base_size = 12) + theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

cell_overlay_plot <- function(obj, overlay, cells, point_size = 1.4) {
  if (!all(c("pixel_x", "pixel_y") %in% names(cells))) stop("Cell overlay requires pixel_x and pixel_y columns in the cells-on-spot CSV.", call. = FALSE)
  df <- cells
  df$x <- as.numeric(df$pixel_x) * overlay$scale_factor
  df$y <- overlay$height - as.numeric(df$pixel_y) * overlay$scale_factor
  df <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
  if (!nrow(df)) stop("The cells-on-spot file contains no finite pixel coordinates.", call. = FALSE)
  color_col <- if ("predicted_cell_type" %in% names(df)) "predicted_cell_type" else if ("cell_type" %in% names(df)) "cell_type" else NULL
  cell_id_col <- if ("cell_id" %in% names(df)) "cell_id" else NULL
  df$hover <- if (!is.null(color_col)) {
    sprintf("Cell: %s<br>Type: %s", if (!is.null(cell_id_col)) df[[cell_id_col]] else seq_len(nrow(df)), df[[color_col]])
  } else {
    sprintf("Cell: %s", if (!is.null(cell_id_col)) df[[cell_id_col]] else seq_len(nrow(df)))
  }
  if (!is.null(color_col)) df[[color_col]] <- factor(df[[color_col]])
  aes_args <- if (is.null(color_col)) aes(x = x, y = y, text = hover) else aes(x = x, y = y, color = .data[[color_col]], text = hover)
  p <- ggplot(df, aes_args) +
    annotation_raster(overlay$raster, xmin = 0, xmax = overlay$width, ymin = 0, ymax = overlay$height, interpolate = TRUE) +
    geom_point(size = point_size, alpha = 0.9, shape = 16) +
    coord_fixed(xlim = c(0, overlay$width), ylim = c(0, overlay$height), expand = FALSE) +
    labs(title = "iVisio cell overlay", subtitle = sprintf("Cells plotted: %d%s", nrow(df), if (!is.null(color_col)) sprintf(" | Cell types: %d", nlevels(df[[color_col]])) else ""), x = NULL, y = NULL, color = "Cell type") +
    theme_void(base_size = 12) + theme(plot.title = element_text(face = "bold"), legend.position = "right")
  if (!is.null(color_col)) p <- p + scale_color_manual(values = ivisio_discrete_palette(nlevels(df[[color_col]])))
  p
}

# Summarizes the predicted/observed cell-type composition of a cells-on-spot
# data frame as counts and proportions, for bar/donut style composition charts.
cell_type_composition <- function(cells) {
  color_col <- if ("predicted_cell_type" %in% names(cells)) "predicted_cell_type" else if ("cell_type" %in% names(cells)) "cell_type" else NULL
  if (is.null(color_col)) return(NULL)
  tab <- table(as.character(cells[[color_col]]))
  data.frame(cell_type = names(tab), n_cells = as.integer(tab), fraction = as.numeric(tab) / sum(tab), stringsAsFactors = FALSE)
}

# Generic continuous-value overlay: plots any named numeric vector (keyed by
# Visium barcode) on the tissue image. This powers gene expression, STIE spot
# probabilities, and reference-deconvolution proportions with one consistent,
# interactive visualization.
spatial_values_overlay_plot <- function(overlay, barcodes, values, legend_name, title, point_size = 2.2, min_cutoff = "q01", max_cutoff = "q99") {
  if (is.null(overlay)) stop("Spatial overlay coordinates are unavailable. Reload the Visium files.", call. = FALSE)
  pos <- overlay$positions
  pos$barcode <- as.character(pos$barcode)
  barcodes <- as.character(barcodes)
  matched <- !is.na(pos$barcode) & pos$barcode %in% barcodes
  if (!any(matched)) {
    b_norm <- normalize_visium_barcode(barcodes)
    pos_norm <- normalize_visium_barcode(pos$barcode)
    matched <- !is.na(pos_norm) & pos_norm %in% b_norm
    if (!any(matched)) stop("No tissue-position barcodes matched the provided values. Check that the data come from the same Visium library.", call. = FALSE)
    idx <- match(pos_norm[matched], b_norm)
  } else {
    idx <- match(pos$barcode[matched], barcodes)
  }
  pos <- pos[matched, , drop = FALSE]
  df <- pos
  df$value <- as.numeric(values[idx])
  df <- df[is.finite(df$value), , drop = FALSE]
  if (!nrow(df)) stop("Matched tissue spots contain no finite values to plot.", call. = FALSE)
  parse_cut <- function(x, values, fallback) {
    x <- trimws(x %||% "")
    if (!nzchar(x)) return(fallback(values))
    if (grepl("^q[0-9]{2}$", x, ignore.case = TRUE)) return(as.numeric(quantile(values, as.numeric(substr(x, 2, 3)) / 100, na.rm = TRUE, names = FALSE)))
    val <- suppressWarnings(as.numeric(x)); if (is.finite(val)) val else fallback(values)
  }
  lo <- parse_cut(min_cutoff, df$value, function(z) quantile(z, 0.01, na.rm = TRUE, names = FALSE))
  hi <- parse_cut(max_cutoff, df$value, function(z) quantile(z, 0.99, na.rm = TRUE, names = FALSE))
  if (!is.finite(lo)) lo <- min(df$value, na.rm = TRUE)
  if (!is.finite(hi)) hi <- max(df$value, na.rm = TRUE)
  if (hi <= lo) hi <- lo + 1e-9
  df$value_plot <- pmin(pmax(df$value, lo), hi)
  df$hover <- sprintf("Spot: %s<br>%s: %.3f", df$barcode, legend_name, df$value)
  ggplot(df, aes(x = x, y = y, color = value_plot, text = hover)) +
    annotation_raster(overlay$raster, xmin = 0, xmax = overlay$width, ymin = 0, ymax = overlay$height, interpolate = TRUE) +
    geom_point(size = point_size, alpha = 0.9, shape = 16) +
    scale_color_viridis_c(option = "magma", name = legend_name, limits = c(lo, hi)) +
    coord_fixed(xlim = c(0, overlay$width), ylim = c(0, overlay$height), expand = FALSE) +
    labs(title = title, subtitle = sprintf("Dots: %d tissue spots", nrow(df)), x = NULL, y = NULL) +
    theme_void(base_size = 12) + theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

spatial_overlay_plot <- function(obj, overlay, feature, assay, point_size = 2.2, min_cutoff = "q01", max_cutoff = "q99") {
  if (is.null(overlay)) stop("Spatial overlay coordinates are unavailable. Reload the Visium files.", call. = FALSE)
  if (!feature %in% available_features(obj, assay)) stop(sprintf("Feature '%s' is not present in assay '%s'.", feature, assay), call. = FALSE)
  vals <- tryCatch(GetAssayData(obj, assay = assay, layer = "data"), error = function(e) GetAssayData(obj, assay = assay, slot = "data"))
  expr_barcodes <- colnames(vals)
  pos <- overlay$positions
  pos$barcode <- as.character(pos$barcode)
  matched <- !is.na(pos$barcode) & pos$barcode %in% expr_barcodes
  if (!any(matched)) {
    # Some pipelines append/remove a suffix such as '-1'; try a normalized match once.
    expr_norm <- normalize_visium_barcode(expr_barcodes)
    pos_norm <- normalize_visium_barcode(pos$barcode)
    matched <- !is.na(pos_norm) & pos_norm %in% expr_norm
    if (!any(matched)) stop("No tissue-position barcodes matched the expression matrix. Check that the H5 matrix and tissue-position CSV come from the same Visium library.", call. = FALSE)
    expr_index <- match(pos_norm[matched], expr_norm)
  } else {
    expr_index <- match(pos$barcode[matched], expr_barcodes)
  }
  pos <- pos[matched, , drop = FALSE]
  # Never pass NA indices to a sparse matrix; all indices are valid after filtering.
  values <- as.numeric(vals[feature, expr_index, drop = TRUE])
  df <- pos
  df$value <- values
  df <- df[is.finite(df$value), , drop = FALSE]
  if (!nrow(df)) stop("Matched tissue spots contain no finite expression values for this feature.", call. = FALSE)
  parse_cut <- function(x, values, fallback) {
    x <- trimws(x %||% "")
    if (!nzchar(x)) return(fallback(values))
    if (grepl("^q[0-9]{2}$", x, ignore.case = TRUE)) return(as.numeric(quantile(values, as.numeric(substr(x, 2, 3)) / 100, na.rm = TRUE, names = FALSE)))
    val <- suppressWarnings(as.numeric(x)); if (is.finite(val)) val else fallback(values)
  }
  lo <- parse_cut(min_cutoff, df$value, function(z) quantile(z, 0.01, na.rm = TRUE, names = FALSE))
  hi <- parse_cut(max_cutoff, df$value, function(z) quantile(z, 0.99, na.rm = TRUE, names = FALSE))
  if (!is.finite(lo)) lo <- min(df$value, na.rm = TRUE)
  if (!is.finite(hi)) hi <- max(df$value, na.rm = TRUE)
  if (hi <= lo) hi <- lo + 1e-9
  df$value_plot <- pmin(pmax(df$value, lo), hi)
  df$hover <- sprintf("Spot: %s<br>%s: %.3f", df$barcode, feature, df$value)
  ggplot(df, aes(x = x, y = y, color = value_plot, text = hover)) +
    annotation_raster(overlay$raster, xmin = 0, xmax = overlay$width, ymin = 0, ymax = overlay$height, interpolate = TRUE) +
    geom_point(size = point_size, alpha = 0.9, shape = 16) +
    scale_color_viridis_c(option = "magma", name = feature, limits = c(lo, hi)) +
    coord_fixed(xlim = c(0, overlay$width), ylim = c(0, overlay$height), expand = FALSE) +
    labs(title = sprintf("iVisio spatial expression: %s", feature), subtitle = sprintf("Dots: %d tissue spots | cutoff: %s to %s", nrow(df), min_cutoff, max_cutoff), x = NULL, y = NULL) +
    theme_void(base_size = 12) + theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

read_visium <- function(root, slice, assay) {
  matrix_path <- find_visium_matrix(root)
  spatial_dir <- find_spatial_dir(root)
  if (is.null(matrix_path)) stop("No filtered_feature_bc_matrix.h5 or filtered_feature_bc_matrix directory was found.", call. = FALSE)
  if (is.null(spatial_dir)) stop("No spatial/ directory was found. Upload the H5 matrix plus tissue image, scale factors, and tissue positions.", call. = FALSE)

  if (!grepl("\\.h5$", matrix_path, ignore.case = TRUE)) {
    return(Load10X_Spatial(data.dir = root, assay = assay, slice = slice, filter.matrix = TRUE, to.upper = FALSE))
  }

  # Manual loading is used here because public datasets often provide only
  # tissue_hires_image.png rather than the standard tissue_lowres_image.png.
  # It also accepts the user's exact GSM-prefixed filenames after staging.
  counts <- Read10X_h5(matrix_path, use.names = TRUE, unique.features = TRUE)
  if (is.list(counts)) {
    preferred <- intersect(c("Gene Expression", "GeneExpression", assay), names(counts))
    counts <- if (length(preferred)) counts[[preferred[1]]] else counts[[1]]
  }
  counts <- as(counts, "dgCMatrix")
  image_path <- first_existing(c(file.path(spatial_dir, "tissue_hires_image.png"), file.path(spatial_dir, "tissue_lowres_image.png")))
  if (is.null(image_path)) stop("A tissue_hires_image.png or tissue_lowres_image.png file is required.", call. = FALSE)
  image <- suppressWarnings(Read10X_Image(image.dir = spatial_dir, image.name = basename(image_path), filter.matrix = TRUE))
  obj <- CreateSeuratObject(counts = counts, assay = assay, project = slice)
  obj[[slice]] <- suppressWarnings(image)
  obj
}

read_user_metadata <- function(path) {
  md <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (!"barcode" %in% names(md)) stop("Metadata CSV must contain a barcode column.", call. = FALSE)
  md <- as.data.frame(md)
  rownames(md) <- make.unique(as.character(md$barcode))
  md
}

read_stie_matrix <- function(path, row_name = "gene") {
  x <- as.data.frame(readr::read_csv(path, show_col_types = FALSE, progress = FALSE), check.names = FALSE)
  if (!ncol(x) || nrow(x) < 1) stop("The STIE CSV is empty.", call. = FALSE)
  first <- names(x)[1]
  rownames(x) <- make.unique(as.character(x[[first]]))
  x[[first]] <- NULL
  out <- as.matrix(x)
  suppressWarnings(storage.mode(out) <- "numeric")
  if (anyNA(out)) stop("The STIE signature CSV must contain only numeric expression values after its first identifier column.", call. = FALSE)
  out
}

reference_catalog <- data.frame(
  resource = c("10x adult mouse brain FFPE", "10x adult mouse kidney FFPE", "10x human breast cancer FFPE", "Human DLPFC HumanPilot10x", "10x CytAssist mouse brain sections 1/2", "SpatialResearch breast cancer", "HDST breast cancer H&E", "Mouse hippocampus snRNA-seq SCP1", "Allen adult mouse cortex taxonomy", "Breast cancer scRNA-seq GSE176078", "Lieber spatialLIBD DLPFC markers"),
  species = c("mouse", "mouse", "human", "human", "mouse", "human", "human", "mouse", "mouse", "human", "human"),
  modality = c("Visium FFPE", "Visium FFPE", "Visium FFPE", "Visium", "Visium CytAssist", "legacy spatial", "H&E image", "snRNA-seq", "scRNA-seq", "scRNA-seq", "marker genes"),
  source = c("10x Genomics", "10x Genomics", "10x Genomics/figshare", "Globus jhpce#HumanPilot10x / spatialLIBD", "10x Genomics", "spatialresearch.org", "HDST publication files", "Broad SCP1", "Allen Institute", "NCBI GEO GSE176078", "LieberInstitute/spatialLIBD"),
  stringsAsFactors = FALSE
)

read_reference_expression <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    x <- readRDS(path)
    if (inherits(x, "Seurat")) {
      assay <- DefaultAssay(x)
      mat <- tryCatch(GetAssayData(x, assay = assay, layer = "data"), error = function(e) GetAssayData(x, assay = assay, slot = "data"))
      groups <- if ("cell_type" %in% colnames(x@meta.data)) x$cell_type else if ("celltype" %in% colnames(x@meta.data)) x$celltype else NULL
      return(list(matrix = as.matrix(mat), groups = groups))
    }
    if (inherits(x, "SingleCellExperiment")) {
      mat <- SummarizedExperiment::assay(x, "logcounts")
      groups <- if ("cell_type" %in% colnames(SummarizedExperiment::colData(x))) SummarizedExperiment::colData(x)$cell_type else NULL
      return(list(matrix = as.matrix(mat), groups = groups))
    }
    if (is.matrix(x) || is.data.frame(x)) return(list(matrix = as.matrix(x), groups = NULL))
    stop("Unsupported RDS reference object.", call. = FALSE)
  }
  if (ext == "h5") {
    x <- Read10X_h5(path, use.names = TRUE, unique.features = TRUE)
    if (is.list(x)) x <- x[[1]]
    return(list(matrix = as.matrix(x), groups = NULL))
  }
  list(matrix = read_stie_matrix(path), groups = NULL)
}

build_reference_signature <- function(reference, groups) {
  mat <- reference$matrix
  if (is.null(groups) || length(groups) != ncol(mat)) stop("Reference RDS must contain a cell_type or celltype metadata column with one label per cell.", call. = FALSE)
  groups <- as.character(groups)
  keep <- !is.na(groups) & nzchar(groups)
  mat <- mat[, keep, drop = FALSE]; groups <- groups[keep]
  lev <- sort(unique(groups))
  sig <- sapply(lev, function(g) rowMeans(mat[, groups == g, drop = FALSE]))
  if (is.null(dim(sig))) sig <- matrix(sig, ncol = 1, dimnames = list(rownames(mat), lev[1]))
  rownames(sig) <- rownames(mat); colnames(sig) <- lev
  sig
}

simple_nnls_deconvolution <- function(st_expr, signature, max_iter = 500, tol = 1e-5) {
  genes <- intersect(colnames(st_expr), rownames(signature))
  if (length(genes) < 10) stop("Fewer than 10 genes overlap the Visium matrix and reference signature.", call. = FALSE)
  y <- as.matrix(st_expr[, genes, drop = FALSE])
  a <- as.matrix(signature[genes, , drop = FALSE])
  a[a < 0] <- 0; y[y < 0] <- 0
  scale_a <- sqrt(colSums(a^2)); scale_a[scale_a == 0] <- 1
  a <- sweep(a, 2, scale_a, "/")
  p <- matrix(1 / ncol(a), nrow = nrow(y), ncol = ncol(a), dimnames = list(rownames(y), colnames(a)))
  step <- 1 / max(1, max(colSums(crossprod(a))))
  for (iter in seq_len(max_iter)) {
    old <- p
    grad <- (p %*% crossprod(a) - y %*% a) * step
    p <- pmax(0, p - grad)
    rs <- rowSums(p); rs[rs == 0] <- 1; p <- p / rs
    if (max(abs(p - old)) < tol) break
  }
  p
}

add_metadata_safely <- function(obj, metadata) {
  common <- intersect(colnames(obj), rownames(metadata))
  if (!length(common)) stop("No metadata barcodes matched the Visium object.", call. = FALSE)
  add <- metadata[common, setdiff(names(metadata), "barcode"), drop = FALSE]
  obj <- AddMetaData(obj, metadata = add)
  obj
}

spatial_assay <- function(obj) {
  candidates <- c("Spatial", "SCT", "RNA")
  hit <- candidates[candidates %in% Assays(obj)]
  if (!length(hit)) Assays(obj)[1] else hit[1]
}

available_features <- function(obj, assay = NULL) {
  assay <- assay %||% DefaultAssay(obj)
  tryCatch(rownames(obj[[assay]]), error = function(e) character())
}

get_counts_matrix <- function(obj, assay) {
  # Seurat v5 uses layers; Seurat v4 uses slots. Try both so the app works
  # across commonly installed Seurat versions.
  out <- tryCatch(GetAssayData(obj, assay = assay, layer = "counts"), error = function(e) NULL)
  if (is.null(out)) out <- GetAssayData(obj, assay = assay, slot = "counts")
  out
}

spatial_plot <- function(obj, feature, image_name, point_size = 1.6, alpha = 0.9, min_cutoff = "q01", max_cutoff = "q99") {
  assay <- spatial_assay(obj)
  if (!assay %in% Assays(obj)) stop(sprintf("Expression assay '%s' is unavailable.", assay), call. = FALSE)
  DefaultAssay(obj) <- assay
  image_names <- Images(obj)
  if (!length(image_names)) stop("No Visium image is attached to the object.", call. = FALSE)
  image_name <- if (!is.null(image_name) && image_name %in% image_names) image_name else image_names[1]
  if (!feature %in% available_features(obj, assay)) stop(sprintf("Feature '%s' is not present in assay '%s'.", feature, assay), call. = FALSE)
  args <- list(object = obj, features = feature, images = image_name, pt.size.factor = point_size,
               alpha = c(alpha, alpha), crop = FALSE, combine = FALSE)
  if (!is.null(min_cutoff) && nzchar(trimws(min_cutoff))) args$min.cutoff <- trimws(min_cutoff)
  if (!is.null(max_cutoff) && nzchar(trimws(max_cutoff))) args$max.cutoff <- trimws(max_cutoff)
  plots <- do.call(SpatialFeaturePlot, args)
  p <- if (is.list(plots)) plots[[1]] else plots
  p + ggtitle(sprintf("Spatial expression: %s", feature)) + theme(plot.title = element_text(face = "bold"))
}

ui <- page_navbar(
  title = "iVisio | Spatial Omics Platform",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2166AC"),
  nav_panel("1. Load", icon = icon("folder-open"),
            layout_sidebar(
              sidebar = sidebar(
                h4("Space Ranger input"),
                fileInput("matrix_h5", "1. Filtered expression matrix (.h5)", accept = c(".h5")),
                fileInput("tissue_image", "2. Tissue image (.png)", accept = c(".png", ".jpg", ".jpeg")),
                fileInput("scale_factors", "3. Scale factors (.json)", accept = c(".json")),
                fileInput("tissue_positions", "4. Tissue positions (.csv)", accept = c(".csv", ".tsv")),
                helpText("Upload each Visium file separately. This avoids a large multi-file request and works with GSM-prefixed filenames."),
                textInput("slice", "Slice name", "slice1"),
                textInput("assay_name", "Assay name", "Spatial"),
                fileInput("metadata", "Optional metadata CSV", accept = ".csv"),
                textInput("project", "Project name", "Visium_Project"),
                actionButton("load", "Load Visium data", class = "btn-primary", icon = icon("play"))
              ),
              card(card_header("Dataset status"), verbatimTextOutput("status")),
              card(card_header("Run log"), tags$pre(class = "log", textOutput("log")))
            )
  ),
  nav_panel("2. QC", icon = icon("filter"),
            layout_sidebar(
              sidebar = sidebar(
                numericInput("min_counts", "Minimum counts", 500, min = 0),
                numericInput("max_counts", "Maximum counts", 100000, min = 1),
                numericInput("min_features", "Minimum detected features", 100, min = 0),
                numericInput("max_mt", "Maximum mitochondrial percentage", 25, min = 0, max = 100),
                textInput("mito_pattern", "Mitochondrial gene pattern", "^MT-", width = "100%"),
                actionButton("run_qc", "Calculate QC metrics", class = "btn-success"),
                actionButton("filter_qc", "Apply QC filters", class = "btn-warning")
              ),
              layout_columns(card(card_header("QC distributions"), plotOutput("qc_plot", height = 520) %>% withSpinner()),
                             card(card_header("QC scatter"), plotlyOutput("qc_scatter", height = 520) %>% withSpinner()))
            )
  ),
  nav_panel("3. Normalize", icon = icon("sliders"),
            layout_sidebar(
              sidebar = sidebar(
                selectInput("normalization", "Normalization", c("LogNormalize" = "log", "SCTransform" = "sct"), selected = "log"),
                numericInput("variable_features", "Variable features", 2000, min = 500, max = 10000),
                numericInput("n_pcs", "Principal components", 20, min = 5, max = 50),
                checkboxInput("fast_norm", "Fast mode", TRUE),
                helpText("Fast mode uses fewer features/PCs and approximate PCA where supported."),
                actionButton("run_norm", "Normalize and run PCA", class = "btn-success")
              ),
              layout_columns(card(card_header("PCA elbow"), plotOutput("elbow", height = 500) %>% withSpinner()),
                             card(card_header("PCA overview"), plotlyOutput("pca_plot", height = 500) %>% withSpinner()))
            )
  ),
  nav_panel("4. Clusters", icon = icon("project-diagram"),
            layout_sidebar(
              sidebar = sidebar(
                numericInput("dims_low", "First PC", 1, min = 1),
                numericInput("dims_high", "Last PC", 30, min = 2),
                numericInput("resolution", "Cluster resolution", 0.5, min = 0.05, max = 5, step = 0.05),
                numericInput("umap_neighbors", "UMAP neighbors", 30, min = 5, max = 100),
                selectInput("umap_metric", "UMAP metric", c("cosine", "euclidean", "manhattan")),
                checkboxInput("fast_clusters", "Fast clustering mode", TRUE),
                actionButton("run_clusters", "Neighbors, UMAP and clusters", class = "btn-success")
              ),
              layout_columns(
                card(card_header("UMAP"), plotlyOutput("umap_plot", height = 560) %>% withSpinner()),
                card(card_header("Clusters on tissue"), plotlyOutput("cluster_spatial_plot", height = 560) %>% withSpinner()),
                col_widths = c(6, 6)
              ),
              card(card_header("Cluster summary"), DTOutput("cluster_table"))
            )
  ),
  nav_panel("5. Spatial map", icon = icon("map"),
            layout_sidebar(
              sidebar = sidebar(
                selectInput("image_name", "Image", choices = NULL),
                selectInput("spatial_mode", "Overlay mode", choices = c("Expression feature" = "expression", "Metadata clusters/labels" = "labels"), selected = "expression"),
                selectizeInput("spatial_feature", "Feature", choices = NULL, options = list(placeholder = "Type a gene or feature"), multiple = FALSE),
                selectizeInput("spatial_group", "Metadata column for labels", choices = NULL, options = list(placeholder = "Select clusters or labels"), multiple = FALSE),
                numericInput("spatial_point_size", "Spot size factor", 1.6, min = 0.1, max = 5, step = 0.1),
                textInput("min_cutoff", "Minimum cutoff", "q01"),
                textInput("max_cutoff", "Maximum cutoff", "q99"),
                helpText("Defaults use the 1st and 99th percentiles. You may enter values such as 0, 2, q05, or q95."),
                actionButton("refresh_spatial", "Refresh spatial map", class = "btn-primary")
              ),
              card(card_header("Spatial feature map"), plotlyOutput("spatial_plot", height = 760) %>% withSpinner())
            )
  ),
  nav_panel("6. Explore", icon = icon("search"),
            layout_sidebar(
              sidebar = sidebar(
                selectizeInput("explore_feature", "Feature", choices = NULL, options = list(placeholder = "Type a gene or feature"), multiple = FALSE),
                selectInput("explore_reduction", "Reduction", choices = c("umap", "pca")),
                actionButton("explore_refresh", "Plot feature", class = "btn-primary"),
                hr(),
                textInput("signature_name", "Signature name", "My_signature"),
                textAreaInput("signature_genes", "Genes, one per line or comma-separated", "COL1A1\nCOL1A2\nDCN", rows = 5),
                actionButton("score_signature", "Score signature", class = "btn-success")
              ),
              layout_columns(card(card_header("Feature expression"), plotlyOutput("feature_plot", height = 650) %>% withSpinner()),
                             card(card_header("Signature scores"), DTOutput("signature_table")))
            )
  ),
  nav_panel("7. Markers / DE", icon = icon("microscope"),
            layout_sidebar(
              sidebar = sidebar(
                selectInput("marker_group", "Group column", choices = "seurat_clusters"),
                numericInput("marker_min_pct", "Minimum percentage", 0.1, min = 0, max = 1, step = 0.05),
                numericInput("marker_logfc", "Log fold-change threshold", 0.25, min = 0, step = 0.05),
                checkboxInput("fast_de", "Fast vectorized DE", TRUE),
                numericInput("de_max_genes", "Maximum genes to test", 2000, min = 100, max = 20000, step = 100),
                helpText("Fast DE uses vectorized group statistics and is recommended for Visium spot matrices."),
                actionButton("run_markers", "Find markers for all groups", class = "btn-success"),
                hr(),
                selectInput("de_group1", "Group 1", choices = NULL),
                selectInput("de_group2", "Group 2", choices = NULL),
                actionButton("run_de", "Run pairwise differential expression", class = "btn-warning")
              ),
              card(card_header("Differential results"), DTOutput("marker_table"))
            )
  ),
  nav_panel("8. STIE", icon = icon("microscope"),
            layout_sidebar(
              sidebar = sidebar(
                helpText("STIE integrates spot expression, nuclear morphology, and optional cell-type signatures. It requires the STIE R package plus a cells-on-spot CSV."),
                fileInput("stie_cells", "Cells-on-spot morphology CSV", accept = ".csv"),
                fileInput("stie_signature", "Optional gene × cell-type signature CSV", accept = ".csv"),
                selectInput("stie_mode", "STIE mode", c("Deconvolution with known signature" = "deconv", "Unsupervised cell-type clustering" = "cluster")),
                textInput("stie_features", "Morphology columns, comma-separated", "area, perimeter, circularity"),
                numericInput("stie_lambda", "Morphology shrinkage lambda", 0, min = 0, step = 0.1),
                numericInput("stie_steps", "STIE iterations", 10, min = 1, max = 100),
                numericInput("stie_min_cells", "Minimum cells per type", 2, min = 0),
                checkboxInput("stie_equal_prior", "Use equal cell-type prior", TRUE),
                actionButton("run_stie", "Run STIE", class = "btn-success"),
                actionButton("install_stie", "Install STIE dependencies", class = "btn-outline-primary"),
                verbatimTextOutput("stie_install_status"),
                helpText("Required cells-on-spot columns: cell_id, spot, plus the selected morphology columns. The spot values must match Visium barcodes."),
                helpText("The installer runs: remotes::install_github('zhushijia/STIE') and install.packages(c('quadprog','magick')). Restart R/Shiny after installation if requested."),
                hr(),
                selectInput("stie_prob_column", "Spot probability to display", choices = NULL)
              ),
              layout_columns(
                card(card_header("STIE spot probabilities"), plotlyOutput("stie_spatial", height = 560) %>% withSpinner()),
                card(card_header("Cells on tissue image"), plotlyOutput("stie_cells_spatial", height = 560) %>% withSpinner()),
                col_widths = c(6, 6)
              ),
              layout_columns(
                card(card_header("Predicted cell-type composition"), plotlyOutput("stie_composition_plot", height = 380) %>% withSpinner()),
                card(card_header("STIE results"), DTOutput("stie_table")),
                col_widths = c(5, 7)
              )
            )
  ),
  nav_panel("9. Reference / Deconvolution", icon = icon("database"),
            layout_sidebar(
              sidebar = sidebar(
                selectInput("reference_catalog_choice", "Reference resource", choices = setNames(seq_len(nrow(reference_catalog)), reference_catalog$resource)),
                fileInput("reference_file", "Reference RDS, H5, or expression CSV", accept = c(".rds", ".h5", ".csv")),
                fileInput("reference_signature", "Optional precomputed signature CSV", accept = ".csv"),
                actionButton("load_reference", "Load/build reference signature", class = "btn-primary"),
                actionButton("run_reference_deconv", "Run reference deconvolution", class = "btn-success"),
                helpText("Reference files are not silently downloaded. Download them from the cited source, then upload a local RDS/H5/CSV or a gene × cell-type signature CSV."),
                helpText("For an RDS, include cell_type or celltype metadata. For a signature CSV, use genes in the first column and cell types in the remaining columns."),
                hr(),
                selectInput("reference_cell_type", "Cell type to view on tissue", choices = NULL)
              ),
              layout_columns(card(card_header("Reference catalog"), DTOutput("reference_catalog_table")),
                             card(card_header("Signature preview"), DTOutput("reference_signature_table"))),
              layout_columns(
                card(card_header("Cell type proportion on tissue"), plotlyOutput("reference_spatial_plot", height = 480) %>% withSpinner()),
                card(card_header("Average composition across spots"), plotlyOutput("reference_composition_plot", height = 480) %>% withSpinner()),
                col_widths = c(6, 6)
              ),
              card(card_header("Deconvolution output"), DTOutput("reference_deconv_table"))
            )
  ),
  nav_panel("10. Spatially variable", icon = icon("bullseye"),
            layout_sidebar(
              sidebar = sidebar(
                numericInput("svg_n", "Number of features", 100, min = 10, max = 1000),
                actionButton("run_svg", "Find spatially variable features", class = "btn-success")
              ),
              layout_columns(card(card_header("Spatially variable features"), DTOutput("svg_table")),
                             card(card_header("Top spatial feature"), plotlyOutput("svg_plot", height = 520) %>% withSpinner()))
            )
  ),
  nav_panel("11. Export", icon = icon("download"),
            layout_sidebar(
              sidebar = sidebar(
                textInput("out_prefix", "Output prefix", "visium_analysis"),
                downloadButton("download_rds", "Download Seurat object"),
                downloadButton("download_markers", "Download marker table"),
                downloadButton("download_svg", "Download SVG table"),
                downloadButton("download_reference_deconv", "Download reference deconvolution"),
                downloadButton("download_reference_signature", "Download reference signature"),
                downloadButton("download_metadata", "Download spot metadata"),
                actionButton("make_report", "Create HTML report", class = "btn-primary"),
                downloadButton("download_report", "Download HTML report")
              ),
              card(card_header("Export status"), verbatimTextOutput("export_status"))
            )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(obj = NULL, log = character(), markers = NULL, de = NULL,
                       svg = NULL, signatures = list(), report = NULL, data_root = NULL,
                       stie = NULL, reference = NULL, reference_signature = NULL, reference_deconv = NULL,
                       spatial_overlay = NULL, stie_install = "Not installed or not checked.")

  output$log <- renderText(paste(rv$log, collapse = "\n"))
  output$status <- renderPrint({
    if (is.null(rv$obj)) return(cat("No object loaded.\n"))
    obj <- rv$obj
    cat("Project:", obj@project.name, "\n")
    cat("Spots:", ncol(obj), "\nFeatures:", nrow(obj[[spatial_assay(obj)]]), "\n")
    cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")
    cat("Reductions:", paste(Reductions(obj), collapse = ", "), "\n")
    cat("Images:", paste(Images(obj), collapse = ", "), "\n")
  })

  observeEvent(input$load, {
    req(input$matrix_h5, input$tissue_image, input$scale_factors, input$tissue_positions)
    tryCatch({
      log_message(rv, "Assembling the four uploaded Visium files...")
      root <- tempfile("visium_")
      dir.create(root, recursive = TRUE)
      dir.create(file.path(root, "spatial"), recursive = TRUE)
      file.copy(input$matrix_h5$datapath, file.path(root, "filtered_feature_bc_matrix.h5"), overwrite = TRUE)
      image_ext <- tolower(tools::file_ext(input$tissue_image$name))
      if (!image_ext %in% c("png", "jpg", "jpeg")) stop("The tissue image must be PNG, JPG, or JPEG.", call. = FALSE)
      image_target <- if (image_ext == "png") "tissue_hires_image.png" else paste0("tissue_hires_image.", image_ext)
      file.copy(input$tissue_image$datapath, file.path(root, "spatial", image_target), overwrite = TRUE)
      file.copy(input$scale_factors$datapath, file.path(root, "spatial", "scalefactors_json.json"), overwrite = TRUE)
      positions_target <- if (grepl("tissue_positions_list", input$tissue_positions$name, ignore.case = TRUE)) "tissue_positions_list.csv" else "tissue_positions.csv"
      file.copy(input$tissue_positions$datapath, file.path(root, "spatial", positions_target), overwrite = TRUE)
      obj <- read_visium(root, input$slice, input$assay_name)
      obj@project.name <- input$project
      if (!is.null(input$metadata)) obj <- add_metadata_safely(obj, read_user_metadata(input$metadata$datapath))
      rv$obj <- obj
      rv$data_root <- root
      rv$spatial_overlay <- build_spatial_overlay(root)
      log_message(rv, sprintf("Loaded %d spots and %d features.", ncol(obj), nrow(obj[[spatial_assay(obj)]])))
      updateSelectInput(session, "image_name", choices = Images(obj), selected = Images(obj)[1])
      feats <- available_features(obj, spatial_assay(obj))
      updateSelectizeInput(session, "spatial_feature", choices = feats, selected = feats[1], server = TRUE)
      updateSelectizeInput(session, "explore_feature", choices = feats, selected = feats[1], server = TRUE)
      metadata_choices <- unique(c("seurat_clusters", colnames(obj@meta.data)))
      updateSelectizeInput(session, "spatial_group", choices = metadata_choices, selected = if ("seurat_clusters" %in% metadata_choices) "seurat_clusters" else metadata_choices[1], server = TRUE)
      updateSelectInput(session, "marker_group", choices = metadata_choices)
    }, error = function(e) log_message(rv, paste("LOAD ERROR:", conditionMessage(e))))
  })

  observeEvent(input$run_qc, {
    obj <- require_object(rv)
    tryCatch({
      assay <- spatial_assay(obj); DefaultAssay(obj) <- assay
      counts <- get_counts_matrix(obj, assay)
      if (!"nCount_Spatial" %in% colnames(obj@meta.data)) obj$nCount_Spatial <- Matrix::colSums(counts)
      if (!"nFeature_Spatial" %in% colnames(obj@meta.data)) obj$nFeature_Spatial <- Matrix::colSums(counts > 0)
      mt <- grep(input$mito_pattern, rownames(obj[[assay]]), value = TRUE, ignore.case = TRUE)
      if (length(mt)) {
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = input$mito_pattern, assay = assay)
      } else {
        log_message(rv, sprintf("No genes matched the mitochondrial pattern '%s'; percent.mt was not computed.", input$mito_pattern))
      }
      rv$obj <- obj; log_message(rv, "QC metrics calculated.")
    }, error = function(e) log_message(rv, paste("QC ERROR:", conditionMessage(e))))
  })

  output$qc_plot <- renderPlot({
    obj <- require_object(rv); md <- obj@meta.data
    metrics <- intersect(c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), names(md))
    validate(need(length(metrics), "Run QC metrics first."))
    df <- tidyr::pivot_longer(cbind(spot = rownames(md), md[, metrics, drop = FALSE]), -spot)
    ggplot(df, aes(x = value)) + geom_histogram(bins = 50, fill = "#2C7FB8", color = "white") +
      facet_wrap(~name, scales = "free", ncol = 1) + theme_minimal(base_size = 12) + labs(x = NULL, y = "Spots")
  })

  output$qc_scatter <- renderPlotly({
    obj <- require_object(rv); md <- obj@meta.data
    validate(need(all(c("nCount_Spatial", "nFeature_Spatial") %in% names(md)), "Run QC metrics first."))
    plot_ly(md, x = ~nCount_Spatial, y = ~nFeature_Spatial, color = ~if ("percent.mt" %in% names(md)) percent.mt else nCount_Spatial,
            type = "scatter", mode = "markers", text = rownames(md), hoverinfo = "text+x+y")
  })

  observeEvent(input$filter_qc, {
    obj <- require_object(rv)
    tryCatch({
      md <- obj@meta.data
      keep <- md$nCount_Spatial >= input$min_counts & md$nCount_Spatial <= input$max_counts &
        md$nFeature_Spatial >= input$min_features
      if ("percent.mt" %in% names(md)) keep <- keep & md$percent.mt <= input$max_mt
      keep[is.na(keep)] <- FALSE
      rv$obj <- subset(obj, cells = rownames(md)[keep])
      log_message(rv, sprintf("QC filters applied; %d spots remain.", ncol(rv$obj)))
    }, error = function(e) log_message(rv, paste("FILTER ERROR:", conditionMessage(e))))
  })

  observeEvent(input$run_norm, {
    obj <- require_object(rv)
    started <- Sys.time()
    tryCatch({
      withProgress(message = "Normalizing Visium spots", value = 0, {
        assay <- spatial_assay(obj)
        if (!assay %in% Assays(obj)) stop(sprintf("Active expression assay '%s' is missing.", assay), call. = FALSE)
        counts <- get_counts_matrix(obj, assay)
        if (nrow(counts) < 3 || ncol(counts) < 3) stop("The expression matrix must contain at least three genes and three spots.", call. = FALSE)
        max_pcs <- max(2, min(input$n_pcs, nrow(counts) - 1, ncol(counts) - 1))
        nfeatures <- min(input$variable_features, nrow(counts))

        if (identical(input$normalization, "sct")) {
          if (!requireNamespace("sctransform", quietly = TRUE)) {
            stop("SCTransform requires the sctransform package. Install it with install.packages('sctransform').", call. = FALSE)
          }
          incProgress(0.1, detail = "Running SCTransform")
          obj <- SCTransform(
            object = obj,
            assay = assay,
            new.assay.name = "SCT",
            variable.features.n = nfeatures,
            return.only.var.genes = FALSE,
            vst.flavor = "v2",
            conserve.memory = TRUE,
            verbose = FALSE
          )
          DefaultAssay(obj) <- "SCT"
          pca_assay <- "SCT"
        } else {
          incProgress(0.1, detail = "Running LogNormalize")
          obj <- NormalizeData(
            object = obj,
            assay = assay,
            normalization.method = "LogNormalize",
            scale.factor = 10000,
            verbose = FALSE
          )
          obj <- FindVariableFeatures(
            object = obj,
            assay = assay,
            selection.method = "vst",
            nfeatures = nfeatures,
            verbose = FALSE
          )
          incProgress(0.45, detail = "Selecting variable genes")
          vars <- VariableFeatures(obj, assay = assay)
          if (length(vars) < 2) stop("Fewer than two variable genes were found after LogNormalize.", call. = FALSE)
          obj <- ScaleData(object = obj, assay = assay, features = vars, verbose = FALSE)
          DefaultAssay(obj) <- assay
          pca_assay <- assay
        }

        pca_features <- VariableFeatures(obj, assay = pca_assay)
        if (length(pca_features) < 2) stop("Fewer than two variable genes are available for PCA.", call. = FALSE)
        incProgress(0.75, detail = "Running PCA")
        obj <- RunPCA(
          object = obj,
          assay = pca_assay,
          features = pca_features,
          npcs = max_pcs,
          verbose = FALSE
        )
        rv$obj <- obj
        log_message(rv, sprintf("%s normalization and PCA completed using assay '%s' and %d PCs in %.1f seconds.", if (input$normalization == "sct") "SCTransform" else "LogNormalize", pca_assay, max_pcs, as.numeric(difftime(Sys.time(), started, units = "secs"))))
      })
    }, error = function(e) {
      log_message(rv, paste("NORMALIZATION ERROR:", conditionMessage(e)))
      log_message(rv, "Try LogNormalize first if SCTransform is unavailable or fails on the installed sctransform version.")
    })
  })

  output$elbow <- renderPlot({ obj <- require_object(rv); validate(need("pca" %in% Reductions(obj), "Run normalization first.")); ElbowPlot(obj, ndims = min(input$n_pcs, 50)) })
  output$pca_plot <- renderPlotly({ obj <- require_object(rv); validate(need("pca" %in% Reductions(obj), "Run normalization first.")); p <- DimPlot(obj, reduction = "pca", group.by = if ("seurat_clusters" %in% colnames(obj@meta.data)) "seurat_clusters" else NULL); ggplotly(p) })

  observeEvent(input$run_clusters, {
    obj <- require_object(rv)
    started <- Sys.time()
    tryCatch({
      validate(need("pca" %in% Reductions(obj), "Run normalization and PCA before clustering."))
      emb <- Embeddings(obj, "pca")
      if (nrow(emb) < 5) stop("At least five spots are required for clustering.", call. = FALSE)
      max_dim <- min(ncol(emb), nrow(emb) - 1)
      low <- max(1, min(input$dims_low, max_dim))
      high <- max(low, min(input$dims_high, max_dim))
      dims <- seq.int(low, high)
      k_use <- min(max(5, input$umap_neighbors), nrow(obj) - 1)

      withProgress(message = "Clustering Visium spots", value = 0, {
        incProgress(0.1, detail = sprintf("Finding neighbors with PCs %d-%d", low, high))
        obj <- FindNeighbors(
          object = obj,
          reduction = "pca",
          dims = dims,
          k.param = k_use,
          nn.method = if (isTRUE(input$fast_clusters)) "annoy" else "rann",
          n.trees = if (isTRUE(input$fast_clusters)) 20 else 50,
          verbose = FALSE
        )
        incProgress(0.45, detail = "Detecting clusters")
        obj <- FindClusters(object = obj, resolution = input$resolution, algorithm = 1, random.seed = 1234, verbose = FALSE)
        incProgress(0.65, detail = "Computing UMAP")
        obj <- RunUMAP(
          object = obj,
          reduction = "pca",
          dims = dims,
          metric = input$umap_metric,
          n.neighbors = k_use,
          min.dist = 0.3,
          n.components = 2,
          n.epochs = if (isTRUE(input$fast_clusters)) 200 else 500,
          verbose = FALSE,
          seed.use = 1234
        )
        rv$obj <- obj
        incProgress(1, detail = "Finished")
      })
      ids <- sort(unique(as.character(obj$seurat_clusters)))
      updateSelectInput(session, "marker_group", choices = unique(c("seurat_clusters", colnames(obj@meta.data))), selected = "seurat_clusters")
      updateSelectInput(session, "de_group1", choices = ids, selected = ids[1])
      updateSelectInput(session, "de_group2", choices = ids, selected = if (length(ids) > 1) ids[2] else ids[1])
      log_message(rv, sprintf("Neighbors, UMAP, and clusters completed: %d clusters in %.1f seconds.", length(ids), as.numeric(difftime(Sys.time(), started, units = "secs"))))
    }, error = function(e) {
      log_message(rv, paste("CLUSTER ERROR:", conditionMessage(e)))
      log_message(rv, "Check that PCA completed and that the selected PC range is within the available dimensions.")
    })
  })

  output$umap_plot <- renderPlotly({ obj <- require_object(rv); validate(need("umap" %in% Reductions(obj), "Run clustering first.")); ggplotly(DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)) })
  output$cluster_table <- renderDT({ obj <- require_object(rv); validate(need("seurat_clusters" %in% colnames(obj@meta.data), "Run clustering first.")); datatable(as.data.frame(table(cluster = obj$seurat_clusters)), options = list(pageLength = 15), rownames = FALSE) })

  output$spatial_plot <- renderPlotly({
    obj <- require_object(rv)
    validate(need(!is.null(rv$spatial_overlay), "Spatial coordinates are not loaded. Reload the four Visium files."))
    p <- tryCatch({
      assay <- spatial_assay(obj)
      if (identical(input$spatial_mode, "labels")) {
        validate(need(nzchar(input$spatial_group %||% ""), "Select a metadata column."))
        spatial_label_overlay_plot(obj, rv$spatial_overlay, input$spatial_group, input$spatial_point_size)
      } else {
        validate(need(nzchar(input$spatial_feature %||% ""), "Select a gene or feature."))
        spatial_overlay_plot(obj, rv$spatial_overlay, input$spatial_feature, assay, input$spatial_point_size, input$min_cutoff, input$max_cutoff)
      }
    }, error = function(e) {
      ggplot() + theme_void() + labs(title = paste("Spatial map error:", conditionMessage(e)))
    })
    ggplotly(p, tooltip = "text") %>% layout(legend = list(orientation = "v"))
  })

  output$cluster_spatial_plot <- renderPlotly({
    obj <- require_object(rv)
    validate(need(!is.null(rv$spatial_overlay), "Spatial coordinates are not loaded. Reload the four Visium files."))
    validate(need("seurat_clusters" %in% colnames(obj@meta.data), "Run clustering first."))
    p <- tryCatch(
      spatial_label_overlay_plot(obj, rv$spatial_overlay, "seurat_clusters", input$spatial_point_size %||% 1.6),
      error = function(e) ggplot() + theme_void() + labs(title = paste("Cluster map error:", conditionMessage(e)))
    )
    ggplotly(p, tooltip = "text")
  })

  output$feature_plot <- renderPlotly({
    obj <- require_object(rv); validate(need(input$explore_feature %in% available_features(obj, DefaultAssay(obj)), "Feature is not available in the active assay."));
    p <- FeaturePlot(obj, features = input$explore_feature, reduction = input$explore_reduction); ggplotly(p)
  })

  observeEvent(input$score_signature, {
    obj <- require_object(rv)
    genes <- unique(trimws(unlist(strsplit(input$signature_genes, "[,\\n;]+"))))
    genes <- genes[nzchar(genes)]; assay <- spatial_assay(obj); present <- intersect(genes, available_features(obj, assay))
    if (length(present) < 2) return(log_message(rv, "SIGNATURE ERROR: fewer than two genes were found in the spatial assay."))
    tryCatch({
      DefaultAssay(obj) <- assay; score_name <- paste0("Signature_", safe_filename(input$signature_name))
      obj <- AddModuleScore(obj, features = list(present), name = score_name, assay = assay)
      score_col <- paste0(score_name, "1"); rv$obj <- obj; rv$signatures[[input$signature_name]] <- list(genes = present, column = score_col)
      log_message(rv, paste("Signature scored:", input$signature_name, "using", length(present), "genes."))
    }, error = function(e) log_message(rv, paste("SIGNATURE ERROR:", conditionMessage(e))))
  })

  output$signature_table <- renderDT({
    if (!length(rv$signatures)) return(datatable(data.frame(Message = "No signatures scored yet."), options = list(dom = "t")))
    datatable(do.call(rbind, lapply(names(rv$signatures), function(n) data.frame(signature = n, genes = paste(rv$signatures[[n]]$genes, collapse = ", "), score_column = rv$signatures[[n]]$column))), options = list(pageLength = 10), rownames = FALSE)
  })

  marker_assay <- function(obj) {
    if ("SCT" %in% Assays(obj) && length(VariableFeatures(obj, assay = "SCT")) > 1) "SCT" else spatial_assay(obj)
  }

  get_groups <- function(obj, group_col) {
    if (is.null(group_col) || !nzchar(group_col) || !group_col %in% colnames(obj@meta.data)) {
      stop("Choose a valid group column. Run clustering first to create seurat_clusters.", call. = FALSE)
    }
    values <- as.character(obj@meta.data[[group_col]])
    values[is.na(values) | !nzchar(values)] <- NA_character_
    tab <- table(values, useNA = "no")
    if (length(tab) < 2) stop(sprintf("Group column '%s' must contain at least two groups.", group_col), call. = FALSE)
    if (any(tab < 3)) stop(sprintf("Every group in '%s' needs at least three spots; smallest group has %d.", group_col, min(tab)), call. = FALSE)
    values
  }

  fast_de_table <- function(obj, assay, groups, group1, group2, min_pct, logfc_threshold, max_genes = 2000) {
    expr <- tryCatch(GetAssayData(obj, assay = assay, layer = "data"), error = function(e) GetAssayData(obj, assay = assay, slot = "data"))
    expr <- as.matrix(expr)
    genes <- VariableFeatures(obj, assay = assay)
    if (!length(genes)) genes <- rownames(expr)
    genes <- intersect(genes, rownames(expr))
    if (length(genes) > max_genes) {
      vars <- apply(expr[genes, , drop = FALSE], 1, var)
      genes <- genes[order(vars[genes], decreasing = TRUE)[seq_len(max_genes)]]
    }
    x <- expr[genes, groups == group1, drop = FALSE]
    y <- expr[genes, groups == group2, drop = FALSE]
    if (ncol(x) < 3 || ncol(y) < 3) stop("Each DE group must contain at least three spots.", call. = FALSE)
    mean_x <- rowMeans(x); mean_y <- rowMeans(y)
    pct_x <- rowMeans(x > 0); pct_y <- rowMeans(y > 0)
    var_x <- apply(x, 1, var); var_y <- apply(y, 1, var)
    se <- sqrt(var_x / ncol(x) + var_y / ncol(y)); se[se == 0] <- Inf
    t_stat <- (mean_x - mean_y) / se
    df <- (var_x / ncol(x) + var_y / ncol(y))^2 / ((var_x / ncol(x))^2 / max(1, ncol(x) - 1) + (var_y / ncol(y))^2 / max(1, ncol(y) - 1))
    p_val <- 2 * stats::pt(-abs(t_stat), pmax(df, 1))
    p_val[!is.finite(p_val)] <- 1
    out <- data.frame(p_val = p_val, avg_log2FC = log2((mean_x + 1e-6) / (mean_y + 1e-6)), pct.1 = pct_x, pct.2 = pct_y, row.names = genes)
    out$p_val_adj <- p.adjust(out$p_val, method = "BH")
    out <- out[out$pct.1 >= min_pct | out$pct.2 >= min_pct, , drop = FALSE]
    out <- out[abs(out$avg_log2FC) >= logfc_threshold, , drop = FALSE]
    out[order(out$p_val_adj, -abs(out$avg_log2FC)), , drop = FALSE]
  }

  observeEvent(input$marker_group, {
    if (is.null(rv$obj) || is.null(input$marker_group) || !input$marker_group %in% colnames(rv$obj@meta.data)) return()
    values <- unique(as.character(rv$obj@meta.data[[input$marker_group]]))
    values <- sort(values[!is.na(values) & nzchar(values)])
    if (length(values)) {
      updateSelectInput(session, "de_group1", choices = values, selected = values[1])
      updateSelectInput(session, "de_group2", choices = values, selected = if (length(values) > 1) values[2] else values[1])
    }
  }, ignoreInit = FALSE)

  observeEvent(input$run_markers, {
    obj <- require_object(rv)
    started <- Sys.time()
    tryCatch({
      groups <- get_groups(obj, input$marker_group)
      assay <- marker_assay(obj)
      DefaultAssay(obj) <- assay
      Idents(obj) <- groups
      rv$markers <- NULL
      rv$de <- NULL
      withProgress(message = "Finding markers", value = 0, {
        incProgress(0.15, detail = sprintf("Using assay '%s'", assay))
        if (isTRUE(input$fast_de)) {
          levels <- sort(unique(groups[!is.na(groups)]))
          marker_list <- lapply(levels, function(g) {
            other <- if (length(levels) == 2) levels[levels != g] else "__other__"
            comparison <- if (length(levels) == 2) groups else ifelse(groups == g, g, other)
            tbl <- fast_de_table(obj, assay, comparison, g, other, input$marker_min_pct, input$marker_logfc, input$de_max_genes)
            if (nrow(tbl)) tbl$cluster <- g
            tbl
          })
          rv$markers <- do.call(rbind, marker_list)
          if (is.null(rv$markers)) rv$markers <- data.frame()
        } else {
          rv$markers <- FindAllMarkers(
            object = obj, assay = assay, only.pos = TRUE, test.use = "wilcox",
            min.pct = input$marker_min_pct, logfc.threshold = input$marker_logfc,
            max.cells.per.ident = if (isTRUE(input$fast_clusters)) 5000 else Inf,
            random.seed = 1234, verbose = FALSE
          )
        }
        incProgress(1, detail = "Finished")
      })
      if (is.null(rv$markers) || !nrow(rv$markers)) stop("No marker genes passed the current thresholds. Lower min.pct or logFC threshold.", call. = FALSE)
      rv$markers$feature <- rownames(rv$markers)
      rv$markers <- rv$markers[order(rv$markers$p_val_adj, -abs(rv$markers$avg_log2FC)), , drop = FALSE]
      log_message(rv, sprintf("Marker analysis completed using '%s' (%s): %d rows in %.1f seconds.", assay, if (isTRUE(input$fast_de)) "fast vectorized" else "Seurat Wilcoxon", nrow(rv$markers), as.numeric(difftime(Sys.time(), started, units = "secs"))))
    }, error = function(e) log_message(rv, paste("MARKER ERROR:", conditionMessage(e))))
  })

  observeEvent(input$run_de, {
    obj <- require_object(rv)
    started <- Sys.time()
    tryCatch({
      groups <- get_groups(obj, input$marker_group)
      ids <- sort(unique(groups[!is.na(groups)]))
      if (is.null(input$de_group1) || !input$de_group1 %in% ids) stop("Select a valid Group 1.", call. = FALSE)
      if (is.null(input$de_group2) || !input$de_group2 %in% ids) stop("Select a valid Group 2.", call. = FALSE)
      if (identical(input$de_group1, input$de_group2)) stop("Group 1 and Group 2 must be different.", call. = FALSE)
      assay <- marker_assay(obj)
      DefaultAssay(obj) <- assay
      Idents(obj) <- groups
      rv$de <- NULL
      withProgress(message = "Running differential expression", value = 0, {
        incProgress(0.15, detail = sprintf("Comparing %s vs %s", input$de_group1, input$de_group2))
        if (isTRUE(input$fast_de)) {
          rv$de <- fast_de_table(obj, assay, groups, input$de_group1, input$de_group2, input$marker_min_pct, input$marker_logfc, input$de_max_genes)
        } else {
          rv$de <- FindMarkers(
            object = obj, assay = assay, ident.1 = input$de_group1, ident.2 = input$de_group2,
            test.use = "wilcox", min.pct = input$marker_min_pct,
            logfc.threshold = input$marker_logfc, random.seed = 1234, verbose = FALSE
          )
        }
        incProgress(1, detail = "Finished")
      })
      if (is.null(rv$de) || !nrow(rv$de)) stop("No differential genes passed the current thresholds. Lower min.pct or logFC threshold.", call. = FALSE)
      rv$de$feature <- rownames(rv$de)
      rv$de <- rv$de[order(rv$de$p_val_adj, -abs(rv$de$avg_log2FC)), , drop = FALSE]
      log_message(rv, sprintf("Pairwise DE completed using '%s' (%s): %d rows in %.1f seconds.", assay, if (isTRUE(input$fast_de)) "fast vectorized" else "Seurat Wilcoxon", nrow(rv$de), as.numeric(difftime(Sys.time(), started, units = "secs"))))
    }, error = function(e) log_message(rv, paste("DE ERROR:", conditionMessage(e))))
  })
  output$marker_table <- renderDT({
    tbl <- rv$de %||% rv$markers
    validate(need(!is.null(tbl), "Run a marker or pairwise DE analysis.")); datatable(as.data.frame(tbl), filter = "top", extensions = "Buttons", options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")), rownames = TRUE)
  })

  output$stie_install_status <- renderText(rv$stie_install)

  observeEvent(input$install_stie, {
    rv$stie_install <- "Installing quadprog, magick, remotes, and STIE..."
    tryCatch({
      cran_pkgs <- c("quadprog", "magick", "remotes")
      missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing_cran)) install.packages(missing_cran, repos = "https://cloud.r-project.org")
      remotes::install_github("zhushijia/STIE", upgrade = "never", dependencies = TRUE, quiet = TRUE)
      rv$stie_install <- "STIE and its dependencies installed. If R reports that a package is in use, restart R/Shiny and run the STIE analysis."
      log_message(rv, "STIE installation completed.")
    }, error = function(e) {
      rv$stie_install <- paste("STIE installation failed:", conditionMessage(e))
      log_message(rv, rv$stie_install)
    })
  })

  normalize_stie_barcode <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("-[0-9]+$", "", x)
    x
  }

  prepare_stie_cells <- function(cells, obj, morphology_text) {
    required <- c("cell_id", "spot")
    missing_required <- setdiff(required, names(cells))
    if (length(missing_required)) stop(sprintf("Cells-on-spot CSV is missing required column(s): %s.", paste(missing_required, collapse = ", ")), call. = FALSE)
    if (anyDuplicated(as.character(cells$cell_id))) stop("The cell_id column must contain unique cell identifiers.", call. = FALSE)
    cells$cell_id <- as.character(cells$cell_id)
    cells$spot <- as.character(cells$spot)
    features <- trimws(unlist(strsplit(morphology_text %||% "", "[,;]+")))
    features <- unique(features[nzchar(features)])
    if (!length(features)) stop("Enter at least one morphology column in the selected morphology columns field.", call. = FALSE)
    missing_features <- setdiff(features, names(cells))
    if (length(missing_features)) stop(sprintf("Selected morphology column(s) are missing from the CSV: %s.", paste(missing_features, collapse = ", ")), call. = FALSE)
    non_numeric <- features[!vapply(cells[, features, drop = FALSE], is.numeric, logical(1))]
    if (length(non_numeric)) stop(sprintf("Morphology column(s) must be numeric: %s.", paste(non_numeric, collapse = ", ")), call. = FALSE)
    visium_barcodes <- colnames(obj)
    exact <- cells$spot %in% visium_barcodes
    if (!any(exact)) {
      norm_visium <- normalize_stie_barcode(visium_barcodes)
      norm_cells <- normalize_stie_barcode(cells$spot)
      index <- match(norm_cells, norm_visium)
      exact <- !is.na(index)
      cells$spot <- ifelse(exact, visium_barcodes[index], cells$spot)
    }
    cells <- cells[exact, , drop = FALSE]
    if (!nrow(cells)) stop("No cells-on-spot rows matched the loaded Visium barcodes. Verify that the morphology CSV and H5 matrix come from the same library.", call. = FALSE)
    list(cells = cells, features = features)
  }

  observeEvent(input$run_stie, {
    obj <- require_object(rv)
    started <- Sys.time()
    tryCatch({
      if (!requireNamespace("STIE", quietly = TRUE)) stop("STIE is not installed. Click 'Install STIE dependencies', restart R/Shiny if requested, and run again.", call. = FALSE)
      if (is.null(input$stie_cells)) stop("Upload a cells-on-spot morphology CSV.", call. = FALSE)
      if (is.null(input$stie_signature)) stop("Upload an initial gene-by-cell-type signature CSV. STIE requires a signature for both supported modes.", call. = FALSE)
      cells <- as.data.frame(readr::read_csv(input$stie_cells$datapath, show_col_types = FALSE, progress = FALSE), check.names = FALSE)
      prepared <- prepare_stie_cells(cells, obj, input$stie_features)
      cells <- prepared$cells
      features <- prepared$features
      signature <- read_stie_matrix(input$stie_signature$datapath)
      assay <- if ("SCT" %in% Assays(obj) && length(VariableFeatures(obj, assay = "SCT")) > 1) "SCT" else spatial_assay(obj)
      expr <- tryCatch(GetAssayData(obj, assay = assay, layer = "data"), error = function(e) GetAssayData(obj, assay = assay, slot = "data"))
      st_expr <- t(as.matrix(expr))
      common_spots <- intersect(rownames(st_expr), as.character(cells$spot))
      if (length(common_spots) < 10) stop(sprintf("Only %d matching spots were found between Seurat and cells-on-spot data; at least 10 are required.", length(common_spots)), call. = FALSE)
      cells <- cells[as.character(cells$spot) %in% common_spots, , drop = FALSE]
      st_expr <- st_expr[common_spots, , drop = FALSE]
      signature <- signature[intersect(rownames(signature), colnames(st_expr)), , drop = FALSE]
      if (nrow(signature) < 10) stop("Fewer than 10 signature genes overlap the Visium expression matrix.", call. = FALSE)
      st_expr <- st_expr[, rownames(signature), drop = FALSE]
      withProgress(message = "Running STIE", value = 0, {
        incProgress(0.1, detail = sprintf("Using %d spots, %d cells, and %d genes", nrow(st_expr), nrow(cells), nrow(signature)))
        result <- STIE::STIE(
          ST_expr = st_expr,
          Signature = signature,
          cells_on_spot = cells,
          features = features,
          lambda = input$stie_lambda,
          steps = input$stie_steps,
          morphology_steps = max(1, ceiling(input$stie_steps / 3)),
          known_signature = identical(input$stie_mode, "deconv"),
          known_cell_types = FALSE,
          min_cells = input$stie_min_cells,
          equal_prior = isTRUE(input$stie_equal_prior),
          verbose = FALSE
        )
        incProgress(0.9, detail = "Attaching STIE probabilities")
        probs <- as.data.frame(result$PE_on_spot)
        probs <- probs[intersect(rownames(probs), colnames(obj)), , drop = FALSE]
        names(probs) <- paste0("STIE_spot_", safe_filename(names(probs)))
        obj <- AddMetaData(obj, metadata = probs)
        rv$obj <- obj
        # Merge STIE's per-cell type predictions back into the cells-on-spot
        # table so the tissue cell overlay can color individual cells by
        # their predicted type, not just show spot-level probabilities.
        cells_with_types <- cells
        if (!is.null(result$cell_types)) {
          ct <- result$cell_types
          cells_with_types$predicted_cell_type <- as.character(ct[match(cells_with_types$cell_id, names(ct))])
        }
        rv$stie <- list(result = result, probabilities = probs, assay = assay, cells = cells_with_types)
        incProgress(1, detail = "Finished")
      })
      updateSelectInput(session, "stie_prob_column", choices = names(probs), selected = names(probs)[1])
      log_message(rv, sprintf("STIE completed using %d matching spots in %.1f seconds.", length(common_spots), as.numeric(difftime(Sys.time(), started, units = "secs"))))
    }, error = function(e) {
      log_message(rv, paste("STIE ERROR:", conditionMessage(e)))
    })
  })

  output$stie_cells_spatial <- renderPlotly({
    obj <- require_object(rv)
    validate(need(!is.null(rv$stie), "Run STIE first."))
    validate(need(!is.null(rv$stie$cells), "STIE cell data are unavailable."))
    p <- tryCatch(
      cell_overlay_plot(obj, rv$spatial_overlay, rv$stie$cells, max(0.8, (input$spatial_point_size %||% 1.6) * 0.7)),
      error = function(e) ggplot() + theme_void() + labs(title = paste("Cell overlay error:", conditionMessage(e)))
    )
    ggplotly(p, tooltip = "text")
  })

  output$stie_table <- renderDT({
    validate(need(!is.null(rv$stie), "Run STIE first."))
    result <- rv$stie$result
    cell_types <- result$cell_types
    if (is.null(cell_types)) return(datatable(data.frame(Message = "STIE returned no cell-level assignments."), options = list(dom = "t")))
    tbl <- data.frame(cell_id = names(cell_types), predicted_cell_type = as.character(cell_types), stringsAsFactors = FALSE)
    datatable(tbl, filter = "top", options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  })

  output$stie_composition_plot <- renderPlotly({
    validate(need(!is.null(rv$stie), "Run STIE first."))
    comp <- cell_type_composition(rv$stie$cells)
    validate(need(!is.null(comp), "STIE did not return per-cell type predictions to summarize."))
    comp <- comp[order(-comp$n_cells), , drop = FALSE]
    comp$cell_type <- factor(comp$cell_type, levels = comp$cell_type)
    p <- ggplot(comp, aes(x = cell_type, y = n_cells, fill = cell_type,
                          text = sprintf("%s<br>Cells: %d (%.1f%%)", cell_type, n_cells, 100 * fraction))) +
      geom_col() +
      scale_fill_manual(values = ivisio_discrete_palette(nrow(comp)), guide = "none") +
      labs(x = NULL, y = "Cells", title = "Predicted cell-type composition") +
      theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
    ggplotly(p, tooltip = "text")
  })

  output$stie_spatial <- renderPlotly({
    obj <- require_object(rv)
    validate(need(!is.null(rv$stie), "Run STIE first."))
    probs <- rv$stie$probabilities
    validate(need(ncol(probs) > 0, "STIE returned no spot probabilities."))
    feature <- if (!is.null(input$stie_prob_column) && nzchar(input$stie_prob_column) && input$stie_prob_column %in% names(probs)) input$stie_prob_column else names(probs)[1]
    p <- tryCatch(
      spatial_values_overlay_plot(
        rv$spatial_overlay, rownames(probs), probs[[feature]], legend_name = feature,
        title = sprintf("iVisio STIE probability: %s", feature),
        point_size = input$spatial_point_size %||% 1.6, min_cutoff = input$min_cutoff %||% "q01", max_cutoff = input$max_cutoff %||% "q99"
      ),
      error = function(e) ggplot() + theme_void() + labs(title = paste("STIE spatial map error:", conditionMessage(e)))
    )
    ggplotly(p, tooltip = "text")
  })

  output$reference_catalog_table <- renderDT({
    datatable(reference_catalog, filter = "top", options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })

  observeEvent(input$load_reference, {
    tryCatch({
      if (is.null(input$reference_file) && is.null(input$reference_signature)) stop("Upload a reference RDS/H5/CSV or a precomputed signature CSV.", call. = FALSE)
      if (!is.null(input$reference_signature)) {
        rv$reference_signature <- read_stie_matrix(input$reference_signature$datapath)
        rv$reference <- list(source = input$reference_signature$name, type = "precomputed signature")
      } else {
        rv$reference <- read_reference_expression(input$reference_file$datapath)
        rv$reference$source <- input$reference_file$name
        if (!is.null(rv$reference$groups)) {
          rv$reference_signature <- build_reference_signature(rv$reference, rv$reference$groups)
        } else {
          stop("The uploaded reference does not contain cell-type labels. Upload a precomputed signature CSV instead.", call. = FALSE)
        }
      }
      log_message(rv, sprintf("Reference signature loaded: %d genes x %d cell types.", nrow(rv$reference_signature), ncol(rv$reference_signature)))
    }, error = function(e) log_message(rv, paste("REFERENCE ERROR:", conditionMessage(e))))
  })

  output$reference_signature_table <- renderDT({
    validate(need(!is.null(rv$reference_signature), "Load a reference first."))
    sig <- as.data.frame(rv$reference_signature)
    sig$gene <- rownames(sig)
    datatable(head(sig, 100), filter = "top", options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  })

  observeEvent(input$run_reference_deconv, {
    obj <- require_object(rv)
    tryCatch({
      validate(need(!is.null(rv$reference_signature), "Load a reference signature first."))
      assay <- if ("SCT" %in% Assays(obj) && length(VariableFeatures(obj, assay = "SCT")) > 1) "SCT" else spatial_assay(obj)
      expr <- tryCatch(GetAssayData(obj, assay = assay, layer = "data"), error = function(e) GetAssayData(obj, assay = assay, slot = "data"))
      st_expr <- t(as.matrix(expr))
      withProgress(message = "Running reference deconvolution", value = 0, {
        incProgress(0.2, detail = sprintf("Using assay '%s'", assay))
        probs <- simple_nnls_deconvolution(st_expr, rv$reference_signature, max_iter = 300)
        names(probs) <- paste0("Reference_", safe_filename(names(probs)))
        probs <- probs[intersect(rownames(probs), colnames(obj)), , drop = FALSE]
        obj <- AddMetaData(obj, metadata = as.data.frame(probs))
        rv$obj <- obj
        rv$reference_deconv <- probs
        incProgress(1, detail = "Finished")
      })
      updateSelectInput(session, "reference_cell_type", choices = names(rv$reference_deconv), selected = names(rv$reference_deconv)[1])
      log_message(rv, sprintf("Reference deconvolution completed for %d spots and %d cell types.", nrow(rv$reference_deconv), ncol(rv$reference_deconv)))
    }, error = function(e) log_message(rv, paste("DECONVOLUTION ERROR:", conditionMessage(e))))
  })

  output$reference_deconv_table <- renderDT({
    validate(need(!is.null(rv$reference_deconv), "Run reference deconvolution first."))
    tbl <- as.data.frame(rv$reference_deconv); tbl$barcode <- rownames(tbl)
    datatable(head(tbl, 500), filter = "top", options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  })

  output$reference_spatial_plot <- renderPlotly({
    validate(need(!is.null(rv$reference_deconv), "Run reference deconvolution first."))
    validate(need(!is.null(rv$spatial_overlay), "Spatial coordinates are not loaded."))
    probs <- rv$reference_deconv
    cell_type <- if (!is.null(input$reference_cell_type) && nzchar(input$reference_cell_type) && input$reference_cell_type %in% names(probs)) input$reference_cell_type else names(probs)[1]
    p <- tryCatch(
      spatial_values_overlay_plot(
        rv$spatial_overlay, rownames(probs), probs[[cell_type]], legend_name = cell_type,
        title = sprintf("iVisio reference proportion: %s", cell_type),
        point_size = input$spatial_point_size %||% 1.6, min_cutoff = "0", max_cutoff = "1"
      ),
      error = function(e) ggplot() + theme_void() + labs(title = paste("Reference spatial map error:", conditionMessage(e)))
    )
    ggplotly(p, tooltip = "text")
  })

  output$reference_composition_plot <- renderPlotly({
    validate(need(!is.null(rv$reference_deconv), "Run reference deconvolution first."))
    probs <- rv$reference_deconv
    comp <- data.frame(cell_type = names(probs), mean_proportion = colMeans(probs, na.rm = TRUE), stringsAsFactors = FALSE)
    comp <- comp[order(-comp$mean_proportion), , drop = FALSE]
    comp$cell_type <- factor(comp$cell_type, levels = comp$cell_type)
    p <- ggplot(comp, aes(x = cell_type, y = mean_proportion, fill = cell_type,
                          text = sprintf("%s<br>Mean proportion: %.1f%%", cell_type, 100 * mean_proportion))) +
      geom_col() +
      scale_fill_manual(values = ivisio_discrete_palette(nrow(comp)), guide = "none") +
      scale_y_continuous(labels = function(v) sprintf("%d%%", round(100 * v))) +
      labs(x = NULL, y = "Mean proportion across spots", title = "Average tissue composition") +
      theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
    ggplotly(p, tooltip = "text")
  })

  observeEvent(input$run_svg, {
    obj <- require_object(rv)
    tryCatch({
      assay <- spatial_assay(obj); DefaultAssay(obj) <- assay
      svg_features <- FindSpatiallyVariableFeatures(obj, assay = assay, selection.method = "markvariogram", features = head(VariableFeatures(obj), input$svg_n), verbose = FALSE)
      rv$svg <- if (is.data.frame(svg_features)) svg_features else data.frame(feature = as.character(svg_features), row.names = NULL)
      log_message(rv, "Spatially variable feature analysis completed.")
    }, error = function(e) log_message(rv, paste("SVG ERROR:", conditionMessage(e))))
  })
  output$svg_table <- renderDT({ validate(need(!is.null(rv$svg), "Run spatial variability analysis first.")); datatable(as.data.frame(rv$svg), options = list(pageLength = 25, scrollX = TRUE), rownames = TRUE) })
  output$svg_plot <- renderPlot({ obj <- require_object(rv); validate(need(!is.null(rv$svg) && nrow(rv$svg) > 0, "Run spatial variability analysis first.")); top <- if ("feature" %in% names(rv$svg)) rv$svg$feature[1] else rownames(rv$svg)[1]; spatial_plot(obj, top, input$image_name, input$spatial_point_size, 0.9, input$min_cutoff, input$max_cutoff) })

  output$export_status <- renderPrint({
    cat("Output prefix:", safe_filename(input$out_prefix), "\n")
    if (is.null(rv$obj)) return(cat("Nothing to export. Load Visium data first.\n"))
    cat("Spots:", ncol(rv$obj), "\nMarkers:", !is.null(rv$markers), "\nDE:", !is.null(rv$de), "\nSVG:", !is.null(rv$svg), "\nReport:", !is.null(rv$report), "\n")
    if (!is.null(rv$report)) cat("Report file:", rv$report, "\n")
  })
  output$download_rds <- downloadHandler(filename = function() paste0(safe_filename(input$out_prefix), ".rds"), content = function(file) { req(rv$obj); saveRDS(rv$obj, file) })
  output$download_markers <- downloadHandler(filename = function() paste0(safe_filename(input$out_prefix), ".markers.csv"), content = function(file) { if (is.null(rv$markers) && is.null(rv$de)) stop("Run marker or DE analysis first.", call. = FALSE); write.csv(rv$markers %||% rv$de, file, row.names = TRUE) })
  output$download_svg <- downloadHandler(filename = function() paste0(safe_filename(input$out_prefix), ".spatially_variable.csv"), content = function(file) { req(rv$svg); write.csv(rv$svg, file, row.names = TRUE) })
  output$download_reference_deconv <- downloadHandler(filename = function() paste0(safe_filename(input$out_prefix), ".reference_deconvolution.csv"), content = function(file) { req(rv$reference_deconv); write.csv(rv$reference_deconv, file, row.names = TRUE) })
  output$download_reference_signature <- downloadHandler(filename = function() paste0(safe_filename(input$out_prefix), ".reference_signature.csv"), content = function(file) { req(rv$reference_signature); write.csv(rv$reference_signature, file, row.names = TRUE) })
  output$download_metadata <- downloadHandler(filename = function() paste0(safe_filename(input$out_prefix), ".spot_metadata.csv"), content = function(file) { req(rv$obj); write.csv(rv$obj@meta.data, file, row.names = TRUE) })

  observeEvent(input$make_report, {
    req(rv$obj)
    tryCatch({
      prefix <- safe_filename(input$out_prefix)
      report_dir <- file.path(tempdir(), paste0("ivisio_report_", as.integer(Sys.time())))
      dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
      path <- file.path(report_dir, paste0(prefix, "_report.html"))
      rmd <- file.path(report_dir, "ivisio_report.Rmd")
      summary <- paste0("Project: ", rv$obj@project.name, "\nSpots: ", ncol(rv$obj), "\nFeatures: ", nrow(rv$obj[[spatial_assay(rv$obj)]]), "\nAssays: ", paste(Assays(rv$obj), collapse = ", "))
      writeLines(c("---", "title: 'iVisio Spatial Omics Report'", "output:", "  html_document:", "    self_contained: true", "params:", "  summary: ''", "---", "", "## Summary", "```{r echo=FALSE}", "cat(params$summary)", "```", "", "## Session information", "```{r echo=FALSE}", "sessionInfo()", "```"), rmd)
      if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("The rmarkdown package is required to create the HTML report.", call. = FALSE)
      rmarkdown::render(rmd, output_file = basename(path), output_dir = report_dir, params = list(summary = summary), quiet = TRUE, envir = new.env(parent = globalenv()))
      if (!file.exists(path)) stop("rmarkdown completed but the HTML report file was not created.", call. = FALSE)
      rv$report <- normalizePath(path, mustWork = TRUE); log_message(rv, paste("Report created:", rv$report))
    }, error = function(e) log_message(rv, paste("REPORT ERROR:", conditionMessage(e))))
  })
  output$download_report <- downloadHandler(filename = function() paste0(safe_filename(input$out_prefix), "_report.html"), content = function(file) { req(rv$report); ok <- file.copy(rv$report, file, overwrite = TRUE); if (!ok) stop("The HTML report could not be copied to the download location.", call. = FALSE) })
}

shinyApp(ui, server)

# Suggested installation commands, run once in a separate R session:
# install.packages(c("shiny", "bslib", "shinycssloaders", "DT", "plotly", "ggplot2", "dplyr", "readr", "rmarkdown", "tidyr"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("Seurat", "SeuratObject", "sctransform"))
# Run with: shiny::runApp("app.R")

# Important: browser uploads of folders are flattened by many Shiny deployments. For reliable
# production use, zip the complete Space Ranger output, unpack it server-side, then call read_visium().
# Shiny's internal limit is disabled, but configure your reverse proxy as well. For nginx, use:
#   client_max_body_size 0;
# or use a large finite value such as:
#   client_max_body_size 1024G;
# For Apache, raise LimitRequestBody accordingly. Posit Connect, shinyapps.io,
# Kubernetes ingress, and cloud load balancers have separate platform limits.
# This app intentionally does not perform package installation or execute uploaded files.

# End of app.R
