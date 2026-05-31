sourced_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (!is.null(sourced_file)) dirname(normalizePath(sourced_file, mustWork = FALSE)) else getwd()
dir_hUSI <- normalizePath(dirname(script_dir), winslash = "/", mustWork = FALSE)
load(file.path(dir_hUSI, 'sc/SenOCLR_mouse_map 250729.rdata'))

cal_mUSI <- function(mat_mouse){
  genes <- intersect(names(w_m), rownames(mat_mouse))
  if(length(genes) < 50) stop("Too few orthologs shared.")
  cor_vec <- apply(mat_mouse[genes, ], 2,
                   \(s) cor(w_m[genes], s, method = "spearman"))
  setNames(minmax(cor_vec), colnames(mat_mouse))
}
