# Seurat objects written by R itself, for scripts/test-seurat.mjs.
#
#   Rscript scripts/make-seurat-fixtures.R scripts/fixtures
#
# The reader's other tests write R's format by hand, which proves the reader
# against its author's idea of the format. These are the format as saveRDS()
# writes it, from the objects Seurat itself builds, and beside each the names
# and sums R gives for it (seurat-expected.json), so the test compares the lab
# with R and not with itself. Two things came from one real object the lab
# refused (a Seurat 5 object, November 2026):
#   - a function among a logged command's parameters. Its environment is a
#     package namespace, which R writes as a bare string vector, not as an item;
#   - a v5 assay (Assay5). Its layers carry no names: the assay's `features`
#     and `cells` maps do.
# Small on purpose: 40 genes, 24 cells. Needs Seurat 5 and jsonlite.

suppressMessages({ library(Seurat); library(SeuratObject); library(Matrix) })
out <- commandArgs(TRUE)[1]
if (is.na(out)) stop("usage: Rscript scripts/make-seurat-fixtures.R <fixtures dir>")
set.seed(20261114)

ng <- 40L; nc <- 24L
counts <- Matrix(rpois(ng * nc, 1.2), ng, nc, sparse = TRUE)
rownames(counts) <- sprintf("GENE%02d", seq_len(ng))
colnames(counts) <- sprintf("%s_CELL%02d", rep(c("S1", "S2"), each = nc / 2), seq_len(nc))

build <- function(version) {
  options(Seurat.object.assay.version = version)
  obj <- CreateSeuratObject(counts, min.cells = 0, min.features = 0)
  obj$sample <- rep(c("S1", "S2"), each = nc / 2)
  obj$cell_type <- factor(rep(c("A", "B", "C"), length.out = nc))
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 10, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 5, verbose = FALSE, approx = FALSE)
  emb <- Embeddings(obj, "pca")[, 1:2]
  colnames(emb) <- c("umap_1", "umap_2")
  obj[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "umap_", assay = "RNA")
  obj
}

# What a logged command holds when a function was one of its arguments: a
# closure from a package (its environment is that package's namespace), a
# builtin, a special, and the two kinds of package environment themselves.
with_functions <- function(obj) {
  cmd <- obj@commands[[1]]
  cmd@params$fn <- stats::median
  cmd@params$builtin <- sum
  cmd@params$special <- quote
  cmd@params$namespace <- asNamespace("stats")
  cmd@params$package <- as.environment("package:stats")
  cmd@params$again <- stats::sd      # the same namespace a second time: a reference
  obj@commands[[1]] <- cmd
  obj
}

sparse_sums <- function(m) list(
  genes = rownames(m), cells = colnames(m), dim = dim(m),
  col_sums = as.numeric(Matrix::colSums(m)), row_sums = as.numeric(Matrix::rowSums(m)),
  nnz = length(m@x))

expected <- list(
  written_by = list(R = R.version.string,
                    Seurat = as.character(packageVersion("Seurat")),
                    SeuratObject = as.character(packageVersion("SeuratObject"))))

save <- function(obj, name) {
  saveRDS(obj, file.path(out, name))
  layers <- list()
  a <- obj[["RNA"]]
  for (l in Layers(a)) {
    m <- LayerData(a, layer = l)
    if (inherits(m, "dgCMatrix")) layers[[l]] <- sparse_sums(m)
  }
  expected[[name]] <<- list(
    class = class(obj)[1], assay_class = class(a)[1],
    cells = rownames(obj@meta.data), layers = layers,
    meta_columns = colnames(obj@meta.data),
    commands = names(obj@commands), reductions = names(obj@reductions))
}

# 1. A v5 object as Seurat 5 makes it: joined layers, names in the assay's maps.
v5 <- build("v5")
save(v5, "seurat-v5-plain.rds")

# 2. The same, with functions among a command's parameters, and one more sparse
#    layer that holds THREE of the forty genes, handed over in another order
#    than the assay's: SeuratObject stores it in the assay's order.
v5f <- with_functions(v5)
a <- v5f[["RNA"]]
LayerData(a, layer = "sub") <- counts[c(31, 4, 17), ]
v5f[["RNA"]] <- a
save(v5f, "seurat-v5-functions.rds")

# 3. Layers split per sample (what merge() and split() leave, before JoinLayers).
v5s <- build("v5")
a <- v5s[["RNA"]]
a <- split(a, f = v5s$sample)
v5s[["RNA"]] <- a
save(v5s, "seurat-v5-split.rds")

# 4. meta.data in another order than the layers' cells (set past the validity
#    check, as a hand edit would): no layer is then the object's matrix.
v5r <- build("v5")
v5r@meta.data <- v5r@meta.data[rev(seq_len(nc)), ]
save(v5r, "seurat-v5-reordered.rds")

# 5. The older assay (v3/v4: names on each matrix), with the same functions.
v4f <- with_functions(build("v3"))
save(v4f, "seurat-v4-functions.rds")

jsonlite::write_json(expected, file.path(out, "seurat-expected.json"),
                     auto_unbox = TRUE, digits = NA, pretty = TRUE)
cat("wrote", length(expected) - 1, "objects and seurat-expected.json to", out, "\n")
