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

build <- function(version, m = counts, assay = "RNA") {
  options(Seurat.object.assay.version = version)
  obj <- CreateSeuratObject(m, assay = assay, min.cells = 0, min.features = 0)
  obj$sample <- rep(c("S1", "S2"), each = nc / 2)[seq_len(ncol(m))]
  obj$cell_type <- factor(rep(c("A", "B", "C"), length.out = ncol(m)))
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 10, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 5, verbose = FALSE, approx = FALSE)
  emb <- Embeddings(obj, "pca")[, 1:2]
  colnames(emb) <- c("umap_1", "umap_2")
  obj[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "umap_", assay = assay)
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
  genes = I(rownames(m)), cells = I(colnames(m)), dim = dim(m),
  col_sums = as.numeric(Matrix::colSums(m)), row_sums = as.numeric(Matrix::rowSums(m)),
  nnz = length(m@x))

expected <- list(
  written_by = list(R = R.version.string,
                    Seurat = as.character(packageVersion("Seurat")),
                    SeuratObject = as.character(packageVersion("SeuratObject"))))

save <- function(obj, name, assay = "RNA") {
  saveRDS(obj, file.path(out, name))
  layers <- list()
  a <- obj[[assay]]
  for (l in Layers(a)) {
    m <- LayerData(a, layer = l)
    if (inherits(m, "dgCMatrix")) layers[[l]] <- sparse_sums(m)
  }
  expected[[name]] <<- list(
    class = class(obj)[1], assay_class = class(a)[1],
    cells = I(rownames(obj@meta.data)), layers = layers,
    meta_columns = I(colnames(obj@meta.data)),
    commands = I(as.character(names(obj@commands))),
    reductions = I(as.character(names(obj@reductions))))
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

# ---- what an independent check of the first fix found (2026-11-14) ---------
# Names made from numbers. as.character(1:24) is kept by R as the numbers and a
# promise to print them (a "deferred string"), and is saved that way. A reader
# that hands back the numbers has no names, and took a split layer of 12 cells
# for an object of 24.
int_cells <- counts
colnames(int_cells) <- as.character(seq_len(nc))
int_both <- int_cells
rownames(int_both) <- as.character(seq_len(ng))

# 6. Integer-made cell names, layers split per sample.
v5is <- build("v5", m = int_cells)
a <- v5is[["RNA"]]; a <- split(a, f = v5is$sample); v5is[["RNA"]] <- a
save(v5is, "seurat-v5-intnames-split.rds")

# 7. Integer-made cell names and gene names, joined.
save(build("v5", m = int_both), "seurat-v5-intnames.rds")

# 8. A second assay that covers only half the cells (twelve of twenty-four),
#    with ONE layer: it is not split, and no note may say it is.
c1 <- build("v5", m = counts[, 1:12])
adt <- Matrix(rpois(5 * 12, 3), 5, 12, sparse = TRUE)
rownames(adt) <- sprintf("ADT%d", 1:5); colnames(adt) <- colnames(counts)[1:12]
c1[["ADT"]] <- CreateAssay5Object(adt)
c2 <- build("v5", m = counts[, 13:24])
mg <- merge(c1, c2)
a <- mg[["RNA"]]; a <- JoinLayers(a); mg[["RNA"]] <- a
stopifnot(length(Layers(mg[["ADT"]])) == 1, ncol(mg[["ADT"]]) == 12, ncol(mg) == 24)
save(mg, "seurat-v5-partial-assay.rds")

# 9. Split layers in an assay that is NOT called RNA, and the same object after
#    the very lines of R the lab's refusal gives for it.
sp <- build("v5", assay = "Spatial")
a <- sp[["Spatial"]]; a <- split(a, f = sp$sample); sp[["Spatial"]] <- a
save(sp, "seurat-v5-spatial-split.rds", assay = "Spatial")
obj <- sp
obj[["Spatial"]] <- JoinLayers(obj[["Spatial"]])
save(obj, "seurat-v5-spatial-joined.rds", assay = "Spatial")

# 10. Deferred strings by themselves: what R prints for each number is the
#     expected value, computed afresh AFTER the unprinted vector was saved.
ints <- c(1:20, NA, -5L, 2147483647L, 0L)
dbls <- c(0, 1, 2.5, 0.8, 1e5, 123456, 1e15, 1e-4, 0.001, 1/3, 2^-22, -1.5, NA, NaN, Inf, -Inf,
          100000, 1200000, 0.1 + 0.2, 1e22, 1e-290, 1e300, 1e-250, .Machine$double.xmax,
          123456789012345678, 0.1, 1e-5, 12345.6789, 99999.5, 999999.5, 1e5 + 0.1, 2^53, -1e5,
          1e-15, 1.5e-10, 100, 1e3, 1e4, 3e5, 0.30000000000000004, 999999999999999.9,
          9999999999999999, 0.99999999999999989, 1e100, 1.5e100, 1e-100, 123456.7)
rnd <- c(runif(500), runif(500) * 1e6, exp(rnorm(500, 0, 20)), round(runif(500) * 1000, 2),
         as.numeric(sample.int(1e6, 500)), 2^(-60:60), -exp(rnorm(200, 0, 5)),
         seq(0, 5, by = 0.1), (0:100) / 8)
# The hard ones: doubles within four steps of a power of ten (the logarithm of a
# number a hair under 1e-07 is -7 exactly), and numbers whose sixteenth digit is
# at a tie. From 1e-294 up; below that R's own digits hang on the last bit of a
# power of ten near the bottom of the range, and the lab does not follow them.
step <- function(v, k) {
  b <- readBin(writeBin(v, raw(), endian = "little"), "integer", 2, endian = "little")
  lo <- b[1] + k
  if (is.na(lo) || abs(lo) > 2147483000) return(v)
  readBin(writeBin(c(as.integer(lo), b[2]), raw(), endian = "little"), "double", 1, endian = "little")
}
ks <- -290:290
pow <- c(suppressWarnings(as.numeric(paste0("1e", ks))), 10^ks)
near <- unlist(lapply(pow, function(v) vapply(-4:4, function(k) step(v, k), 0)))
ties <- as.numeric(paste0(
  sample(1:9, 3000, TRUE), ".",
  vapply(1:3000, function(i) paste(sample(0:9, 13, TRUE), collapse = ""), ""),
  sample(c("5", "49", "51", "50"), 3000, TRUE), "e", sample(-290:290, 3000, TRUE)))
hard <- c(near, ties)
hard <- hard[is.finite(hard) & hard > 0]

saveRDS(list(int = as.character(ints), dbl = as.character(dbls), rnd = as.character(rnd),
             hard = as.character(hard)),
        file.path(out, "r-deferred-strings.rds"))
expected[["r-deferred-strings.rds"]] <- list(
  int = I(as.character(ints)), dbl = I(as.character(dbls)), rnd = I(as.character(rnd)),
  hard = I(as.character(hard)))

jsonlite::write_json(expected, file.path(out, "seurat-expected.json"),
                     auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
cat("wrote", length(expected) - 1, "objects and seurat-expected.json to", out, "\n")
