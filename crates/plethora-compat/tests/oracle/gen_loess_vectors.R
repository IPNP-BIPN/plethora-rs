# Golden vectors pinning `plethora-compat::loess` to the local R.
#
# `loess()` exposes its internal k-d tree as `fit$kd`, which lets the port be
# checked in pieces instead of only end to end: the tree first, then the vertex
# fits, then the interpolation. A whole-pipeline check that fails tells you
# nothing about where; these tell you exactly where.
#
# Run:
#     Rscript crates/plethora-compat/tests/oracle/gen_loess_vectors.R \
#         > crates/plethora-compat/tests/data/loess_vectors.tsv
#
# One record per dataset, tab-separated, with comma-separated vectors inside a
# field. %.17g round-trips an IEEE double exactly.

options(stringsAsFactors = FALSE)

g17 <- function(v) paste(vapply(v, function(z) sprintf("%.17g", z), character(1)), collapse = ",")
i_ <- function(v) paste(as.integer(v), collapse = ",")

# Datasets chosen for the shapes that stress the splitter rather than for
# realism: repeated x values decide which side a tie falls on, an odd point
# count decides the median, and a tiny sample must stay a single leaf.
datasets <- list()

set.seed(1)
x <- seq(0.20, 0.72, by = 0.01)
datasets[["gc_model_shape"]] <- list(
    x = x,
    y = log(30) - 8 * (x - 0.42)^2 + rnorm(length(x), 0, 0.02)
)

set.seed(2)
x <- seq(0.20, 0.73, by = 0.01)
datasets[["gc_model_even"]] <- list(
    x = x,
    y = log(25) - 6 * (x - 0.45)^2 + rnorm(length(x), 0, 0.05)
)

# Heavy ties: every x value appears three times.
set.seed(3)
x <- rep(seq(0.25, 0.70, by = 0.05), each = 3)
datasets[["repeated_x"]] <- list(x = x, y = 2 + 3 * x + rnorm(length(x), 0, 0.1))

# Small enough that the tree is a single leaf and predict() is one local fit.
datasets[["single_leaf"]] <- list(x = c(0.1, 0.2, 0.3, 0.4, 0.5), y = c(1, 2, 1.5, 3, 2.5))

# Unevenly spaced, so cell widths differ and the interpolation parameter varies.
set.seed(4)
x <- sort(c(runif(30, 0, 0.3), runif(30, 0.7, 1)))
datasets[["clustered"]] <- list(x = x, y = sin(6 * x) + rnorm(60, 0, 0.05))

# A long run, to push the tree past a couple of levels.
set.seed(5)
x <- sort(runif(200, 0.2, 0.8))
datasets[["large"]] <- list(x = x, y = exp(-((x - 0.5)^2) / 0.02) + rnorm(200, 0, 0.02))

cat("# name\tn\tnc\tnv\tvert\ta\txi\tvval\tx\ty\tfitted\n")

for (nm in names(datasets)) {
    d <- datasets[[nm]]
    fit <- loess(d$y ~ d$x)
    kd <- fit$kd
    p <- kd$parameter

    cat(paste(
        nm,
        length(d$x),
        p[["nc"]],
        p[["nv"]],
        g17(kd$vert),
        i_(kd$a[seq_len(p[["nc"]])]),
        g17(kd$xi[seq_len(p[["nc"]])]),
        g17(kd$vval[seq_len((p[["d"]] + 1) * p[["nv"]])]),
        g17(d$x),
        g17(d$y),
        g17(predict(fit)),
        sep = "\t"
    ), "\n", sep = "")
}
