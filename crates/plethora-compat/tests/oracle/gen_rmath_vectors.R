# Golden vectors pinning `plethora-compat::rmath` to the local R.
#
# Two behaviours of gc_correction.R ride on this file.
#
#   1. `mutate(X, percent.gc = round(percent.gc, 2))` decides which GC bin every
#      domain lands in, and so which correction factor it gets. The GC file is
#      full of three-decimal values, so exact halves such as 0.345 are common,
#      and they are precisely where a naive round disagrees with R's.
#
#   2. `write.table(X, sep = "\t", row.names = FALSE, quote = FALSE)` decides how
#      the resulting doubles are spelled in the output file. Matching upstream
#      byte for byte means matching that spelling, not just the value.
#
# Run:
#     Rscript crates/plethora-compat/tests/oracle/gen_rmath_vectors.R \
#         > crates/plethora-compat/tests/data/rmath_vectors.tsv

options(stringsAsFactors = FALSE, scipen = 0)

# %.17g round-trips an IEEE double exactly, so the Rust side can compare bit
# patterns instead of tolerances.
g17 <- function(v) vapply(v, function(z) sprintf("%.17g", z), character(1))

values <- c(
    # Every three-decimal value, which is the exact domain of the GC file for
    # the 1000 bp baseline domains: percent GC is a count over 1000.
    seq(0, 1, by = 0.001),
    # Halfway cases at two decimals, the ones round-half-to-even decides.
    c(0.005, 0.015, 0.025, 0.045, 0.055, 0.125, 0.135, 0.145, 0.155),
    c(0.345, 0.355, 0.365, 0.375, 0.385, 0.395),
    c(2.5, 1.5, 0.5, -0.5, -1.5, -2.5),
    # Values whose binary representation sits just under or just over a half.
    c(0.1 + 0.2, 1 / 3, 2 / 3, 1e-8, 1 - 1e-8),
    # Negative, large and tiny, to exercise the sign and magnitude branches.
    c(-0.345, -0.355, -1234.5678, 1234.5678),
    c(1e15, 1e16, 1e-300, 1e300, 0, -0),
    c(123456789.123456789, 0.049999999999999996)
)
values <- unique(values)

cat("# section\tinput\trounded2\tas_written\n")

for (x in values) {
    cat(paste("round", g17(x), g17(round(x, 2)), sep = "\t"), "\n", sep = "")
}

# How write.table spells a double. Written to a real connection rather than
# reasoned about, because the answer depends on R's internal formatting rules
# rather than on any single documented function.
written <- c(
    0, 1, -1, 0.5, 0.1, 1 / 3, 2 / 3, 1e-5, 1e-4, 1e5, 1e15, 1e16, 1e-15,
    28.2794, 2.55338, 0.66548, 0.345, 1.0, 100, 1e300, 1e-300,
    123456789.123456789, 0.049999999999999996, -0.0001234, 3.14159265358979,
    1234567890123456, 0.30000000000000004
)
tmp <- tempfile()
write.table(data.frame(v = written), tmp, sep = "\t", row.names = FALSE,
            quote = FALSE, col.names = FALSE)
spelled <- readLines(tmp)
unlink(tmp)

for (i in seq_along(written)) {
    cat(paste("write", g17(written[i]), "", spelled[i], sep = "\t"), "\n", sep = "")
}
