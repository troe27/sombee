library(hdp)
library(ggplot2)
library(RColorBrewer)
library(lsa)
library(lattice)

SBS52_COLUMNS <- c(
  "C>A,A-C","C>A,A-T","C>A,C-A","C>A,C-C","C>A,C-G","C>A,C-T","C>A,G-C","C>A,T-A","C>A,T-C","C>A,T-T",
  "C>G,A-C","C>G,A-T","C>G,C-A","C>G,C-C","C>G,C-G","C>G,C-T","C>G,G-C","C>G,T-A","C>G,T-C","C>G,T-T",
  "C>T,A-A","C>T,A-C","C>T,A-G","C>T,A-T","C>T,C-A","C>T,C-C","C>T,C-G","C>T,C-T","C>T,G-A","C>T,G-C","C>T,G-G","C>T,G-T","C>T,T-A","C>T,T-C","C>T,T-G","C>T,T-T",
  "T>A,A-C","T>A,A-T","T>A,C-A","T>A,C-C","T>A,C-G","T>A,C-T","T>A,G-C","T>A,T-A","T>A,T-C","T>A,T-T",
  "T>G,A-C","T>G,C-A","T>G,C-C","T>G,C-T","T>G,T-C","T>G,T-T"
)

SBS52_CATEGORIES <- c(
  "A[C>A]C","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]C","T[C>A]A","T[C>A]C","T[C>A]T",
  "A[C>G]C","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]C","T[C>G]A","T[C>G]C","T[C>G]T",
  "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
  "A[T>A]C","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]C","T[T>A]A","T[T>A]C","T[T>A]T",
  "A[T>G]C","C[T>G]A","C[T>G]C","C[T>G]T","T[T>G]C","T[T>G]T"
)

args <- commandArgs(TRUE)
chlist_file <- args[1]
hdp_input_path <- args[2]
output_path <- args[3]
prefix <- args[4]

chlist <- vector("list", 10)
for (i in 1:10) {
  chlist[[i]] <- readRDS(paste0(chlist_file, i, ".rds"))
}

mut_example_multi <- hdp_multi_chain(chlist)
mut_example_multi <- hdp_extract_components(mut_example_multi)

hdp_exposures <- mut_example_multi@comp_dp_distn[["mean"]][
  2:nrow(mut_example_multi@comp_dp_distn[["mean"]]),
  ,
  drop = FALSE
]

input_for_hdp <- read.csv(
  hdp_input_path,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

input_for_hdp[] <- lapply(input_for_hdp, as.numeric)

missing_cols <- setdiff(SBS52_COLUMNS, colnames(input_for_hdp))
extra_cols <- setdiff(colnames(input_for_hdp), SBS52_COLUMNS)

if (length(missing_cols) > 0 || length(extra_cols) > 0) {
  stop(
    paste0(
      "HDP input columns do not match the expected SBS52 schema.\n",
      "Missing: ", paste(missing_cols, collapse = ", "), "\n",
      "Unexpected: ", paste(extra_cols, collapse = ", ")
    )
  )
}

input_for_hdp <- input_for_hdp[, SBS52_COLUMNS, drop = FALSE]
input_for_hdp_sum <- rowSums(input_for_hdp)
input_for_hdp <- input_for_hdp[input_for_hdp_sum > 0, , drop = FALSE]

rownames(hdp_exposures)[
  (nrow(hdp_exposures) - nrow(input_for_hdp) + 1):nrow(hdp_exposures)
] <- rownames(input_for_hdp)

write.csv(
  hdp_exposures,
  paste0(output_path, "/", prefix, "_HDP_exposure.csv"),
  quote = FALSE
)

hdp_sigs <- data.frame(
  t(mut_example_multi@comp_categ_distn[["mean"]]),
  check.names = FALSE
)

if (nrow(hdp_sigs) == length(SBS52_CATEGORIES)) {
  rownames(hdp_sigs) <- SBS52_CATEGORIES
} else {
  stop(
    paste0(
      "Expected ", length(SBS52_CATEGORIES), " mutation categories, but hdp_sigs has ",
      nrow(hdp_sigs), " rows. Check comp_categ_distn dimensions."
    )
  )
}

write.csv(
  hdp_sigs,
  paste0(output_path, "/", prefix, "_HDP_sigs.csv"),
  quote = FALSE
)
