options(stringsAsFactors = FALSE)

args <- commandArgs(TRUE)
mut_file <- args[1]
iter <- as.numeric(args[2])
out <- as.character(args[3])

library(hdp)

SBS52_COLUMNS <- c(
  "C>A,A-C","C>A,A-T","C>A,C-A","C>A,C-C","C>A,C-G","C>A,C-T","C>A,G-C","C>A,T-A","C>A,T-C","C>A,T-T",
  "C>G,A-C","C>G,A-T","C>G,C-A","C>G,C-C","C>G,C-G","C>G,C-T","C>G,G-C","C>G,T-A","C>G,T-C","C>G,T-T",
  "C>T,A-A","C>T,A-C","C>T,A-G","C>T,A-T","C>T,C-A","C>T,C-C","C>T,C-G","C>T,C-T","C>T,G-A","C>T,G-C","C>T,G-G","C>T,G-T","C>T,T-A","C>T,T-C","C>T,T-G","C>T,T-T",
  "T>A,A-C","T>A,A-T","T>A,C-A","T>A,C-C","T>A,C-G","T>A,C-T","T>A,G-C","T>A,T-A","T>A,T-C","T>A,T-T",
  "T>G,A-C","T>G,C-A","T>G,C-C","T>G,C-T","T>G,T-C","T>G,T-T"
)

print(paste0(out, iter, ".rds"))

input_for_hdp <- read.csv(
  mut_file,
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
input_for_hdp_sum <- rowSums(input_for_hdp)

median_muts <- 3000
input_for_hdp <- round(median_muts * input_for_hdp / input_for_hdp_sum)

ppindex <- c(0, rep(1, nrow(input_for_hdp)))
cpindex <- c(1, rep(2, nrow(input_for_hdp)))

hdp_mut <- hdp_init(
  ppindex = ppindex,
  cpindex = cpindex,
  hh = rep(1, 52),
  alphaa = rep(1, length(unique(cpindex))),
  alphab = rep(5, length(unique(cpindex)))
)

hdp_mut <- hdp_setdata(
  hdp_mut,
  dpindex = 2:numdp(hdp_mut),
  input_for_hdp
)

hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc = 10, seed = iter * 200)

chlist <- hdp_posterior(
  hdp_activated,
  burnin = 20000,
  n = 100,
  space = 1000,
  cpiter = 3,
  seed = iter * 1e3
)

saveRDS(chlist, paste0(out, iter, ".rds"))
print(paste0("HDP chain saved to ", out, iter, ".rds"))
