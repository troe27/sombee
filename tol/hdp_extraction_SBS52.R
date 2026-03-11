library(hdp)
library(ggplot2)
library(RColorBrewer)
library(lsa)
library(lattice)

# chlist_file<-commandArgs(T)[1]
# hdp_input_path<-commandArgs(T)[2]
# output_path<-commandArgs(T)[3]
# prefix<-commandArgs(T)[4]

args <- commandArgs(TRUE)
chlist_file <- args[1]
hdp_input_path <- args[2]
output_path <- args[3]
prefix <- args[4]

#------------extract HDP results---------
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

# old:
# input_for_hdp = read.table(hdp_input_path,check.names = F,header=T)
# input_for_hdp_sum <- apply(input_for_hdp,1,sum)
# input_for_hdp <- input_for_hdp[apply(input_for_hdp,1,sum)>0,]

input_for_hdp <- read.csv(
  hdp_input_path,
  check.names = FALSE,
  row.names = 1
)

input_for_hdp[] <- lapply(input_for_hdp, as.numeric)

# old:
# input_for_hdp_sum <- apply(input_for_hdp,1,sum)
input_for_hdp_sum <- rowSums(input_for_hdp)

# old:
# input_for_hdp <- input_for_hdp[apply(input_for_hdp,1,sum)>0,]
input_for_hdp <- input_for_hdp[input_for_hdp_sum > 0, , drop = FALSE]

rownames(hdp_exposures)[
  (nrow(hdp_exposures) - nrow(input_for_hdp) + 1):nrow(hdp_exposures)
] <- rownames(input_for_hdp)

write.csv(
  hdp_exposures,
  paste0(output_path, "/", prefix, "_HDP_exposure.csv"),
  quote = FALSE
)

# old:
# hdp_sigs=data.frame(t(mut_example_multi@comp_categ_distn[["mean"]][1:dim(mut_example_multi@comp_categ_distn[["mean"]])[1],]))
# categories = paste0(substr(colnames(input_for_hdp),5,5),'[',substr(colnames(input_for_hdp),1,3),']',substr(colnames(input_for_hdp),7,7))
# rownames(hdp_sigs)<-colnames(categories)

hdp_sigs <- data.frame(
  t(mut_example_multi@comp_categ_distn[["mean"]]),
  check.names = FALSE
)

categories <- c(
  "A[C>A]C","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]C","T[C>A]A","T[C>A]C","T[C>A]T",
  "A[C>G]C","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]C","T[C>G]A","T[C>G]C","T[C>G]T",
  "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
  "A[T>A]C","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]C","T[T>A]A","T[T>A]C","T[T>A]T",
  "A[T>G]C","C[T>G]A","C[T>G]C","C[T>G]T","T[T>G]C","T[T>G]T"
)

print(dim(mut_example_multi@comp_categ_distn[["mean"]]))
print(dim(hdp_sigs))

if (nrow(hdp_sigs) == length(categories)) {
  rownames(hdp_sigs) <- categories
} else {
  stop(
    paste0(
      "Expected ", length(categories), " mutation categories, but hdp_sigs has ",
      nrow(hdp_sigs), " rows. Check comp_categ_distn dimensions."
    )
  )
}

write.csv(
  hdp_sigs,
  paste0(output_path, "/", prefix, "_HDP_sigs.csv"),
  quote = FALSE
)
