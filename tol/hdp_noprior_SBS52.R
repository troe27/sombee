options(stringsAsFactors = F)
mut_file<-commandArgs(T)[1]
iter <-as.numeric(commandArgs(T)[2])
out <-as.character(commandArgs(T)[3])


library(hdp)
print(paste0(out,iter,".rds"))
#input_for_hdp<-read.table(mut_file)

#---------------
x <- read.csv(mut_file, header = FALSE, skip = 1, row.names = 1, check.names = FALSE)
x[] <- lapply(x, as.numeric)

stopifnot(ncol(x) == 52)

colnames(x) <- c(
  "C>A_A-C","C>A_A-T","C>A_C-A","C>A_C-C","C>A_C-G","C>A_C-T","C>A_G-C","C>A_T-A","C>A_T-C","C>A_T-T",
  "C>G_A-C","C>G_A-T","C>G_C-A","C>G_C-C","C>G_C-G","C>G_C-T","C>G_G-C","C>G_T-A","C>G_T-C","C>G_T-T",
  "C>T_A-A","C>T_A-C","C>T_A-G","C>T_A-T","C>T_C-A","C>T_C-C","C>T_C-G","C>T_C-T","C>T_G-A","C>T_G-C","C>T_G-G","C>T_G-T","C>T_T-A","C>T_T-C","C>T_T-G","C>T_T-T",
  "T>A_A-C","T>A_A-T","T>A_C-A","T>A_C-C","T>A_C-G","T>A_C-T","T>A_G-C","T>A_T-A","T>A_T-C","T>A_T-T",
  "T>G_A-C","T>G_C-A","T>G_C-C","T>G_C-T","T>G_T-C","T>G_T-T"
)
input_for_hdp <- x
#---------------

input_for_hdp_sum <- apply(input_for_hdp,1,sum)
input_for_hdp <- input_for_hdp[apply(input_for_hdp,1,sum)>0,]
input_for_hdp_sum <- apply(input_for_hdp,1,sum)

median_muts = 3000
input_for_hdp = round(median_muts*input_for_hdp/input_for_hdp_sum)


ppindex <- c(0, rep(1, nrow(input_for_hdp)))
cpindex <- c(1, rep(2, nrow(input_for_hdp)))


hdp_mut <- hdp_init(ppindex = ppindex, # index of parental node
		    cpindex = cpindex, # index of the CP to use
		    hh = rep(1, 52), # prior is uniform over 52 categories
                    alphaa = rep(1, length(unique(cpindex))), # shape hyperparameters for 2 CPs
                    alphab = rep(5, length(unique(cpindex))))  # rate hyperparameters for 2 CPs

hdp_mut <- hdp_setdata(hdp_mut,  
		       dpindex = 2:numdp(hdp_mut), # index of nodes to add data to
		       input_for_hdp) # input data (mutation counts, sample rows match up with specified dpindex)

hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=iter*200)

chlist <- hdp_posterior(hdp_activated,
                        burnin=20000,
			n=100,
			space=1000,
			cpiter=3,
			seed=iter*1e3)
saveRDS(chlist, paste0(out,iter,".rds"))
print(paste0("HDP chain saved to",out,iter,".rds"))

