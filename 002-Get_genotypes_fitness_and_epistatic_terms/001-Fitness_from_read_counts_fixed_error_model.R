### tRNA Phylogeny Dataset: Use hierarchical fixed effect model to get fitness values and error estimate for each genotype ###

# JDE, September 2017

library(plotrix)
library(ggplot2)
library(stringr)
library(gridExtra)


### Load replicates read counts per genotype after read quality filtering ###
setwd("002-Get_genotypes_fitness_and_epistatic_terms/")
file_path <- "../001-Data/"
filename = paste0(file_path, "tRNA_Phylo_varcounts.tab")

# P is the dataset that will contain the information of each genotype and the number of read counts per replicates in inputs and outputs
P = read.table(filename, header = T) 

# Convert DNA into RNA nomenclature
P$seq = gsub("T", "U", P$seq)
P$ntd_vars = gsub("T", "U", P$ntd_vars)
P$ntd_wt = gsub("T", "U", P$ntd_wt)

### Necessary variables ###
oligo_vars_P_dna = c("G1A", "T2C", "T2G", "G6A", "G6T", "C27T", "A43C", "C46T", "C66T", "C66A", "A69G", "A70G", "A70T", "C71T")
oligo_vars_P = gsub("T", "U", oligo_vars_P_dna)

### Functions  ###
## Merge ntd_wt, pos_vars and ntd_vars to get an unique identifier for each mutant
GetMutantIDs <- function(df) {
  new_id = apply(X = df, MARGIN = 1, FUN = function(x){
    paste(paste0(unlist(strsplit(as.character(x[1]), ",")), unlist(strsplit(as.character(x[2]), ",")), unlist(strsplit(as.character(x[3]), ","))), collapse = ",")
  })
  return(new_id)
}

## Get mutants only with the designed variation
GetIndexById <- function(df, id_idx, oligo_vars) {
  new_idx = apply(X = df, MARGIN = 1, FUN = function(x){
    all(unlist(strsplit(as.character(x[id_idx]), ",")) %in% oligo_vars)
  })
}


### Filter the dataset ###

## Create unique mutation identifier
P$id = GetMutantIDs(P[, c("ntd_wt", "pos_vars", "ntd_vars")])
idx_oligo_vars = GetIndexById(P, grepl("id", colnames(P)), oligo_vars = oligo_vars_P )
idx_zeros = apply(P[, grepl("OUT", colnames(P))], 1, function(x) all(x != 0))


## Keep the WT
WT = P[P$num_vars == 0,]

## Filter based on minimum number of input reads
idx_reads = P$P_IN_1 > 8 & P$P_IN_2 > 8

# Final P dataset
P = rbind(WT, P[idx_oligo_vars & idx_zeros & idx_reads,])


### Calculate fitness and error of it using stratified replicates, fixed error model and error propagetion of this ###

## Calculate frequency of output reads per variant over the total number of reads
OD_P_out = c(1.39, 1.21, 1.30, 1.17, 1.23, 1.14)  
total_reads_out = apply(P[, grepl("P_OUT", colnames(P))], 2, sum)

P_out_freq = do.call("rbind", lapply(1:nrow(P), function(x){
  (P[x, grepl("OUT", colnames(P))]/total_reads_out) * OD_P_out
}))

## Calculate an error for each output using a poisson assumption and an average error using a fixed model
P_out_se = do.call("rbind", lapply(1:nrow(P), function(x){
  sqrt( (1/P[x, grepl("OUT", colnames(P))]) + (1/total_reads_out))
}))

P_out_se$P_OUT_1 = sqrt(1/(P_out_se$P_OUT_11^(-2) + P_out_se$P_OUT_12^(-2) + P_out_se$P_OUT_13^(-2)))
P_out_se$P_OUT_2 = sqrt(1/(P_out_se$P_OUT_21^(-2) + P_out_se$P_OUT_22^(-2) + P_out_se$P_OUT_23^(-2)))

# Weigth the frequency of output based on the error of each replicate
P_out_freq$P_OUT_1 = (P_out_freq$P_OUT_11*(1/(P_out_se$P_OUT_11^2)) + P_out_freq$P_OUT_12*(1/(P_out_se$P_OUT_12^2)) + P_out_freq$P_OUT_13*(1/(P_out_se$P_OUT_13^2)))/
  ((1/(P_out_se$P_OUT_11^2)) + (1/(P_out_se$P_OUT_12^2)) + (1/(P_out_se$P_OUT_13^2)))
P_out_freq$P_OUT_2 = (P_out_freq$P_OUT_21*(1/(P_out_se$P_OUT_21^2)) + P_out_freq$P_OUT_22*(1/(P_out_se$P_OUT_22^2)) + P_out_freq$P_OUT_23*(1/(P_out_se$P_OUT_23^2)))/
  ((1/(P_out_se$P_OUT_21^2)) + (1/(P_out_se$P_OUT_22^2)) + (1/(P_out_se$P_OUT_23^2)))

## Calculate the frequency of each variant in the input and its error
OD_P_in = c(0.015, 0.015, 0.015, 0.015, 0.015, 0.015)  
total_reads_in = apply(P[, grepl("P_IN", colnames(P))], 2, sum)

P_in_freq = do.call("rbind", lapply(1:nrow(P), function(x){
  (P[x, grepl("IN", colnames(P))]/total_reads_in) * OD_P_in
}))

P_in_se = do.call("rbind", lapply(1:nrow(P), function(x){
  sqrt( (1/P[x, grepl("IN", colnames(P))]) + (1/total_reads_in)) 
}))

## Calculate fitness values for each of the input replicates
P_fit = data.frame(P_1 = log2(P_out_freq$P_OUT_1/P_in_freq$P_IN_1), P_2 = log2(P_out_freq$P_OUT_2/P_in_freq$P_IN_2))
P_se = data.frame(P_1 = (1/log(2))*sqrt( (P_out_se$P_OUT_1)^2 + (P_in_se$P_IN_1)^2), 
                  P_2 = (1/log(2))*sqrt( (P_out_se$P_OUT_2)^2 + (P_in_se$P_IN_2)^2))

P_fig_reps = cbind(data.frame(id = P$id), P_fit)

P_fit$fitness_var = (P_fit$P_1*(1/(P_se$P_1^2)) + P_fit$P_2*(1/(P_se$P_2^2)))/((1/(P_se$P_1^2))+(1/(P_se$P_2^2)))
P_se$se_var = sqrt(1/(P_se$P_1^(-2)+ P_se$P_2^(-2)))

P$fitness = log(P_fit$fitness_var/P_fit$fitness_var[P$num_vars == 0])
P$SE = sqrt((P_se$se_var/P_fit$fitness_var)^2 + (P_se$se_var[P$num_vars == 0]/P_fit$fitness_var[P$num_vars==0])^2)


### Additional features ###

# Add long id (nts states)
pos_mutations = sort(as.numeric(as.character(unique(P$pos_vars[P$num_vars == 1]))))
P$long_id = do.call("c", lapply(P$seq, function(x){
  paste0(str_sub(as.character(x),start = pos_mutations, end = pos_mutations), collapse = "")
}))


# Add long id (nts states) using 0, 1 or 2
nt2num = t(as.matrix(do.call("rbind", lapply(pos_mutations, function(x){
  vars_covered = c(unique(P$ntd_wt[P$num_vars == 1 & P$pos_vars == x]), unique(P$ntd_vars[P$num_vars == 1 & P$pos_vars == x]))
  vars_notcovered = c("G", "A", "U", "C")[!c("G", "A", "U", "C") %in% vars_covered]
  vars = c(vars_covered, vars_notcovered)
  df = setNames(data.frame( 0,  1,  2,  3), vars)
  df[, c("G", "A", "U", "C")]
}))))
colnames(nt2num) <- pos_mutations


long_id_M = do.call("rbind",strsplit(P$long_id, ""))
colnames(long_id_M) <-  pos_mutations

num_id_M = do.call("cbind", lapply(1:length(pos_mutations), function(x){
  nt2num[long_id_M[,x],x]
}))
colnames(num_id_M) <- pos_mutations
rownames(num_id_M) <- c()

P$var = apply(num_id_M, 1, function(x) {
  paste(x, collapse ="")
})


### Save dataframe for analysis ###
save(list = c("P"), file = "../001-Data/001-GenotypesFitness.RData")


### Figures ###
P_fig_reps$P_1_rel2wt = P_fig_reps$P_1/P_fig_reps$P_1[P_fig_reps$id == "-0-"]
P_fig_reps$P_2_rel2wt = P_fig_reps$P_2/P_fig_reps$P_2[P_fig_reps$id == "-0-"]
P_fig_reps$IN_1 = P$P_IN_1
P_fig_reps$IN_2 = P$P_IN_2

P_fig_reps$P_1_rel2wt_log = log(P_fig_reps$P_1_rel2wt)
P_fig_reps$P_2_rel2wt_log = log(P_fig_reps$P_2_rel2wt)

## Figure 1d: Correlation between input replicates fitness values.
r = round(cor(P_fig_reps$P_1_rel2wt_log, P_fig_reps$P_2_rel2wt_log, method = "spearman"),2)
p1d <- ggplot(P_fig_reps, aes(x = P_1_rel2wt_log, y = P_2_rel2wt_log)) + theme_classic() +
  stat_bin_hex(bins = 30) + scale_fill_gradient("Genotypes counts", low = "grey85", high = "grey15") +
  xlab("Fitness (Input replicate 1)") + ylab("Fitness (Input replicate 2)") +
  theme(legend.position ="") + annotate("text", x=-0.75, y=-0.1, label=paste0("rs = ",r)) 
p1d


## Figure: Variation of fitness depending on number of input reads. 
pXa <- ggplot(P_fig_reps, aes(x=log10(IN_1), y=P_1_rel2wt_log)) + geom_point(alpha=0.5) + theme_classic() +
  geom_vline(xintercept = log10(9), color="red") +
  xlab("log10(input reads)") + ggtitle("Input rep 1") + ylab("fitness") 
pXb <- ggplot(P_fig_reps, aes(x=log10(IN_2), y=P_2_rel2wt_log)) + geom_point(alpha=0.5) + theme_classic() +
  geom_vline(xintercept = log10(9), color="red") +
  xlab("log10(input reads)") + ggtitle("Input rep 2") +ylab("fitness")
pX <- grid.arrange(pXa,pXb, nrow=2)




