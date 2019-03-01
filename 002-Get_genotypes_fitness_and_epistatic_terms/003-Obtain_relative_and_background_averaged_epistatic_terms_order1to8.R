### Obtain all epistatic terms from order 1 to 8, background relative and background averaged terms ###
# (A) Summary statistics of background relative terms (Extended Data Figure 1d)
# (B) Summary statistics of background average terms (Extended Data Figure 7b)

# JDE, October 2017

library(grDevices)
library(ggplot2)
library(reshape)
library(reshape2)
library(lemon)
library(gridExtra)
library(stringr)
library(plyr)


# Load data
setwd("002-Get_genotypes_fitness_and_epistatic_terms/")
load(file = "../001-Data/001-GenotypesFitness.RData")
P$id[P$class_var=="WT"] = ""

# Initialize vars
colors = data.frame(green = "#33A02C", yellow="gold", orange = "#FF7F00", red = "#E31A1C",  purple= "#992ca0", blue = "#1F78B4", brown="#b45b1f", grey="#6f6f6f")


# Functions
ttest <- function(av, se, df = 5, mu = 0) {
  tstat <- (av-mu)/se
  # Using T-dist
  pval = 2*pt(abs(tstat), df, lower=FALSE)
  return(pval)
}

## Generate all

# Get conversion from variant id (binary) to fitness and SE
P_var2fit = P[,c("fitness", "SE")]
rownames(P_var2fit) <- P$var

L_CD_fit = lapply(1:8, function(x){
  comb = rep(0:x, choose(x, 0:x))
  fit_names = paste("fit_", comb, "_", 1:length(comb), sep="")
  se_names = paste("se_", comb, "_", 1:length(comb), sep="")
  load(file = paste0("~/Google Drive/PhD/Projects/tRNAs/tR-CCU-J/Sequencing_Results/HeatSalt/Data/04-ModifiedDatasets/03-Phylo/CompleteLandscapes/CompleteLands_n", x, ".RData"))
  
  fit = setNames(data.frame(t(apply(CL[,4:ncol(CL)], 1, function(x){
    P_var2fit[x,"fitness"]
  }))),fit_names)
  
  se = setNames(data.frame(t(apply(CL[,4:ncol(CL)], 1, function(x){
    P_var2fit[x,"SE"]
  }))),se_names)
  
  CL = cbind(CL, fit, se)
  CL
})

# Get Ids from vars 
P_var2id = data.frame(id=P$id)
rownames(P_var2id) <- P$var


# Calculate all individual terms of epistasis
L_CD_ep = lapply(1:8, function(x){
  orders = x:0
  sums_v = rep(T, length(orders))
  if (orders[1] %% 2 == 0) {
    sums_v[!(orders) %% 2 == 0] = F
  } else{
    sums_v[(orders) %% 2 == 0] = F
  }
  
  temp = do.call("cbind", lapply(orders, function(y){
    order_regex = paste0("fit_", y, "_")
    if (y == x | y == 0) {
      L_CD_fit[[x]][, grepl(order_regex, colnames(L_CD_fit[[x]]))]
    } else {
      rowSums(L_CD_fit[[x]][, grepl(order_regex, colnames(L_CD_fit[[x]]))])
    }
  }))
  
  if (sum(sums_v) > 1) {
    ep_term = rowSums(temp[,sums_v])
  } else {
    ep_term = temp[,sums_v]
  }
  
  if (sum(!sums_v) > 1) {
    ep_term = ep_term - rowSums(temp[,!sums_v])
  } else {
    ep_term = ep_term - temp[,!sums_v]
  }
  L_CD_fit[[x]]$ep_term = ep_term
  L_CD_fit[[x]]$ep_term_SE = sqrt(rowSums(L_CD_fit[[x]][, grepl("^se_", colnames(L_CD_fit[[x]]))]^2))
  L_CD_fit[[x]]$pval = ttest(av = L_CD_fit[[x]]$ep_term, se = L_CD_fit[[x]]$ep_term_SE, mu=0, df=5)
  L_CD_fit[[x]]$qval_order = p.adjust(L_CD_fit[[x]]$pval, method = "fdr")
  L_CD_fit[[x]]$id_bckg = P_var2id[as.character(L_CD_fit[[x]][,grepl("^id_0_", colnames(L_CD_fit[[x]]))]), "id"]
  L_CD_fit[[x]]$id_mut = P_var2id[as.character(L_CD_fit[[x]][,grepl(paste0("^id_", x, "_"), colnames(L_CD_fit[[x]]))]), "id"]
  L_CD_fit[[x]]$Mutant = do.call("c", lapply(1:nrow(L_CD_fit[[x]]), function(y){
    paste(unlist(strsplit(as.character(L_CD_fit[[x]]$id_mut[y]), ","))[!(unlist(strsplit(as.character(L_CD_fit[[x]]$id_mut[y]), ",")) %in% unlist(strsplit(as.character(L_CD_fit[[x]]$id_bckg[y]),",")))], collapse = ",")
  }))
  L_CD_fit[[x]]
})

## Join in a single DS all mutation combinations with their epistatic terms
CD_ep = do.call("rbind", lapply(1:8, function(x){
  cbind(data.frame(n = rep(x, nrow(L_CD_ep[[x]]))), L_CD_ep[[x]][,!grepl("(^id_[1-9])|(^fit_[1-9])|(^se_)", colnames(L_CD_ep[[x]]))] )
}))
CD_ep$qval_all = p.adjust(CD_ep$pval, method = "fdr")

## (A) Estended Data Figure 1d: Summary complete landscapes 
CD_summary = do.call("rbind", lapply(1:8, function(x){
  n_bckg = count(CD_ep[CD_ep$n == x, c("Mutant")])$freq
  data.frame(n = x, n_comb=length(n_bckg), n_landscapes = sum(CD_ep$n == x), min_bckg = min(n_bckg), median_bckg = median(n_bckg), max_bckg = max(n_bckg))
}))
CD_summary


## General tendecies of epistasis at any order
L_Mutants = setNames(lapply(1:8, function(x){
  unique(L_CD_ep[[x]]$Mutant)
}), 1:8)


EpGlobal_FDRall = do.call("rbind", lapply(1:8, function(x){
  temp = CD_ep[CD_ep$n == x,]
  DS = do.call("rbind", lapply(L_Mutants[[x]], function(y){
    n_bckg = sum(temp$Mutant == y)
    mean_ep_term = mean(temp$ep_term[temp$Mutant == y])
    global_SE = (1/n_bckg)*sqrt(sum(temp$ep_term_SE[temp$Mutant == y]^2))
    positive_fdr10 = sum(temp$ep_term[temp$Mutant == y] > 0 & temp$qval_all[temp$Mutant == y] < 0.1)
    negative_fdr10 = sum(temp$ep_term[temp$Mutant == y] < 0 & temp$qval_all[temp$Mutant == y] < 0.1)
    data.frame(Mutant = y, mean_ep_term = mean_ep_term, n_bckg = n_bckg,
               positive_fdr10 = positive_fdr10, 
               negative_fdr10 = negative_fdr10,
               neutral_fdr10 = n_bckg - (positive_fdr10 + negative_fdr10),
               global_SE = global_SE, 
               global_pval = ttest(av = mean_ep_term, se = global_SE, df = 5))
  }))
  DS$global_qval_order = p.adjust(DS$global_pval, method = "fdr")
  cbind(data.frame(n = rep(x, nrow(DS))), DS)
}))

EpGlobal_FDRall$global_qval_all = p.adjust(EpGlobal_FDRall$global_pval, method = "fdr")


EpGlobal_FDRall$global_state = do.call("c", lapply(1:nrow(EpGlobal_FDRall), function(x){
  if (EpGlobal_FDRall$global_qval_all[x] >= 0.1) {
    "Neutral"
  } else if (EpGlobal_FDRall$mean_ep_term[x] < 0) {
    "Negative"
  } else {
    "Positive"
  }
}))
EpGlobal_FDRall$n = factor(EpGlobal_FDRall$n, levels = 1:8)

SummaryEpGlobal_FDRall = do.call("rbind", lapply(1:8, function(x){
  data.frame(n = x, sig_global_ep = sum(EpGlobal_FDRall$global_qval_all[EpGlobal_FDRall$n == x] < 0.1), n_comb = nrow(EpGlobal_FDRall[EpGlobal_FDRall$n == x,]))
}))

