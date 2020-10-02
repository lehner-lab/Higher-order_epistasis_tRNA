## HOep - CV approach for introducing significant coefficients of the same order at a time
# Run on the cluster

# JDE, March 2018

library(grDevices)
library(ggplot2)
library(reshape)
library(reshape2)
library(lemon)
library(gridExtra)
library(stringr)
library(plyr)
library(RColorBrewer)
library(parallel)
library(matrixStats)
library(caret)

# Load data
setwd("~/Google Drive/PhD/Projects/tRNAs/tR-CCU-J/Sequencing_Results/HeatSalt/Manuscript/NatureSub2/01-Log_Biallelic/HigherOrderEpistasis/")
load("~/Google Drive/PhD/Projects/tRNAs/tR-CCU-J/Sequencing_Results/HeatSalt/Data/04-ModifiedDatasets/03-Phylo/PhyloDS_NoAC_LOG_FixErrorModel.RData")
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

walsh_matrix = function(order){
  if (order < 1) {
    m = 1
  } else {
    m = 1
    for(i in 1:order){
      m = rbind(  cbind(m,m),
                  cbind(m,m*-1)
      )
    }
  }
  m
}

weighting_matrix = function(order){
  if (order < 1) {
    m = 1
  } else {
    m = 1
    zm = 0
    for(i in 1:order){
      m = rbind(  cbind(0.5*m,zm),
                  cbind(zm,m*-1))
      zm = matrix(0,nrow(m),ncol(m))
    }
  }
  m
}

relative_matrix = function(order){
  if (order < 1) {
    m = 1
  } else {
    m = 1
    for(i in 1:order){
      m = rbind(  cbind(m,0*m),
                  cbind(-1*m,m)
      )
    }
  }
  m
}

GoP = function(obs, pred) {
  Rtot = sum((obs-mean(obs))^2)
  Rres = sum((obs-pred)^2)
  R2 = 1/(1+(Rres/Rtot))
  return(R2)
}

VarExp = function(obs, pred) {
  Rtot = sum((obs-mean(obs))^2)
  Rres = sum((obs-pred)^2)
  R2 = 1-(Rres/Rtot)
  return(R2)
}

RMSE = function(obs, pred, n) {
  Rres = sum((obs-pred)^2)
  R2 = sqrt((Rres/n))
  return(R2)
}

partialAv_matrix = function(genotypeToRemove,order){
  m = do.call("cbind",mclapply(genotypeToRemove,function(w){
    pow = order
    quad_seq = NULL
    w2 = w
    for(i in pow:1){
      tmp = 2^i/2
      if(w2 > 2^i/2){
        quad_seq = c(quad_seq,F)
        w2 = w2 - tmp
      }else{
        quad_seq=c(quad_seq,T)
      }
    }
    m = matrix(0,ncol=1,nrow=1)
    for(i in length(quad_seq):1){
      if(quad_seq[i] == 1){
        m2 = cbind(m,matrix(1,nrow=nrow(m),ncol=ncol(m)))
      }else{
        m2 = cbind(matrix(1,nrow=nrow(m),ncol=ncol(m)),m)
      }
      m = rbind(m2,cbind(m,m))
    }
    as.vector(m)
  }))
  matrix(rowProds(m),nrow=2^order)
}


pred_partialAv = function(order,fit,genotypeToRemove){
  
  m = partialAv_matrix(genotypeToRemove, order, fit)
  tmp = rowSums(m)
  tmp[tmp==0] = 1
  weights = 2^order/tmp
  (weighting_matrix(order)*weights) %*% ((walsh_matrix(order) * m) %*% fit)
}

weighted_partial = function(order, partial_M) {
  l = log2(rowSums(abs(relative_matrix(order))))
  rs = rowSums(partial_M)
  d = ((-1)^l)*(2^l)/rs
  d[d == Inf | d == -Inf] <- 0
  diag(d)
}


### Get Data for all the 8-loci landscapes ###
P_var2fit = P[,c("fitness", "SE")]
rownames(P_var2fit) <- P$var

x=8
comb = rep(0:x, choose(x, 0:x))
fit_names = paste("fit_", comb, "_", 1:length(comb), sep="")
se_names = paste("se_", comb, "_", 1:length(comb), sep="")
load(file = "~/Google Drive/PhD/Projects/tRNAs/tR-CCU-J/Sequencing_Results/HeatSalt/Data/04-ModifiedDatasets/03-Phylo/CompleteLandscapes/CompleteLands_n8.RData")

fit = setNames(data.frame(t(apply(CL[,4:ncol(CL)], 1, function(x){
  P_var2fit[x,"fitness"]
}))),fit_names)

se = setNames(data.frame(t(apply(CL[,4:ncol(CL)], 1, function(x){
  P_var2fit[x,"SE"]
}))),se_names)

CL8 = cbind(CL, fit, se)

IDs = CL8[1, grepl("^id_", colnames(CL8))]
temp = CL8[,grepl("fit_", colnames(CL8))]
M8_Fit = t(data.matrix(temp[,order(IDs)])) # Rows are genotypes, columns number of sublandscapes


temp_se = CL8[,grepl("se_", colnames(CL8))]
M8_SE = t(data.matrix(temp_se[,order(IDs)]))

temp_id = CL8[,grepl("^id_", colnames(CL8))]
M8_ID = t(temp_id[,order(IDs)])

### ###

### Calculations ###


o = 8 # Order
n_subland = ncol(M8_Fit)

# Run the Pred, Obs and RMSE all together after setting the folds

M8_Walsh = walsh_matrix(o)
M8_Relative = relative_matrix(o)
M8_Weighting = weighting_matrix(o)

coeff_orders_loci8 = log2(rowSums(abs(M8_Relative)))

L_CoeffsIter_Folds = lapply(1:n_subland, function(x){
  folds = createFolds(1:256, k=10)
})

# Coefficients and rank for each 10fold CV - Usable whenever iterating by order or by coefficients
L_OrderIter_Coeffs = lapply(1:n_subland, function(x){
  cat("\r",x,"\t\t\t\t")
  mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    Zp = partialAv_matrix(genotypeToRemove = L_CoeffsIter_Folds[[x]][[y]], order = o)
    Wp = weighted_partial(order = o, partial_M = Zp)
    W_Coeffs = Wp %*% (Zp*M8_Walsh) %*% M8_Fit[,x]
    W_SE = abs(Wp) %*% sqrt(abs(Zp*M8_Walsh) %*% M8_SE[,x]^2)    
    W_tstat = W_Coeffs/W_SE
    cbind(W_Coeffs,rank(-abs(W_tstat)))
  }, mc.cores = detectCores(all.tests = FALSE, logical = TRUE))
})


L_CoeffsIter_Pred = lapply(1:n_subland, function(x){
  cat("\r",x,"\t\t\t\t")
  mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    coeffs_M = do.call("cbind", lapply(1:256, function(z){
      coeffs_ranks = L_OrderIter_Coeffs[[x]][[y]]
      coeffs_ranks[coeffs_ranks[,2] > z,1] <- 0
      coeffs_ranks[,1]
    }))
    Zp = partialAv_matrix(genotypeToRemove = L_CoeffsIter_Folds[[x]][[y]], order = o)
    Wp = weighted_partial(order = o, partial_M = Zp)
    diag(Wp) = 1/diag(Wp)
    Wp[!is.finite(Wp)] = 0
    solve(M8_Walsh) %*% (Wp %*% coeffs_M)
  }, mc.cores = detectCores(all.tests = FALSE, logical = TRUE))
})

# L_CoeffsIter_Pred2 = lapply(1:n_subland, function(x){
#   cat("\r",x,"\t\t\t\t")
#   mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
#     coeffs_M = do.call("cbind", lapply(1:256, function(z){
#       coeffs_ranks = L_OrderIter_Coeffs[[x]][[y]]
#       coeffs_ranks[coeffs_ranks[,2] > z,1] <- 0
#       coeffs_ranks[,1]
#     }))
#     Wp = M8_Weighting
#     diag(Wp) = 1/diag(Wp)
#     solve(M8_Walsh) %*% (Wp %*% coeffs_M)
#   }, mc.cores = detectCores(all.tests = FALSE, logical = TRUE))
# })

M_CoeffsIter_OrderAdded = do.call("cbind", lapply(1:n_subland, function(x){
  do.call("cbind", lapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    coeffs_ranks = L_OrderIter_Coeffs[[x]][[y]][,2]
    coeff_orders_loci8[order(coeffs_ranks)]
  }))
}))


L_CoeffsIter_RMSE = lapply(1:n_subland, function(x){
  cat("\r",x,"\t\t\t\t")
  do.call("rbind",mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    k = L_CoeffsIter_Folds[[x]][[y]]
    pred = L_CoeffsIter_Pred[[x]][[y]][k,]
    obs = M8_Fit[k,x]
    test_err = sqrt(colMeans((obs - pred)^2))
    pred = L_CoeffsIter_Pred[[x]][[y]][-k,]
    obs = M8_Fit[-k,x]
    train_err = sqrt(colMeans((obs - pred)^2))
    cbind(test_err,train_err)
  }, mc.cores = detectCores(all.tests = FALSE, logical = TRUE)))
})

# L_CoeffsIter_RMSE2 = lapply(1:n_subland, function(x){
#   cat("\r",x,"\t\t\t\t")
#   do.call("rbind",mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
#     k = L_CoeffsIter_Folds[[x]][[y]]
#     pred = L_CoeffsIter_Pred2[[x]][[y]][k,]
#     obs = M8_Fit[k,x]
#     test_err = sqrt(colMeans((obs - pred)^2))
#     pred = L_CoeffsIter_Pred2[[x]][[y]][-k,]
#     obs = M8_Fit[-k,x]
#     train_err = sqrt(colMeans((obs - pred)^2))
#     cbind(test_err,train_err)
#   }, mc.cores = detectCores(all.tests = FALSE, logical = TRUE)))
# })

M_CoeffsIter_minRMSE = do.call("cbind", lapply(1:n_subland, function(x){
  do.call("c",lapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    k = L_CoeffsIter_Folds[[x]][[y]]
    pred = L_CoeffsIter_Pred[[x]][[y]][k,]
    obs = M8_Fit[k,x]
    test_err = sqrt(colMeans((obs - pred)^2))
    which(test_err == min(test_err))
  }))
}))
# M_CoeffsIter_minRMSE2 = do.call("cbind", lapply(1:n_subland, function(x){
#   do.call("c",lapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
#     k = L_CoeffsIter_Folds[[x]][[y]]
#     pred = L_CoeffsIter_Pred2[[x]][[y]][k,]
#     obs = M8_Fit[k,x]
#     test_err = sqrt(colMeans((obs - pred)^2))
#     which(test_err == min(test_err))
#   }))
# }))
v_CoeffsIter_minRMSE = as.vector(M_CoeffsIter_minRMSE)
# v_CoeffsIter_minRMSE2 = as.vector(M_CoeffsIter_minRMSE2)


p6a <- ggplot(data.frame(ncoeffs = v_CoeffsIter_minRMSE), aes(x=ncoeffs)) + theme_classic() +
  geom_histogram(bins=20, fill="grey80") +
  xlab("Number of coefficients for models that give the lowest RMSE\n(different sublandscapes and 10-fold CV tests)")
p6b <- ggplot(data.frame(ncoeffs = apply(M_CoeffsIter_minRMSE,2, median)), aes(x=ncoeffs)) + theme_classic() +
  geom_histogram(bins=20, fill="grey80") +
  xlab("Median number of coefficients for models that give the lowest RMSE\n(different sublandscapes)")
p6 <- grid.arrange(p6a,p6b, nrow = 2)
# ggsave("Figures-VarianceExplained/CV/WWalsh/20180320_CVWhalshWeighted_DistNumberBestModels.pdf", p6, width = 4, height = 4.5)

# p6a <- ggplot(data.frame(ncoeffs = v_CoeffsIter_minRMSE2), aes(x=ncoeffs)) + theme_classic() +
#   geom_histogram(bins=20, fill="grey80") +
#   xlab("Number of coefficients for models that give the lowest RMSE\n(different sublandscapes and 10-fold CV tests)")
# p6b <- ggplot(data.frame(ncoeffs = apply(M_CoeffsIter_minRMSE2,2, median)), aes(x=ncoeffs)) + theme_classic() +
#   geom_histogram(bins=20, fill="grey80") +
#   xlab("Median number of coefficients for models that give the lowest RMSE\n(different sublandscapes)")
# p6 <- grid.arrange(p6a,p6b, nrow = 2)
# 


M_CoeffsIter_sumOrdersAdded = do.call("cbind", lapply(1:length(v_CoeffsIter_minRMSE), function(y){
  table(factor(M_CoeffsIter_OrderAdded[1:v_CoeffsIter_minRMSE[y],y], levels = c(0:6)))
}))


v_CoeffsIter_minRMSE_medianCV = rep(floor(apply(M_CoeffsIter_minRMSE, 2, median)), each=10)
M_CoeffsIter_sumOrdersAdded_medianCV = do.call("cbind", lapply(1:length(v_CoeffsIter_minRMSE_medianCV), function(y){
  table(factor(M_CoeffsIter_OrderAdded[1:v_CoeffsIter_minRMSE_medianCV[y],y], levels = c(0:6)))
}))
DS_M_CoeffsIter_sumOrdersAdded = data.frame(orders = factor(0:6), count = rowMeans(M_CoeffsIter_sumOrdersAdded),
                                            count_ci = apply(M_CoeffsIter_sumOrdersAdded,1,sd) / sqrt(n_subland) * 1.96,
                                            count_rel = rowMeans(M_CoeffsIter_sumOrdersAdded)/choose(8,0:6),
                                            count_rel_ci = (apply(M_CoeffsIter_sumOrdersAdded,1,sd) / sqrt(n_subland) * 1.96)/choose(8,0:6)                                                     )
cols_orders = colorRampPalette(colors = c("grey90", "grey10"))(7)
p2_a <- ggplot(DS_M_CoeffsIter_sumOrdersAdded, aes(x=orders, y=count)) + theme_classic() +
  geom_col(aes(fill=orders)) + geom_errorbar(aes(ymin = count - count_ci, ymax = count + count_ci), width=0.2) +
  scale_fill_manual("", values = cols_orders) + ylab("absolute count") + theme(legend.position = "")
p2_b <- ggplot(DS_M_CoeffsIter_sumOrdersAdded, aes(x=orders, y=count_rel)) + theme_classic() +
  geom_col(aes(fill=orders)) + geom_errorbar(aes(ymin = count_rel - count_rel_ci, ymax = count_rel + count_rel_ci), width=0.2) +
  scale_fill_manual("", values = cols_orders) + ylab("relative count to num of possible coeffs") + theme(legend.position = "")
p2 <- grid.arrange(p2_a, p2_b, nrow=2) 
# ggsave("Figures-VarianceExplained/CV/WWalsh/20180320_CVWhalshWeighted_OrdersAdded.pdf", p2, width = 3, height = 4.3)

# # Using Walsh weighed not partial walsh weighted
# M_CoeffsIter_sumOrdersAdded2 = do.call("cbind", lapply(1:length(v_CoeffsIter_minRMSE2), function(y){
#   table(factor(M_CoeffsIter_OrderAdded[1:v_CoeffsIter_minRMSE2[y],y], levels = c(0:6)))
# }))
# 
# v_CoeffsIter_minRMSE_medianCV2 = rep(floor(apply(M_CoeffsIter_minRMSE2, 2, median)), each=10)
# M_CoeffsIter_sumOrdersAdded_medianCV2 = do.call("cbind", lapply(1:length(v_CoeffsIter_minRMSE_medianCV2), function(y){
#   table(factor(M_CoeffsIter_OrderAdded[1:v_CoeffsIter_minRMSE_medianCV2[y],y], levels = c(0:6)))
# }))
# DS_M_CoeffsIter_sumOrdersAdded2 = data.frame(orders = factor(0:6), count = rowMeans(M_CoeffsIter_sumOrdersAdded2),
#                                             count_ci = apply(M_CoeffsIter_sumOrdersAdded2,1,sd) / sqrt(n_subland) * 1.96,
#                                             count_rel = rowMeans(M_CoeffsIter_sumOrdersAdded2)/choose(8,0:6),
#                                             count_rel_ci = (apply(M_CoeffsIter_sumOrdersAdded2,1,sd) / sqrt(n_subland) * 1.96)/choose(8,0:6)                                                     )
# cols_orders = colorRampPalette(colors = c("grey90", "grey10"))(7)
# p2_a <- ggplot(DS_M_CoeffsIter_sumOrdersAdded2, aes(x=orders, y=count)) + theme_classic() +
#   geom_col(aes(fill=orders)) + geom_errorbar(aes(ymin = count - count_ci, ymax = count + count_ci), width=0.2) +
#   scale_fill_manual("", values = cols_orders) + ylab("absolute count") + theme(legend.position = "")
# p2_b <- ggplot(DS_M_CoeffsIter_sumOrdersAdded2, aes(x=orders, y=count_rel)) + theme_classic() +
#   geom_col(aes(fill=orders)) + geom_errorbar(aes(ymin = count_rel - count_rel_ci, ymax = count_rel + count_rel_ci), width=0.2) +
#   scale_fill_manual("", values = cols_orders) + ylab("relative count to num of possible coeffs") + theme(legend.position = "")
# p2 <- grid.arrange(p2_a, p2_b, nrow=2) 




# median numbero of coefficients added across folds
DS_M_CoeffsIter_sumOrdersAdded_medianCV = data.frame(orders = factor(0:6), count = rowMeans(M_CoeffsIter_sumOrdersAdded_medianCV),
                                                     count_ci = apply(M_CoeffsIter_sumOrdersAdded_medianCV,1,sd) / sqrt(n_subland) * 1.96,
                                                     count_rel = rowMeans(M_CoeffsIter_sumOrdersAdded_medianCV)/choose(8,0:6),
                                                     count_rel_ci = (apply(M_CoeffsIter_sumOrdersAdded_medianCV,1,sd) / sqrt(n_subland) * 1.96)/choose(8,0:6)                                                     )
cols_orders = colorRampPalette(colors = c("grey90", "grey10"))(7)
p2_medianCV_a <- ggplot(DS_M_CoeffsIter_sumOrdersAdded_medianCV, aes(x=orders, y=count)) + theme_classic() +
  geom_col(aes(fill=orders)) + geom_errorbar(aes(ymin = count - count_ci, ymax = count + count_ci), width=0.2) +
  scale_fill_manual("", values = cols_orders) + ylab("absolute count") + theme(legend.position = "")
p2_medianCV_b <- ggplot(DS_M_CoeffsIter_sumOrdersAdded_medianCV, aes(x=orders, y=count_rel)) + theme_classic() +
  geom_col(aes(fill=orders)) + geom_errorbar(aes(ymin = count_rel - count_rel_ci, ymax = count_rel + count_rel_ci), width=0.2) +
  scale_fill_manual("", values = cols_orders) + ylab("relative count to num of possible coeffs") + theme(legend.position = "")
p2b <- grid.arrange(p2_medianCV_a, p2_medianCV_b, nrow=2) 
# ggsave("Figures-VarianceExplained/CV/WWalsh/20180320_CVWhalshWeighted_OrdersAddedMedianCV.pdf", p2b, width = 6, height = 7)

# Accuracy by RMSE
M_CoeffsIter_RMSE_comb1 = do.call("cbind",lapply(1:n_subland, function(x){
  xx = L_CoeffsIter_RMSE[[x]]
  train_err = matrix(xx[,2],nrow=256)
  test_err = matrix(xx[,1],nrow=256)
  test = sqrt(rowSums(test_err^2) / length(L_CoeffsIter_Folds[[x]]))
  train = sqrt(rowSums(train_err^2) / length(L_CoeffsIter_Folds[[x]]))
  cbind(test,train)
}))

RMSE_test_CoeffsIter = M_CoeffsIter_RMSE_comb1[,rep(c(T,F),n_subland)]
RMSE_train_CoeffsIter = M_CoeffsIter_RMSE_comb1[,rep(c(F,T),n_subland)]
DS_RMSE_CoeffsIter_mean = data.frame(coeffs = c(1:256,1:256), rmse = c(rowMeans(RMSE_test_CoeffsIter), rowMeans(RMSE_train_CoeffsIter)),
                                     rmse_ci=c(apply(RMSE_test_CoeffsIter,1,sd) / sqrt(n_subland) * 1.96, apply(RMSE_train_CoeffsIter,1,sd) / sqrt(n_subland) * 1.96),
                                     dataset = rep(c("test", "train"), each=256),
                                     added_o = factor(rep(apply(M_CoeffsIter_OrderAdded,1, median),2)))
p3 <- ggplot(DS_RMSE_CoeffsIter_mean[DS_RMSE_CoeffsIter_mean$coeffs<=250,], aes(x=coeffs, y=rmse)) + theme_classic() +
  geom_errorbar(aes(ymin = rmse -rmse_ci, ymax = rmse + rmse_ci), width = 0.5) +
  geom_point(size=2, aes(fill=dataset), shape=21) +
  geom_point(data=DS_RMSE_CoeffsIter_mean[DS_RMSE_CoeffsIter_mean$rmse == min(DS_RMSE_CoeffsIter_mean$rmse[DS_RMSE_CoeffsIter_mean$dataset == "test"]),], color="#cc0000", shape=21, fill= "#7F3F98", size=3) +
  theme(legend.position = c(0.8,0.8)) +
  scale_fill_manual("", values = c("#7F3F98", "#FFDE17")) +
  xlab("Number of coefficients added") +
  ylab("RMSE")
p3
# ggsave("Figures-VarianceExplained/CV/WWalsh/20180320_CVWhalshWeighted_IterCoeffsTestTrainMean.pdf",p3,  width = 7, height = 5)
p3b <- ggplot(DS_RMSE_CoeffsIter_mean[DS_RMSE_CoeffsIter_mean$dataset =="test" & DS_RMSE_CoeffsIter_mean$coeffs<=20 ,], aes(x=coeffs, y=rmse)) + theme_classic() +
  #geom_line(data=cbind(melt(RMSE_test_CoeffsIter[1:25,]), data.frame(land = rep(1:76, each=25))), aes(x=X1, y=value,group=land), alpha=0.2) +
  geom_errorbar(aes(ymin = rmse -rmse_ci, ymax = rmse + rmse_ci), width = 0.2) +
  geom_point(size=4, shape=21, aes(fill=added_o)) +
  scale_fill_manual("Median order of\n coefficient added", values = colorRampPalette(colors = c("#f2ebf4", "#7F3F98"))(4)) +
  xlab("Number of coefficients added") + scale_x_continuous(limits = c(0.9,20.1)) +
  ylab("RMSE") + theme(legend.position = c(0.7,0.7) )
p3b
# ggsave("Figures-VarianceExplained/CV/WWalsh/20180320_CVWhalshWeighted_IterCoeffsTestMedianOrderCoeffAdded.pdf", p3b, width = 3, height = 2.9)


# # using walsh weighted instead of parital walsh weighed
# M_CoeffsIter_RMSE_comb1_2 = do.call("cbind",lapply(1:n_subland, function(x){
#   xx = L_CoeffsIter_RMSE2[[x]]
#   train_err = matrix(xx[,2],nrow=256)
#   test_err = matrix(xx[,1],nrow=256)
#   test = sqrt(rowSums(test_err^2) / length(L_CoeffsIter_Folds[[x]]))
#   train = sqrt(rowSums(train_err^2) / length(L_CoeffsIter_Folds[[x]]))
#   cbind(test,train)
# }))
# 
# RMSE_test_CoeffsIter2 = M_CoeffsIter_RMSE_comb1_2[,rep(c(T,F),n_subland)]
# RMSE_train_CoeffsIter2 = M_CoeffsIter_RMSE_comb1_2[,rep(c(F,T),n_subland)]
# DS_RMSE_CoeffsIter_mean2 = data.frame(coeffs = c(1:256,1:256), rmse = c(rowMeans(RMSE_test_CoeffsIter2), rowMeans(RMSE_train_CoeffsIter2)),
#                                      rmse_ci=c(apply(RMSE_test_CoeffsIter2,1,sd) / sqrt(n_subland) * 1.96, apply(RMSE_train_CoeffsIter2,1,sd) / sqrt(n_subland) * 1.96),
#                                      dataset = rep(c("test", "train"), each=256),
#                                      added_o = factor(rep(apply(M_CoeffsIter_OrderAdded,1, median),2)))
# p3 <- ggplot(DS_RMSE_CoeffsIter_mean2[DS_RMSE_CoeffsIter_mean$coeffs<=250,], aes(x=coeffs, y=rmse)) + theme_classic() +
#   geom_errorbar(aes(ymin = rmse -rmse_ci, ymax = rmse + rmse_ci), width = 0.5) +
#   geom_point(size=2, aes(fill=dataset), shape=21) +
#   geom_point(data=DS_RMSE_CoeffsIter_mean[DS_RMSE_CoeffsIter_mean$rmse == min(DS_RMSE_CoeffsIter_mean$rmse[DS_RMSE_CoeffsIter_mean$dataset == "test"]),], color="#cc0000", shape=21, fill= "#7F3F98", size=3) +
#   theme(legend.position = c(0.8,0.8)) +
#   scale_fill_manual("", values = c("#7F3F98", "#FFDE17")) +
#   xlab("Number of coefficients added") +
#   ylab("RMSE")
# p3
# 
# p3b <- ggplot(DS_RMSE_CoeffsIter_mean2[DS_RMSE_CoeffsIter_mean$dataset =="test" & DS_RMSE_CoeffsIter_mean$coeffs<=20 ,], aes(x=coeffs, y=rmse)) + theme_classic() +
#   #geom_line(data=cbind(melt(RMSE_test_CoeffsIter[1:25,]), data.frame(land = rep(1:76, each=25))), aes(x=X1, y=value,group=land), alpha=0.2) +
#   geom_errorbar(aes(ymin = rmse -rmse_ci, ymax = rmse + rmse_ci), width = 0.2) +
#   geom_point(size=4, shape=21, aes(fill=added_o)) +
#   scale_fill_manual("Median order of\n coefficient added", values = colorRampPalette(colors = c("#f2ebf4", "#7F3F98"))(4)) +
#   xlab("Number of coefficients added") + scale_x_continuous(limits = c(0.9,20.1)) +
#   ylab("RMSE") + theme(legend.position = c(0.7,0.7) )
# p3b

# save.image("20180505_WeighedVsPartialWeighted.pdf")


# Accuracy by %var exp
M_CoeffsIter_VarExp_comb1 = do.call("cbind",lapply(1:n_subland, function(x){
  predobs = do.call("rbind", lapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    k = L_CoeffsIter_Folds[[x]][[y]]
    cbind(L_CoeffsIter_Pred[[x]][[y]][k,],obs = M8_Fit[k,x])
  }))
  test = 1-(colSums((predobs[,257]-predobs[,1:256])^2))/sum((predobs[,257]-mean(predobs[,257]))^2)
  predobs = do.call("rbind", lapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    k = L_CoeffsIter_Folds[[x]][[y]]
    cbind(L_CoeffsIter_Pred[[x]][[y]][-k,],obs = M8_Fit[-k,x])
  }))
  train = 1-(colSums((predobs[,257]-predobs[,1:256])^2))/sum((predobs[,257]-mean(predobs[,8]))^2)
  cbind(test,train)
}))

VarExp_test_CoeffsIter = M_CoeffsIter_VarExp_comb1[,rep(c(T,F),n_subland)]
VarExp_train_CoeffsIter = M_CoeffsIter_VarExp_comb1[,rep(c(F,T),n_subland)]
DS_VarExp_CoeffsIter_mean = data.frame(coeffs = c(1:256,1:256), rmse = c(rowMeans(VarExp_test_CoeffsIter), rowMeans(VarExp_train_CoeffsIter)),
                                     rmse_ci=c(apply(VarExp_test_CoeffsIter,1,sd) / sqrt(n_subland) * 1.96, apply(VarExp_train_CoeffsIter,1,sd) / sqrt(n_subland) * 1.96),
                                     dataset = rep(c("test", "train"), each=256),
                                     added_o = factor(rep(apply(M_CoeffsIter_OrderAdded,1, median),2)))
p3c <- ggplot(DS_VarExp_CoeffsIter_mean[DS_VarExp_CoeffsIter_mean$dataset =="test" & DS_VarExp_CoeffsIter_mean$coeffs<=20 ,], aes(x=coeffs, y=rmse)) + theme_classic() +
  #geom_line(data=cbind(melt(VarExp_test_CoeffsIter[1:25,]), data.frame(land = rep(1:76, each=25))), aes(x=X1, y=value,group=land), alpha=0.2) +
  geom_errorbar(aes(ymin = rmse -rmse_ci, ymax = rmse + rmse_ci)) +
  geom_point(size=3, shape=21, aes(fill=added_o)) +
  scale_fill_manual("Median order of\ncoefficient added", values = colorRampPalette(colors = c("#f2ebf4", "#7F3F98"))(4)) +
  xlab("Number of coefficients added") +
  ylab("RMSE") + theme(legend.position = c(0.7,0.4) )
p3c
ggsave("Figures-VarianceExplained/CV/WWalsh/20180320_CVWhalshWeighted_IterCoeffsTestMedianOrderCoeffAdded_VarExp.pdf",p3c, width = 6, height = 4)




#Add each order at a time but only using the significant coefficients
L_OrderIter_Pred_minRMSECoeffsIter = lapply(1:n_subland, function(x){
  cat("\r",x,"\t\t\t\t")
  mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    coeffs_M = do.call("cbind",lapply(0:6, function(oo){
      co = L_OrderIter_Coeffs[[x]][[y]]
      co[co[,2] > M_CoeffsIter_minRMSE[y,x] | coeff_orders_loci8 > oo,1] = 0
      co[,1]
    }))
    Zp = partialAv_matrix(genotypeToRemove = L_CoeffsIter_Folds[[x]][[y]], order = o)
    Wp = weighted_partial(order = o, partial_M = Zp)
    diag(Wp) = 1/diag(Wp)
    Wp[!is.finite(Wp)] = 0
    solve(M8_Walsh) %*% (Wp %*% coeffs_M)
  }, mc.cores = detectCores(all.tests = FALSE, logical = TRUE))
})

# Measure accuracy with RMSE
L_OrderIter_RMSE_minRMSECoeffsIter = lapply(1:n_subland, function(x){
  cat("\r",x,"\t\t\t\t")
  do.call("rbind",mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    k = L_CoeffsIter_Folds[[x]][[y]]
    pred = L_OrderIter_Pred_minRMSECoeffsIter[[x]][[y]][k,]
    obs = M8_Fit[k,x]
    test_err = sqrt(colMeans((obs - pred)^2))
    pred = L_OrderIter_Pred_minRMSECoeffsIter[[x]][[y]][-k,]
    obs = M8_Fit[-k,x]
    train_err = sqrt(colMeans((obs - pred)^2))
    cbind(test_err,train_err)
  },mc.cores=detectCores(all.tests = FALSE, logical = TRUE)))
})

L_OrderIter_RMSE_comb1_minRMSECoeffsIter =  do.call("cbind",lapply(1:n_subland, function(x){
  xx = L_OrderIter_RMSE_minRMSECoeffsIter[[x]]
  train_err = matrix(xx[,2],nrow=7)
  test_err = matrix(xx[,1],nrow=7)
  test = sqrt(rowSums(test_err^2) / length(L_CoeffsIter_Folds[[x]]))
  train = sqrt(rowSums(train_err^2) / length(L_CoeffsIter_Folds[[x]]))
  cbind(test,train)
}))

RMSE_minRMSECoeffsIter = L_OrderIter_RMSE_comb1_minRMSECoeffsIter[,rep(c(T,F),n_subland)]
DS_RMSE_minRMSECoeffsIter_mean = data.frame(orders = 0:6, RMSE_mean = rowMeans(RMSE_minRMSECoeffsIter), RMSE_ci = apply(RMSE_minRMSECoeffsIter,1,sd) / sqrt(n_subland) * 1.96)
DS_RMSE_minRMSECoeffsIter_all = melt(RMSE_minRMSECoeffsIter)
DS_RMSE_minRMSECoeffsIter_all$land = rep(1:76, each=7)
p5 <- ggplot(DS_RMSE_minRMSECoeffsIter_mean[DS_RMSE_minRMSECoeffsIter_mean$orders>=0,], aes(x=orders,y=RMSE_mean)) + theme_classic() +
  geom_line(linetype=2, alpha=0.5)  +
  #geom_line(data=DS_RMSE_minRMSECoeffsIter_all, aes(x=X1-1, y=value, group=land), alpha=0.1) +
  geom_errorbar(aes(ymin = RMSE_mean - RMSE_ci, ymax = RMSE_mean + RMSE_ci), width = 0.2) +
  geom_point(size=4, shape=21, aes(fill=factor(orders)), show.legend = F) +
  scale_fill_manual("", values = cols_orders) +
  ylab("RMSE") + xlab("Orders of coefficients used") + 
  scale_y_continuous(breaks = c(seq(0.08, 0.13, 0.01))) +
  scale_x_continuous(breaks = 0:6)
p5
ggsave("Figures-VarianceExplained/CV/WWalsh/20180320_CVWhalshWeighted_OrderIterSignifCoeffs_RMSEvsOrderAdded.pdf", p5, width = 4, height = 3.7)

grid.arrange(p3b, p5 + theme(axis.title.y = element_blank()) , layout_matrix = rbind(c(1,1,1,2,2), c(1,1,1,2,2))  )

# Measure accuracy with %VE
L_OrderIter_VarExp_minRMSECoeffsIter = lapply(1:n_subland, function(x){
  cat("\r",x,"\t\t\t\t")
  do.call("rbind",mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    k = L_CoeffsIter_Folds[[x]][[y]]
    pred = L_OrderIter_Pred_minRMSECoeffsIter[[x]][[y]][k,]
    obs = M8_Fit[k,x]
    test_ve = 1-(colSums((obs-pred)^2))/sum((obs-mean(obs))^2)
    pred = L_OrderIter_Pred_minRMSECoeffsIter[[x]][[y]][-k,]
    obs = M8_Fit[-k,x]
    train_ve = 1-(colSums((obs-pred)^2))/sum((obs-mean(obs))^2)
    cbind(test_ve,train_ve)
  },mc.cores=detectCores(all.tests = FALSE, logical = TRUE)))
})

L_OrderIter_VarExp_comb1_minRMSECoeffsIter =  do.call("cbind",lapply(1:n_subland, function(x){
  # xx = L_OrderIter_VarExp_minRMSECoeffsIter[[x]]
  # train_ev = matrix(xx[,2],nrow=7)
  # test_ev = matrix(xx[,1],nrow=7)
  # test = apply(test_ev, 1, mean)
  # train = apply(train_ev, 1, mean)
  # cbind(test,train)
  predobs = do.call("rbind", lapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    k = L_CoeffsIter_Folds[[x]][[y]]
    cbind(L_OrderIter_Pred_minRMSECoeffsIter[[x]][[y]][k,],obs = M8_Fit[k,x])
  }))
  test = 1-(colSums((predobs[,8]-predobs[,1:7])^2))/sum((predobs[,8]-mean(predobs[,8]))^2)
  predobs = do.call("rbind", lapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    k = L_CoeffsIter_Folds[[x]][[y]]
    cbind(L_OrderIter_Pred_minRMSECoeffsIter[[x]][[y]][-k,],obs = M8_Fit[-k,x])
  }))
  train = 1-(colSums((predobs[,8]-predobs[,1:7])^2))/sum((predobs[,8]-mean(predobs[,8]))^2)
  cbind(test,train)
}))

VarExp_minRMSECoeffsIter = L_OrderIter_VarExp_comb1_minRMSECoeffsIter[,rep(c(T,F),n_subland)]
DS_VarExp_minRMSECoeffsIter_mean = data.frame(orders = 0:6, VarExp_mean = rowMeans(VarExp_minRMSECoeffsIter), VarExp_ci = apply(VarExp_minRMSECoeffsIter,1,sd) / sqrt(n_subland) * 1.96)
DS_VarExp_minRMSECoeffsIter_all = melt(VarExp_minRMSECoeffsIter)
DS_VarExp_minRMSECoeffsIter_all$land = rep(1:76, each=7)
DS_VarExp_minRMSECoeffsIter_all$X1 = rep(1:7, 76)
p5b <- ggplot(DS_VarExp_minRMSECoeffsIter_mean[DS_VarExp_minRMSECoeffsIter_mean$orders >=0,], aes(x=orders,y=VarExp_mean)) + theme_classic() +
  geom_line(linetype=2)  +
  #geom_line(data=DS_VarExp_minRMSECoeffsIter_all, aes(x=X1-1, y=value, group=land), alpha=0.1) +
  geom_errorbar(aes(ymin = VarExp_mean - VarExp_ci, ymax = VarExp_mean + VarExp_ci), width = 0.1) +
  geom_point(size=3, shape=21, aes(fill=factor(orders)), show.legend = F) +
  scale_fill_manual("", values = cols_orders) +
  ylab("% variance explained") + xlab("Orders of coefficients used") +
  scale_x_continuous(breaks = 0:6)
p5b
ggsave("Figures-VarianceExplained/CV/WWalsh/20180320_CVWhalshWeighted_OrderIterSignifCoeffs_VarExpvsOrderAdded.pdf", p5b, width = 6, height = 5)

temp = VarExp_minRMSECoeffsIter
temp[1,] = 0

GainVarExp_minRMSECoeffsIter = data.frame(gain = rowMeans(apply(temp, 2, diff)), order = factor(1:6))
ggplot(GainVarExp_minRMSECoeffsIter, aes(x=order, y=gain)) + theme_classic() +
  geom_col(aes(fill=order)) + 
  scale_fill_manual("Order", values = cols_orders) +
  ylab("% of variance explained") + xlab("") + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
ggplot(GainVarExp_minRMSECoeffsIter, aes(x=order, y=gain)) + theme_classic() +
  geom_col(aes(fill=order)) + 
  scale_fill_manual("Order", values = cols_orders) +
  ylab("% of variance explained") + xlab("") + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())


DS_CoeffsIter_BestModel =  setNames(data.frame(do.call("rbind", lapply(1:n_subland, function(x){
  cat("\r",x,"\t\t\t\t")
  do.call("rbind",mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    k = L_CoeffsIter_Folds[[x]][[y]]
    pred = L_CoeffsIter_Pred[[x]][[y]][k,M_CoeffsIter_minRMSE[y,x]]
    obs = M8_Fit[k,x]
    cbind(pred,obs)
  }, mc.cores = 4))
}))), c("Pred", "Obs"))


min_best = floor(min(as.vector(as.matrix(DS_CoeffsIter_BestModel)))*100)/100
max_best = ceiling(max(as.vector(as.matrix(DS_CoeffsIter_BestModel)))*100)/100
p4 <- ggplot(DS_CoeffsIter_BestModel, aes(x=Pred,y=Obs)) + theme_classic() +
  geom_abline(linetype=2, alpha=0.5) +
  stat_bin_hex(bins = 50) + scale_fill_gradient("Genotypes counts", low = "grey85", high = "grey15") +
  xlab("Predicted ln(fitness)") + ylab("Observed ln(fitness)") +
  theme(legend.position = "") +
  annotate("text", x=-0.65, y=-0.05, label=paste0("RMSE = ", round(RMSE(pred = DS_CoeffsIter_BestModel$Pred, obs = DS_CoeffsIter_BestModel$Obs, n=length(DS_CoeffsIter_BestModel$Obs)),2))) +
  annotate("text", x=-0.65, y=-0.15, label=paste0("%VE = ", round(VarExp(pred = DS_CoeffsIter_BestModel$Pred, obs = DS_CoeffsIter_BestModel$Obs),2))) +
  scale_x_continuous(limits = c(min_best, 0.1)) +
  scale_y_continuous(limits = c(min_best, 0.1))
p4
ggsave("Figures/20180320_CVWhalshWeighted_BestModel_PredvsObs.pdf", p4, width = 6, height = 6)


cor(DS_CoeffsIter_BestModel$Pred, DS_CoeffsIter_BestModel$Obs,method = "spearman")


### Use the same median number of coefficients used in each sublandscape and predict fitness using epistatic coefficients relative to a single background
MedianNumCoeffs_sublands = floor(apply(M_CoeffsIter_minRMSE,2, median))

# Coefficients and rank for each 10fold CV - Usable whenever iterating by order or by coefficients
L_Relative_Coeffs = lapply(1:n_subland, function(x){
  W_Coeffs = M8_Relative %*% M8_Fit[,x]
  W_SE = sqrt(abs(M8_Relative) %*% M8_SE[,x]^2)   
  W_tstat = (W_Coeffs/W_SE)/2^(coeff_orders_loci8) 
  cbind(W_Coeffs,rank(-abs(W_tstat)))
})

M_Relative_OrderAdded = do.call("cbind", lapply(1:n_subland, function(x){
  coeffs_ranks = L_Relative_Coeffs[[x]][,2]
  coeff_orders_loci8[order(coeffs_ranks)]
}))

L_Relative_Pred = lapply(1:n_subland, function(x){
  cat("\r",x,"\t\t\t\t")
  coeffs_M = do.call("cbind", mclapply(1:256, function(z){
    coeffs_ranks = L_Relative_Coeffs[[x]]
    coeffs_ranks[coeffs_ranks[,2] > z,1] <- 0
    coeffs_ranks[,1]
  }, mc.cores = 4))
  solve(M8_Relative) %*%  coeffs_M
})

M_Relative_Pred_Order1 = do.call("cbind", lapply(1:n_subland, function(x){
  coeffs_ranks = L_Relative_Coeffs[[x]]
  coeffs_ranks[coeff_orders_loci8 > 1,1] <- 0
  coeffs_M = coeffs_ranks[,1]
  solve(M8_Relative) %*%  coeffs_M
}))

M_Relative_RMSE = do.call("cbind", lapply(1:n_subland, function(x){
    pred = L_Relative_Pred[[x]]
    obs = M8_Fit[,x]
    test_err = sqrt(colMeans((obs - pred)^2))
}))

DS_RMSE_Relative_mean = data.frame(coeffs = c(1:256,1:256), rmse = rowMeans(M_Relative_RMSE),
                                     rmse_ci= apply(M_Relative_RMSE,1,sd) / sqrt(n_subland) * 1.96,
                                     dataset = rep(c("test", "train"), each=256),
                                     added_o = factor(rep(apply(M_Relative_OrderAdded,1, median),2)))

p3b <- ggplot(DS_RMSE_Relative_mean, aes(x=coeffs, y=rmse)) + theme_classic() +
  #geom_line(data=cbind(melt(M_Relative_RMSE), data.frame(land = rep(1:76, each=256))), aes(x=X1, y=value,group=land), alpha=0.2) +
  geom_errorbar(aes(ymin = rmse -rmse_ci, ymax = rmse + rmse_ci)) +
  geom_point(size=2, shape=21, aes(fill=added_o)) +
  scale_fill_manual("Median order of\n coefficient added", values = colorRampPalette(colors = c("#f2ebf4", "#7F3F98"))(8)) +
  xlab("Number of coefficients added") +
  ylab("RMSE") + theme(legend.position = c(0.9,0.7) )
p3b

DS_Relative_BestModel = setNames(data.frame(do.call("rbind", lapply(1:n_subland, function(x){
  cbind(L_Relative_Pred[[x]][,MedianNumCoeffs_sublands[x]], M8_Fit[,x])
}))), c("Pred", "Obs"))

min_best_rel = floor(min(as.vector(as.matrix(DS_Relative_BestModel)))*100)/100
max_best_rel = ceiling(max(as.vector(as.matrix(DS_Relative_BestModel)))*100)/100
p4b <- ggplot(DS_Relative_BestModel, aes(x=Pred,y=Obs)) + theme_classic() +
  geom_abline(linetype=2, alpha=0.5) +
  stat_bin_hex(bins = 40) + scale_fill_gradient("Genotypes counts", low = "grey85", high = "grey15") +
  xlab("Predicted ln(fitness)") + ylab("Observed ln(fitness)") +
  theme(legend.position = "") +
  annotate("text", x=1.5, y=-0.75, label=paste0("RMSE = ", round(RMSE(pred = DS_Relative_BestModel$Pred, obs = DS_Relative_BestModel$Obs, n=length(DS_Relative_BestModel$Obs)),2))) +
  annotate("text", x=1.5, y=-0.9, label=paste0("%VE = ", round(VarExp(pred = DS_Relative_BestModel$Pred, obs = DS_Relative_BestModel$Obs),2))) +
  scale_x_continuous(limits = c(min_best_rel, max_best_rel)) +
  scale_y_continuous(limits = c(-1, 0.2))
p4b


DS_Relative_ModelOrder1 = data.frame(Pred = as.vector(M_Relative_Pred_Order1), Obs=as.vector(M8_Fit))
min_best_rel_o1 = floor(min(as.vector(as.matrix(DS_Relative_ModelOrder1)))*100)/100
max_best_rel_o1 = ceiling(max(as.vector(as.matrix(DS_Relative_ModelOrder1)))*100)/100
p4d <- ggplot(DS_Relative_ModelOrder1, aes(x=Pred,y=Obs)) + theme_classic() +
  geom_abline(linetype=2, alpha=0.5) +
  stat_bin_hex(bins = 50) + scale_fill_gradient("Genotypes counts", low = "grey85", high = "grey15") +
  xlab("Predicted ln(fitness)") + ylab("Observed ln(fitness)") +
  theme(legend.position = "") +
  annotate("text", x=-0.2, y=-0.85, label=paste0("RMSE = ", round(RMSE(pred = DS_Relative_ModelOrder1$Pred, obs = DS_Relative_ModelOrder1$Obs, n=length(DS_Relative_ModelOrder1$Obs)),2))) +
  annotate("text", x=-0.2, y=-0.95, label=paste0("%VE = ", round(VarExp(pred = DS_Relative_ModelOrder1$Pred, obs = DS_Relative_ModelOrder1$Obs),2))) +
  scale_x_continuous(limits = c(min_best_rel_o1, 0.1)) +
  scale_y_continuous(limits = c(min_best_rel_o1, 0.1))
p4d

# VarExp(pred = DS_Relative_BestModel$Pred, obs = DS_Relative_BestModel$Obs)


### Get the Relative using first only up to 1st order coefficient with and without the relative matrix

L_Order1_Pred = lapply(1:n_subland, function(x){
  cat("\r",x,"\t\t\t\t")
  mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    co = L_OrderIter_Coeffs[[x]][[y]]
    co[coeff_orders_loci8 > 1,1] = 0
    coeffs_M = co[,1]
    Zp = partialAv_matrix(genotypeToRemove = L_CoeffsIter_Folds[[x]][[y]], order = o)
    Wp = weighted_partial(order = o, partial_M = Zp)
    diag(Wp) = 1/diag(Wp)
    Wp[!is.finite(Wp)] = 0
    solve(M8_Walsh) %*% (Wp %*% coeffs_M)
  }, mc.cores = detectCores(all.tests = FALSE, logical = TRUE))
})



L_Order1_RMSE = lapply(1:n_subland, function(x){
  cat("\r",x,"\t\t\t\t")
  do.call("rbind",mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    k = L_CoeffsIter_Folds[[x]][[y]]
    pred = L_Order1_Pred[[x]][[y]][k,]
    obs = M8_Fit[k,x]
    test_err = sqrt(mean((obs - pred)^2))
    pred = L_Order1_Pred[[x]][[y]][-k,]
    obs = M8_Fit[-k,x]
    train_err = sqrt(mean((obs - pred)^2))
    cbind(test_err,train_err)
  }, mc.cores = 4))
})

M_Order1_RMSE_comb1 = do.call("rbind",lapply(1:n_subland, function(x){
  xx = L_Order1_RMSE[[x]]
  train_err = xx[,2]
  test_err = xx[,1]
  test = sqrt(sum(test_err^2) / length(L_CoeffsIter_Folds[[x]]))
  train = sqrt(sum(train_err^2) / length(L_CoeffsIter_Folds[[x]]))
  cbind(test,train)
}))

DS_Order1_AvModel =  setNames(data.frame(do.call("rbind", lapply(1:n_subland, function(x){
  do.call("rbind",mclapply(1:length(L_CoeffsIter_Folds[[x]]), function(y){
    k = L_CoeffsIter_Folds[[x]][[y]]
    pred = L_Order1_Pred[[x]][[y]][k,]
    obs = M8_Fit[k,x]
    cbind(pred,obs)
  }, mc.cores = 4))
}))), c("Pred", "Obs"))
cor(DS_Order1_AvModel$Pred, DS_Order1_AvModel$Obs, method = "spearman")



p4c <- ggplot(DS_Order1_AvModel, aes(x=Pred,y=Obs)) + theme_classic() +
  geom_abline(linetype=2, alpha=0.5) +
  stat_bin_hex(bins = 50) + scale_fill_gradient("Genotypes counts", low = "grey85", high = "grey15") +
  xlab("Predicted ln(fitness)") + ylab("Observed ln(fitness)") +
  theme(legend.position = "") +
  annotate("text", x=-0.65, y=-0.05, label=paste0("RMSE = ", round(RMSE(pred = DS_Order1_AvModel$Pred, obs = DS_Order1_AvModel$Obs, n=length(DS_Order1_AvModel$Obs)),2))) +
  annotate("text", x=-0.65, y=-0.15, label=paste0("%VE = ", round(VarExp(pred = DS_Order1_AvModel$Pred, obs = DS_Order1_AvModel$Obs),2))) +
  scale_x_continuous(limits = c(min_best, 0.1)) +
  scale_y_continuous(limits = c(min_best, 0.1))
p4c

p4_PredVsObs = grid.arrange(p4d, p4c, p4, nrow = 1)
ggsave("Figures-VarianceExplained/CV/WWalsh/20180321_AllPredVsObs_3plots.pdf", p4_PredVsObs, width = 10, height = 3.1)
