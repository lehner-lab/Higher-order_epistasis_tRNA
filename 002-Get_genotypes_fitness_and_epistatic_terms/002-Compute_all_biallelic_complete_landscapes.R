### Computing all biallelic complete sublandscapes from 1 to 8 mutations ###
# (A) Single function to obtain all complete biallelic landscapes of n number of mutations
# (B) Summary statistics of the complete landscapes (Extended Figure 1d)


# JDE and GD, September 2017

# Load data
setwd("002-Get_genotypes_fitness_and_epistatic_terms//")
load(file = "../001-Data/001-GenotypesFitness.RData")

library(parallel)
library(plyr)


## (A) Function to get the IDs of all variants with n mutations and all their intermediate in all backgrounds
getIDs = function(n, var, authorized_pos=1:10){
  dummy = rep(0,10)
  pos_list = t(combn(1:10,n))
  # For each variant in the list (actually a matrix)
  ids = do.call("rbind",mclapply(as.list(1:nrow(pos_list)),function(co){
    co = pos_list[co,]
    tri.allelic.mat = matrix(FALSE, nrow=16 , ncol = 10)
    tri.allelic.mat[,c(2,3,7,9)] = as.matrix(expand.grid(c(T,F),c(T,F),c(T,F),c(T,F)))
    do.call("rbind", lapply(1:nrow(tri.allelic.mat),function(lv){            
      logical.vector=tri.allelic.mat[lv,]
      # all combinations of these n mutants
      if(n == 1){
        all.comb = matrix(dummy,ncol=10)
        all.comb[,co] = 1
        all.comb[,logical.vector & all.comb == 1] = 2
      }else{
        all.comb = do.call("rbind",sapply(1:n,function(level){
          t(apply(combn(co,level),2,function(intermediate){
            dummy2 = dummy
            dummy2[intermediate] = 1
            dummy2[logical.vector & dummy2 == 1] = 2
            dummy2
          }))
        }))
      } 
      back.pos = authorized_pos[!(authorized_pos %in% co)]
      
      # identify all the possible backgrounds in which they may occur
      back = do.call("rbind",sapply(0:(length(authorized_pos)-n),function(pos){
        # backbones with 'pos' mutations
        b = t(combn(back.pos,pos))
        t(apply(b,1,function(back){
          dummy3 = dummy
          dummy3[back] = 1
          dummy3[logical.vector & dummy3 == 1] = 2
          dummy3
        }))
      }))
      
      # keep only backgrounds that have been observed
      back2 = do.call("paste0",as.data.frame(back))
      back = back[back2 %in% var,]
      back2 = back2[back2 %in% var]
      
      # check if all intermediates are present between this triple mutant and a given background
      # First add the id vector of the background to all intermediates
      present1 = t(do.call("cbind",lapply(as.list(1:nrow(back)),function(b){
        b = back[b,]
        b + t(all.comb)
      })))
      # Then paste the 21 positions together to get a vector of IDs
      present = do.call("paste0",as.data.frame(present1))
      # check if these are present and make a logical matrix where each row represents one landscape and addup the TRUEs to get how many variants from this ladnscape are present
      present2 = colSums(matrix(present %in% var, nrow=nrow(all.comb)))
      k = length(present2[present2 == nrow(all.comb)])
      if(k>1){
        # make a matrix with the IDs and keep only the ones where all the variants are present
        present = cbind(back2,t(matrix(present, nrow=nrow(all.comb))))[present2 == nrow(all.comb),]
        # add as a first column the number of mutations and then the positions of these mutations
        present = cbind(rep(n,times=k),rep(paste(co,collapse="-"),times=k), rep(paste(all.comb[nrow(all.comb),co], collapse = "-"),times=k),present)
      }else if(k == 1){
        present = cbind(back2,t(matrix(present, nrow=nrow(all.comb))))[present2 == nrow(all.comb),]
        present = matrix(c(n,paste(co,collapse="-"),paste(all.comb[nrow(all.comb),co], collapse = "-"),present),nrow=1)
      }else{
        NULL
      }
      
    }))
  },mc.cores=4))
}

# Iterate over 1:8 to get all possible complete landscapes from any background
L_DS = lapply(1:8, function(x){
  M = getIDs(n=x, var = P$var)
  M = unique(M)
  comb = rep(0:x, choose(x, 0:x))
  var_names = c(c("num_mut", "pos_mut", "vars_mut"), paste("id_", comb, "_", 1:length(comb), sep=""))
  DS = setNames(data.frame(M), var_names)
})


## (B) Get a summary of the combinatiorial of complete landscapes and the number of bckgrounds where each n combination of mutations are found
DS_summary = do.call("rbind", lapply(1:8, function(x){
  n_bckg = count(L_DS[[x]], c(1:3))$freq
  data.frame(n = x, n_comb=length(n_bckg), n_landscapes = dim(L_DS[[x]])[1], min_bckg = min(n_bckg), median_bckg = median(n_bckg), max_bckg = max(n_bckg))
}))
DS_summary
# write.table(DS_summary, sep = "\t", col.names = T, row.names = F, file = "../001-Data/Summary_complete_landscapes.txt", quote = F)


# Save each dataset 
for (I in 1:8) {
  file_name = paste0("../001-Data/CompleteLands_n", I, ".RData")
  CL = L_DS[[I]]
  save(CL, file = file_name)
}
