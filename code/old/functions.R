######################################################################################
# Functions from https://github.com/GRousselet/onsetsim/blob/main/code/functions.R
# Generate permutation distributions of t values for independent groups
# cond1 and cond2 are Nt trials x Nf time frames matrices
######################################################################################

# Cluster correction
# Form clusters using binary vector: pvals < alpha
# clusters must be at least 2 time frames
cluster.make <- function (x) {
    
    y <- rle(x)
    cmap <- vector(mode = "numeric", length = 0)
    nC <- length(y$values) # number of clusters
    indx <- 0 # cluster counter
    
    for (CL in 1:nC) {
        
    if(y$values[CL] == 0 || y$lengths[CL] == 1){
      val <- 0
    } else {
      indx <- indx + 1
      val <- indx
    }
    cmap <- c(cmap, rep(val, y$lengths[CL]))
    }
    
    cmap
    
}

# ## Save sum for each cluster
# # values = statistics (t values)
# # cmap = output from cluster.make
# cluster.sum <- function(values, cmap){
#   csum <- vector(mode = "numeric", length = max(cmap))
#   if(max(cmap)>0){
#     for(CL in 1:max(cmap)){
#       csum[CL] <- sum(values[cmap==CL])
#     }
#   } else {
#     csum <- 0
#   }
#   csum
# }
# 
# # Cluster test
# # values = statistics (t values)
# # cmap = output from cluster.make
# # boot.th = bootstrap quantile of distribution of max cluster sums
# cluster.test <- function(values, cmap, boot.th){
#   csig <- vector(mode = "logical", length = length(cmap))
#   if(max(cmap)>0){
#     for(CL in 1:max(cmap)){
#       csig[cmap==CL] <- sum(values[cmap==CL]) > boot.th
#     }
#   }
#   csig
# }

sim.counter <- function(S, nsim, inc){
  if(S == 1){
    # print(paste(nsim,"iterations:",S))
    cat(nsim,"iterations:",S)
    #beep(2)
  }
  if(S %% inc == 0){
    # print(paste("iteration",S,"/",nsim))
    cat(" /",S)
    #beep(2)
  }
}

# INPUTS:
#   sigmask = a vector of 0s and 1s indicating points below and above threshold.
#   Xf      = a vector of time points, matching sigmask.
#   rmzero  = option to discard an onset if it belongs to a cluster that starts in the baseline, including zero; default = TRUE.
#
# OUTPUT:
#   onset   = latency of the first cluster in the units of Xf.
# find_onset <- function (sigmask, Xf, rmzero = TRUE) {
#     
#     onset <- NA
#     onset <- try(Xf[which(sigmask)[1]], silent = TRUE)
#     
#     # removing un-realistic clusters that includes zero
#     if (rmzero) {
#         
#         if (is.finite(onset) ) {
#             
#             if (onset == 0) {
#                 
#                 cmap <- cluster.make(sigmask)
#                 cmap[cmap==1] <- 0
#                 onset <- Xf[which(sigmask * (cmap>0)>0)[1]]
#                 
#             }
#             
#         }
#         
#     }
#     
#     # returning the onset
#     return (onset)
#     
# }

# https://www.statology.org/mode-in-r/
# find_mode <- function(x, na.rm = TRUE) {
#   if(na.rm){
#     x <- na.omit(x)
#   }
#   u <- unique(x)
#   tab <- tabulate(match(x, u))
#   u[tab == max(tab)]
# }
# 
# keeporder <- function(x){
#   x <- as.character(x)
#   x <- factor(x, levels=unique(x))
#   x
# }

# code from Rand Wilcox https://osf.io/xhe8u/
# hd<-function(x,q=.5,na.rm=TRUE){
#   #
#   #  Compute the Harrell-Davis estimate of the qth quantile
#   #
#   #  The vector x contains the data,
#   #  and the desired quantile is q
#   #  The default value for q is .5.
#   #
#   if(na.rm)x=elimna(x)
#   n<-length(x)
#   m1<-(n+1)*q
#   m2<-(n+1)*(1-q)
#   vec<-seq(along=x)
#   w<-pbeta(vec/n,m1,m2)-pbeta((vec-1)/n,m1,m2)  # W sub i values
#   y<-sort(x)
#   hd<-sum(w*y)
#   hd
# }

# elimna<-function(m){
#   #
#   # remove any rows of data having missing values
#   #
#   DONE=FALSE
#   if(is.list(m) && is.matrix(m)){
#     z=pool.a.list(m)
#     m=matrix(z,ncol=ncol(m))
#     DONE=TRUE
#   }
#   if(!DONE){
#     if(is.list(m) && is.matrix(m[[1]])){
#       for(j in 1:length(m))m[[j]]=na.omit(m[[j]])
#       e=m
#       DONE=TRUE
#     }}
#   if(!DONE){
#     if(is.list(m) && is.null(dim(m))){ #!is.matrix(m))
#       for(j in 1:length(m))m[[j]]=as.vector(na.omit(m[[j]]))
#       e=m
#       DONE=TRUE
#     }}
#   if(!DONE){
#     m<-as.matrix(m)
#     ikeep<-c(1:nrow(m))
#     for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
#     e<-m[ikeep[ikeep>=1],]
#   }
#   e
# }
