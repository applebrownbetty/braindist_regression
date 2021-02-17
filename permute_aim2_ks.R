.libPaths("/nas/longleaf/home/chalmer/bin/R3.6.0")
library(permuco)

task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#task_ID <- 1

set.seed(task_ID)

nps <- 10000



form=~Sex+Trt+IQ+Age

GLiSh <- function(form,DM,ID1,ID2,df) {
  df <- as.data.frame(df)
  
  y <- as.vector(df[,DM])
  X <- model.matrix(form,data = df)
  
  lmfit <- lm(y~X-1)
  beta <- coef(lmfit)
  
  lvls <- unique(c(df[,ID1],df[,ID2]))
  
#  IDmat <- model.matrix(~factor(df[,ID1],levels=lvls)-1) + 
#    model.matrix(~factor(df[,ID2],levels=lvls)-1)
#  IDvarmat <- IDmat %*% t(IDmat) - 2*diag(nrow(IDmat))
  
#  IDmat <- as(
#            model.matrix(~factor(df[,ID1],levels=lvls)-1) + 
#            model.matrix(~factor(df[,ID2],levels=lvls)-1),
#            "dgCMatrix")
#  IDvarmat <- IDmat %*% Matrix::t(IDmat) - 2*diag(nrow(IDmat))
#  
#  r <- y - X %*% beta
#  s2 <- 1/dim(df)[1]*sum(r^2)
#  l2 <- 1/sum(IDvarmat)*sum(t(r) %*% IDvarmat %*% r)
#  V <- s2*diag(dim(df)[1]) + l2*IDvarmat
##  cholV <- chol(V)
##  cholVinv <- chol2inv(cholV) 
##  return(as.data.frame(as.matrix(cholVinv %*% cbind(y,X))))
#  cholV <- Matrix::chol(V)
#  cholVinv <- Matrix::chol2inv(cholV) 
#  return(as.data.frame(as.matrix(cholVinv %*% cbind(y,X))))

  IDmat <- model.matrix(~factor(df[,ID1],levels=lvls)-1) + 
    model.matrix(~factor(df[,ID2],levels=lvls)-1)
  IDvarmat <- IDmat %*% t(IDmat) - 2*diag(nrow(IDmat))
  
  r <- y - X %*% beta
  s2 <- 1/dim(df)[1]*sum(r^2)
  l2 <- 1/sum(IDvarmat)*sum(t(r) %*% IDvarmat %*% r)
  V <- s2*diag(dim(df)[1]) + l2*IDvarmat
#  cholV <- chol(V)
#  cholVinv <- chol2inv(cholV) 
#  return(as.data.frame(as.matrix(cholVinv %*% cbind(y,X))))
  cholV <- chol(V)
  cholVinv <- chol2inv(cholV) 
  return(as.data.frame(as.matrix(cholVinv %*% cbind(y,X))))
  
}

genpmat <- function (nperms,df) {
  mat <- matrix(nrow=100,ncol=100)
  for(i in 1:4950) {
    mat[df$ID1[i],df$ID2[i]] <- i
  }
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  
  pmat <- matrix(nrow=4950,ncol=nperms)
  
  matpmat <- Pmat(nperms,100)
  
  for(i in 1:nperms) {
    X <- matpmat[,i]
    mat_p <- mat[X,X]
    pmat[,i] <- t(mat_p)[lower.tri(mat_p)]
  }
  return(as.Pmat(pmat))
}

sims <- 4
signals <- seq(0,100,10)
distances <- c("Jaccard","KS","Minkowski_0.5","Canberra","Minkowski_1.4","Euclidean","Minkowski_3")
#meths <- c("draper_stoneman","freedman_lane","manly","terBraak","kennedy","dekker")
meths <- c("freedman_lane")


#simIDs <- (task_ID*10+1):(task_ID*10 + 10)
#for(simID in simIDs) {
simID <- task_ID

  df <- readRDS(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets",paste0("sim4_",sprintf("%06d",simID),".Rds")))
  output_file_name <- paste0("perm4_",sprintf("%06d",simID),".Rds")
  output_file_path <- file.path("/pine/scr/c/h/chalmer/aim2/simulations/permutations/",output_file_name)
  
  dat_results <- data.frame(ID = rep(simID,length(sims)*length(signals)*length(distances)*2*(length(meths)+1)))
  
  pmat <- Pmat(np = nps,n=4950)
  adjpmat <- genpmat(nps,df)
  
  counter = 1
  for(sim in sims) {
    for(signal in signals) {
    #print(counter)
      for(distance in distances) {
        DM <- paste0("Sim",sim,"_wt",sprintf("%03d",signal),ifelse(distance == "Jaccard","_top20_","_corr_"),distance)
        print(DM)
        dftmp <- df[,c(DM,"Age","IQ","Sex","Trt")]
        for(meth in meths) {
          xxx <- lmperm(get(DM) ~ Age+IQ+Sex+Trt, data = dftmp,P=adjpmat,method=meth)
          dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",simID),sim,signal,distance,meth))
          dat_results[counter,paste0(rownames(xxx$table),"_p")] <- xxx$table$`permutation Pr(>|t|)`
          counter <- counter + 1
        }
  
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",simID),sim,signal,distance,"parametric"))
        dat_results[counter,paste0(rownames(xxx$table),"_p")] <- xxx$table$`parametric Pr(>|t|)`
        dat_results[counter,paste0(rownames(xxx$table),"_est")] <- xxx$table$Estimate
        dat_results[counter,paste0(rownames(xxx$table),"_se")] <- xxx$table$`Std. Error`
        dat_results[counter,paste0(rownames(xxx$table),"_t")] <- xxx$table$`t value`
        counter <- counter + 1
  
        dftmp <- GLiSh(~Age+IQ+Sex+Trt,DM,"ID1","ID2",df)
        
        for(meth in meths) {
          xxx <- lmperm(y ~ ., data = dftmp,P=pmat,method=meth)
          dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",simID),sim,signal,distance,paste0(meth,"_adj")))
          dat_results[counter,paste0(rownames(xxx$table),"_p")] <- xxx$table$`permutation Pr(>|t|)`
          counter <- counter + 1
        }
        
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",simID),sim,signal,distance,"parametric_adj"))
        dat_results[counter,paste0(rownames(xxx$table),"_p")] <- xxx$table$`parametric Pr(>|t|)`
        dat_results[counter,paste0(rownames(xxx$table),"_est")] <- xxx$table$Estimate
        dat_results[counter,paste0(rownames(xxx$table),"_se")] <- xxx$table$`Std. Error`
        dat_results[counter,paste0(rownames(xxx$table),"_t")] <- xxx$table$`t value`
        counter <- counter + 1
      }
    }
    saveRDS(dat_results,output_file_path)
}
#}
