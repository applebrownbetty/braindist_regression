.libPaths("/nas/longleaf/home/chalmer/bin/R3.6.0")
library(permuco)
library(tidyverse)
library(MDMR)

task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#task_ID <- 1

df <- readRDS(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets",paste0("sim_",sprintf("%06d",task_ID),".Rds")))
IDmat <- model.matrix(~ID1-1,data=df)+model.matrix(~ID2-1,data=df)
colnames(IDmat) <- paste0("ID_",1:100)
set.seed(task_ID)


##### -- SUBJECT INFO -- #####
np <- 100 # 20 total participants (make this an even number)
nn <- 268 # number of nodes in each subject's network

dat <- tibble(
  ID = as.factor(1:np),
  Sex = ifelse(rbinom(np,1,.5)==1,"Male","Female"),
  Trt = sample(rep(c("Treatment","Placebo"),np/2)),
  IQ = round(rnorm(np,mean=100,sd=15)),
  Age = round(rnorm(np,mean=100,sd=15))
)

output_file_name <- paste0("perm_JAC_",sprintf("%06d",task_ID),".Rds")
output_file_path <- file.path("/pine/scr/c/h/chalmer/aim2/simulations/results/",output_file_name)

GLiSh <- function(form,DM,ID1,ID2,df) {
  df <- as.data.frame(df)
  
  y <- as.vector(df[,DM])
  X <- model.matrix(form,data = df)
  
  lmfit <- lm(y~X-1)
  beta <- coef(lmfit)
  
  lvls <- unique(c(df[,ID1],df[,ID2]))
  IDmat <- as(
    model.matrix(~factor(df[,ID1],levels=lvls)-1) + 
      model.matrix(~factor(df[,ID2],levels=lvls)-1),
    "dgCMatrix")
  IDvarmat <- IDmat %*% Matrix::t(IDmat) - 2*diag(nrow(IDmat))
  r <- y - X %*% beta
  s2 <- 1/dim(df)[1]*sum(r^2)
  l2 <- 1/sum(IDvarmat)*sum(t(r) %*% IDvarmat %*% r)
  V <- s2*diag(dim(df)[1]) + l2*IDvarmat
  cholV <- Matrix::chol(V)
  cholVinv <- Matrix::chol2inv(cholV) 
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

quad <- function(a, b, c)
{
  a <- as.complex(a)
  answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
              (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  if(all(Im(answer) == 0)) answer <- Re(answer)
  if(answer[1] == answer[2]) return(answer[1])
  answer
}

symmat_diag0 <- function(X) {
  X[lower.tri(X)]  <- t(X)[lower.tri(X)]
  diag(X) <- 0
  return(X)
}

VDM <- function(a) {
  n <- max(abs(quad(1,1,-2*length(a))))
  b <- matrix (0,n,n)
  b[lower.tri(b, diag=FALSE)] <- a
  
  return(symmat_diag0(t(b)))
}

sims <- 1:4
signals <- seq(0,100,10)
distances <- c("Jaccard","KS","Minkowski_0.5","Canberra","Minkowski_1.4","Euclidean","Minkowski_3")
form=~Sex+Trt+IQ+Age
nps <- 10000

dat_results <- data.frame(ID = rep(task_ID,length(sims)*length(signals)*length(distances)*6))

pmat <- Pmat(np = nps,n=4950)
adjpmat <- genpmat(nps,df)

counter = 1
for(sim in sims) {
  for(signal in signals) {
  #print(counter)
    for(distance in distances) {
      DM <- paste0("Sim",sim,"_wt",sprintf("%03d",signal),ifelse(distance == "Jaccard","_top20_","_corr_"),distance)
      y <- as.vector(as.data.frame(df)[,DM])
      if(distance=="KS") y <- log(y)
      
      # Permutation
      dftmp <- cbind(y,df[,c("Age","IQ","Sex","Trt")])
      xxx <- lmperm(y ~ Age+IQ+Sex+Trt, data = dftmp,P=adjpmat,method="freedman_lane")
      dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"Permutation"))
      dat_results[counter,paste0(rownames(xxx$table),"_p")] <- xxx$table$`permutation Pr(>|t|)`
      counter <- counter + 1

      # F-Test
      dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"F-Test"))
      dat_results[counter,paste0(rownames(xxx$table),"_p")] <- xxx$table$`parametric Pr(>|t|)`
      dat_results[counter,paste0(rownames(xxx$table),"_est")] <- xxx$table$Estimate
      dat_results[counter,paste0(rownames(xxx$table),"_se")] <- xxx$table$`Std. Error`
      dat_results[counter,paste0(rownames(xxx$table),"_t")] <- xxx$table$`t value`
      counter <- counter + 1

      # FGLS
      dftmp <- GLiSh(~Age+IQ+Sex+Trt,"y","ID1","ID2",cbind(dftmp,df[,c("ID1","ID2")]))
      dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"FGLS"))
      dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- t(summary(lm(y ~ .-1,data=dftmp))$coefficients[c("Age","IQ","Sex","Trt"),"Pr(>|t|)"])
      counter <- counter + 1
      
      # F-Test with ILE
      dftmp <- cbind(y,df[,c("Age","IQ","Sex","Trt")],IDmat)
      dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"F-Test with I.L.E."))
      dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- t(summary(lm(y ~ .-1,data=dftmp))$coefficients[c("Age","IQ","Sex","Trt"),"Pr(>|t|)"])
      counter <- counter + 1
      
      # MDMR Permutation
      D <- VDM(y)
      invisible(capture.output(suppressWarnings(mdmr.res <- mdmr(X = dat[,2:5], D = D,perm.p = TRUE,nperm=10000))))
      x <- mdmr.res$pv
      dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"MDMR Permutation"))
      dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- x[c("Age","IQ","SexMale","TrtTreatment"),"perm.pvals"]
      counter = counter+1
      
      # MDMR Analytic
      invisible(capture.output(suppressWarnings(mdmr.res <- mdmr(X = dat[,2:5], D = D,perm.p = FALSE))))
      x <- mdmr.res$pv
      dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"MDMR Analytic"))
      dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- x[c("Age","IQ","SexMale","TrtTreatment"),"analytic.pvals"]
      counter = counter+1
    }
  }
  saveRDS(dat_results,output_file_path)
}
