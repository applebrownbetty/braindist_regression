.libPaths("/nas/longleaf/home/chalmer/bin/R3.6.0")

# install.packages('pracma',repo="http://cran.rstudio.com/")

library(tidyverse)
library(MASS)
library(pdSpecEst)
library(correlation)
library(pracma)

## Wrap with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#task_ID <- 0

################################################
##### ----- Step 1: simulate dataset ----- #####
################################################

UT <- function(x) {
  return(x[upper.tri(x)])
}

jaccard_wt <- function(x,y) {
  # calculates weighted Jaccard similarity between
  # x and y (positive vectors of equal length)
  if(length(x)!=length(y)) {
    warning("vectors must be of the same length")
  } else if(any(c(x,y) < 0)) {
    warning("all elements must be nonnegative")
  } else {
    return(sum(pmin(x,y)^2)/sum(pmax(x,y)^2))
  }
}

binary_coefs <- function(x,y) {
  # calculates weighted Jaccard similarity between
  # x and y (positive vectors of equal length)
  if(length(x)!=length(y)) {
    warning("vectors must be of the same length")
  } else if(any(!(c(x,y) %in%  c(0,1)))) {
    warning("all elements must be 0 or 1")
  } else {
    M_11 <- sum(x==1 & y==1)
    M_01 <- sum(x==0 & y==1)
    M_10 <- sum(x==1 & y==0)
    M_00 <- sum(x==0 & y==0)
    JI <- (M_11)/(M_01+M_10+M_11)
    JD <- 1-JI
    vec <- c(JD)
    names(vec) <- c("JD")    
    return(t(vec))
  }
}

cont_coefs <- function(x,y) {
  # calculates weighted Jaccard similarity between
  # x and y (positive vectors of equal length)
  if(length(x)!=length(y)) {
    warning("vectors must be of the same length")
  } else if(any((c(x,y) < 0))) {
    warning("all elements should be >= 0")
  } else {
    # jaccard_wt <- jaccard_wt(x,y)
    ks <- ks.test(x,y)$statistic
    euc <- dist(rbind(x,y),method = "euclidean")
    # maxim <- dist(rbind(x,y),method = "maximum")
    # manhat <- dist(rbind(x,y),method = "manhattan")
    # canber <- dist(rbind(x,y),method = "canberra")
    # minkow0.5 <- dist(rbind(x,y),method = "minkowski",p=0.5)
    # minkow1.4 <- dist(rbind(x,y),method = "minkowski",p=1.4)
    # minkow3 <- dist(rbind(x,y),method = "minkowski",p=3)
    # vec <- c(jaccard_wt,ks,euc,maxim,manhat,canber,minkow0.5,minkow1.4,minkow3)
    # names(vec) <- c("Jaccard_wt","KS","Euclidean","Maximum","Manhattan","Canberra","Minkowski_0.5", "Minkowski_1.4", "Minkowski_3")
    vec <- c(ks,euc)
    names(vec) <- c("KS","Euclidean")
    return(t(vec))
  }
}

band1 <- function(n,val) {
  m <- diag(0, n)
  m[abs(row(m) - col(m)) <= 1] <- 1
  m <- m + (1-m)*(1-val)
  return(m)
}

randcorr <- function(n) {
  Q <- pracma::randortho(n)
  D <- diag(50*rbeta(n,.075,2))
  A <- t(Q) %*% D %*% Q
  invsqA <- diag(diag(A)^(-1/2))
  A <- cor_smooth(invsqA %*% A %*% invsqA,verbose=F,tol=10^-6)
  return(A)
}

aim2simdat <- function(seed) {
 
  set.seed(seed)
  
  #########################################
  ##### ------ DEMOGRAPHIC DATA ----- #####
  #########################################
  
  np <- 100 # minimum number of participants
  nn <- 268 # number of nodes in each subject's network
  
  FamID <- rep(0,np+6)
  TwinFamID <- rep(0,np+6)
  counterFamID = 1
  counterID = 1
  while(counterID <= np) {
    famsize <- sample(2:6,1,prob=c(.3,.5,.15,.035,.015))
    FamID[counterID:(counterID + famsize - 1)] <- counterFamID
    TwinFamID[counterID:(counterID + famsize - 1)][sample(1:famsize,2)] <- counterFamID
    counterFamID = counterFamID + 1
    counterID = counterID + famsize
  }
  
  np <- counterID - 1
  nf <- counterFamID - 1
  FamID <- FamID[1:np]
  TwinFamID <- TwinFamID[1:np]
  
  dat <- tibble(
    ID = 1:np,
    FamID = FamID,
    TwinFamID = TwinFamID,
    Age = NA,
    IQ = NA,
    Sex = NA,
    Trt = NA
  )
  
  dat <- dat[order(dat$FamID,dat$TwinFamID),]
  dat$ID <- paste0("ID_",dat$ID)
  dat$FamID <- paste0("FamID_",dat$FamID)
  dat$TwinFamID <- ifelse(dat$TwinFamID == 0, dat$FamID, paste0("TwinFamID_",dat$TwinFamID))
  

  wts <- seq(100,0,-20)
    
  df <-as.data.frame(matrix(nrow=choose(np,2)*length(wts)*6,ncol=1))
  colnames(df) <- "Sim"
  counter = 1
  
  for(sm in 1:6) {
    if(sm==1) {
        cov_corr_family <- 0
        cov_corr_twin <- 0
        
        corr_family <- 0
        corr_twin <- 0
    }
    if(sm==2) {
        cov_corr_family <- 0
        cov_corr_twin <- 0
        
        corr_family <- 0.5
        corr_twin <- 0.8
    }
    if(sm==3) {
        cov_corr_family <- 0.5
        cov_corr_twin <- 0.8
        
        corr_family <- 0
        corr_twin <- 0
    }
    if(sm==4) {
        cov_corr_family <- 0.5
        cov_corr_twin <- 0.8
        
        corr_family <- 0.5
        corr_twin <- 0.8
    }
    if(sm==5) {
        cov_corr_family <- 1
        cov_corr_twin <- 1
        
        corr_family <- 0.5
        corr_twin <- 0.8
    }
    if(sm==6) {
        cov_corr_family <- 0.5
        cov_corr_twin <- 1
        
        corr_family <- 0.5
        corr_twin <- 0.8
    }   
    set.seed(seed)
    for(fid in unique(dat$FamID)) {
    loopids <- dat$ID[dat$FamID == fid]
    Sigma <- diag(length(loopids)) + (1-diag(length(loopids)))*cov_corr_family
    twintmp <- which(dat$TwinFamID[dat$FamID == fid]==paste0("Twin",fid))
    Sigma[twintmp[1], twintmp[2]] <- cov_corr_twin
    Sigma[twintmp[2], twintmp[1]] <- cov_corr_twin
    dimnames(Sigma) <- rep(list(loopids),2)
    quants <- pnorm(MASS::mvrnorm(n = 4, mu = rep(0,length(loopids)), Sigma = Sigma, tol = 1e-06, empirical = FALSE))
    for(id in loopids) {
      dat$Sex[dat$ID == id] = ifelse(qbinom(quants[1,id],1,.5)==1,"Male","Female")
      dat$Trt[dat$ID == id] = ifelse(qbinom(quants[2,id],1,.5)==1,"Treatment","Placebo")
      dat$IQ[dat$ID == id] = round(qnorm(quants[3,id],mean=100,sd=15))
      dat$Age[dat$ID == id] = round(qnorm(quants[4,id],mean=100,sd=15))
    }
  }
  
  # pick out binary and continuous variables - 
  binary_vars <- cbind(select_if(dat, is.character),select_if(dat, is.factor))
  binary_vars <- binary_vars[,!(names(binary_vars) %in% c("ID","FamID","TwinFamID"))]
  contin_vars <- select_if(dat,is.numeric)
  contin_vars <- contin_vars[,!(names(contin_vars) %in% c("ID","FamID","TwinFamID"))]
  
  # 
  amod <- pmin(5.95,pmax(-5.95,(dat$IQ-100)*.15 + ifelse(dat$Trt=="Treatment",2,-2)))
  names(amod) <- dat$ID
  rand_noise <- randcorr(nn)
  rand_noise <- rand_noise[lower.tri(rand_noise)]
  signal_noise <- function(q) return((1-wt/100)*abs(sample(rand_noise,length(q))) + wt/100*qbeta(q, 7+amod[id], 7-amod[id]))
  
  for(wt in wts) {
  
    ##############################################
    ##### ------ CONNECTIVITY MATRICES ----- #####
    ##############################################
  
    CM_weight <- CM_top005 <- matrix(nrow = np,ncol = choose(nn,2),dimnames = list(dat$ID,1:choose(nn,2)))
  
    for(fid in unique(dat$FamID)) {
      famdat <- dat %>% filter(FamID == fid)
      
      Sigma <- diag(dim(famdat)[1]) + (1-diag(dim(famdat)[1]))*corr_family
      twintmp <- which(famdat$TwinFamID==paste0("Twin",fid))
      Sigma[twintmp[1], twintmp[2]] <- corr_twin
      Sigma[twintmp[2], twintmp[1]] <- corr_twin
      dimnames(Sigma) <- rep(list(famdat$ID),2)
      quants <- pnorm(MASS::mvrnorm(n = 3, mu = rep(0,dim(famdat)[1]), Sigma = Sigma, tol = 1e-06, empirical = FALSE))
      for(id in famdat$ID) {      
        tmpmat <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
        tmpmat[,21:35] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[1,id]))
        tmpmat[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[2,id],5,5))
        tmpmat[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[3,id],5,5))
        tmpmat <- cor_smooth(cor(tmpmat),verbose = F,tol=10^-6)
        assign(paste0("CM_",id),tmpmat)
        CM_weight[id,] <- UT(tmpmat)
        CM_top005[id,] <- CM_weight[id,] > quantile(CM_weight[id,],1-.005)
      }
    }
    
    ##################################
    ##### ------ DISTANCES ----- #####
    ##################################
      
    for(i in 1:(dim(dat)[1]-1)) {
      for(j in (i+1):dim(dat)[1]) {
        df[counter,"Sim"] <- sm
        df[counter,"Weight"] <- wt
        df[counter,c("IDA","FamIDA","TwinFamIDA")] <- dat[i,c("ID","FamID","TwinFamID")]
        df[counter,c("IDB","FamIDB","TwinFamIDB")] <- dat[j,c("ID","FamID","TwinFamID")]
        df[counter,names(binary_vars)] <-  1*(binary_vars[which(dat$ID==df$IDA[counter]),] != binary_vars[which(dat$ID==df$IDB[counter]),]) 
        df[counter,names(contin_vars)] <- abs(contin_vars[which(dat$ID==df$IDA[counter]),] - contin_vars[which(dat$ID==df$IDB[counter]),])
#        df[counter,"KS"] <- ks.test(
#          CM_weight[df$IDA[counter],],
#          CM_weight[df$IDB[counter],]
#        )$statistic
        df[counter,"Jaccard"] <- binary_coefs(
          CM_top005[df$IDA[counter],],
          CM_top005[df$IDB[counter],])
        df[counter,"Euclidean"] <- dist(rbind(
          CM_weight[df$IDA[counter],],
          CM_weight[df$IDB[counter],]),
          method = "euclidean")
        df[counter,"logEuclidean"] <- pdDist(
          get(paste0("CM_",df$IDA[counter])), 
          get(paste0("CM_",df$IDB[counter])), 
          metric = "logEuclidean")
        counter = counter +1
      }
    }
  }
}
output_file_name <- paste0("dyn_fam_",sprintf("%06d",seed),".Rds")
output_file_path <- file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets/",output_file_name)
saveRDS(df,output_file_path)
}

aim2simdat(task_ID)

#if(file.exists(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets/",output_file_name))) {
#  print(simid)
#} else {
#  aim2simdat(simid)
#}