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
  
  wts <- seq(100,0,-10)
  
  corr_time_band1 <- 0
  corr_withintask <- 0.7/(1-corr_time_band1)
  corr_withingroup <- 0.3/(1-corr_time_band1)
  tasks <- paste0("task",rep(1:3,each=4),"_",rep(1:4,times=3))
  
  A <- diag(4) + corr_withintask*(1-diag(4))
  B <- corr_withingroup*diag(4)
  errormat <- kronecker(diag(1, 3), A) + kronecker(1-diag(1, 3), B)
  dimnames(errormat) <- rep(list(tasks),2)
  
  
  ##### -- SUBJECT INFO -- #####
  np <- 30 # 40 total participants (make this an even number)
  nn <- 268 # number of nodes in each subject's network
  
  dat <- tibble(
    ID = as.factor(1:np),
    Sex = ifelse(rbinom(np,1,.5)==1,"Male","Female"),
    Trt = ifelse(rbinom(np,1,.5)==1,"Treatment","Placebo"),
    IQ = round(rnorm(np,mean=100,sd=15)),
    Age = round(rnorm(np,mean=100,sd=15))
  )
  
  # pick out binary and continuous variables - 
  binary_vars <- cbind(select_if(dat, is.character),select_if(dat, is.factor))
  binary_vars <- binary_vars[,!(names(binary_vars) %in% "ID")]
  contin_vars <- select_if(dat,is.numeric)
  contin_vars <- contin_vars[,!(names(contin_vars) %in% "ID")]
  
  amod <- pmin(5.95,pmax(-5.95,(dat$IQ-100)*.15 + ifelse(dat$Trt=="Treatment",2,-2)))
  task_time <- matrix(nrow = np,ncol = length(tasks),dimnames = list(1:np,tasks))
  invisible(sapply(paste0(tasks,"_list"), assign, list(), envir = parent.frame()))
  invisible(sapply(paste0(tasks,"_deg"), assign, matrix(nrow = np,ncol = choose(nn,2),dimnames = list(1:np,1:choose(nn,2))) , envir = parent.frame()))
  invisible(sapply(paste0(tasks,"_top005"), assign, matrix(nrow = np,ncol = choose(nn,2),dimnames = list(1:np,1:choose(nn,2))) , envir = parent.frame()))
  invisible(sapply(paste0(tasks,"_top20"), assign, matrix(nrow = np,ncol = choose(nn,2),dimnames = list(1:np,1:choose(nn,2))) , envir = parent.frame()))
  
  loopdat <- tibble(
    ID = rep(1:np,each=length(tasks)),
    Task = rep(rep(1:3,each=4),times=np),
    Group = rep(1:4,times=np*3)
  )
  
  df <-as.data.frame(matrix(nrow=choose(np*4,2)*length(wts)*3,ncol=1))
  colnames(df) <- "IDA"
  counter = 1
  
  for(wt in wts) {
    for(id in 1:np) {
      # fill matrix with noise
      task1_1 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      task1_2 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      task1_3 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      task1_4 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      task2_1 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      task2_2 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      task2_3 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      task2_4 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      task3_1 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      task3_2 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      task3_3 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      task3_4 <- mvrnorm(n=2500,mu=rep(0,nn),Sigma = randcorr(nn))
      
      rand_task <- sample(1:length(tasks))
      errormat_i <- errormat[rand_task,rand_task]
      errormat_i <- errormat_i*band1(length(tasks),corr_time_band1)
      
      quants <- pnorm(MASS::mvrnorm(n = 3, mu = rep(0,length(tasks)), Sigma = errormat, tol = 1e-06, empirical = FALSE))
      square1 <- 1
      square2 <- 2
      square3 <- 3
           
      rand_noise <- randcorr(nn)
      rand_noise <- rand_noise[lower.tri(rand_noise)]
      signal_noise <- function(q) return((1-wt/100)*abs(sample(rand_noise,length(q))) + wt/100*qbeta(q, 7+amod[id], 7-amod[id]))
      
      task1_1[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square1,"task1_1"],5,5))
      task1_2[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square1,"task1_2"],5,5))
      task1_3[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square1,"task1_3"],5,5))
      task1_4[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square1,"task1_4"],5,5))
      task2_1[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square1,"task2_1"],5,5))
      task2_2[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square1,"task2_2"],5,5))
      task2_3[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square1,"task2_3"],5,5))
      task2_4[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square1,"task2_4"],5,5))
      task3_1[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square1,"task3_1"]))
      task3_2[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square1,"task3_2"]))
      task3_3[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square1,"task3_3"]))
      task3_4[,21:35] <-   mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square1,"task3_4"]))
      
      task1_1[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square2,"task1_1"],5,5))
      task1_2[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square2,"task1_2"],5,5))
      task1_3[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square2,"task1_3"],5,5))
      task1_4[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*qbeta(quants[square2,"task1_4"],5,5))
      task2_1[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square2,"task2_1"]))
      task2_2[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square2,"task2_2"]))
      task2_3[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square2,"task2_3"]))
      task2_4[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square2,"task2_4"]))
      task3_1[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square2,"task3_1"]))
      task3_2[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square2,"task3_2"]))
      task3_3[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square2,"task3_3"]))
      task3_4[,111:125] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square2,"task3_4"]))
      
      task1_1[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task1_1"]))
      task1_2[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task1_2"]))
      task1_3[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task1_3"]))
      task1_4[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task1_4"]))
      task2_1[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task2_1"]))
      task2_2[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task2_2"]))
      task2_3[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task2_3"]))
      task2_4[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task2_4"]))
      task3_1[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task3_1"]))
      task3_2[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task3_2"]))
      task3_3[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task3_3"]))
      task3_4[,201:215] <- mvrnorm(n=2500,mu=rep(0,15),Sigma = diag(15) + (1-diag(15))*signal_noise(quants[square3,"task3_4"]))
      
      # turn into correlation matrices
      task1_1 <- cor_smooth(cor(task1_1),verbose = F,tol=10^-6)
      task1_2 <- cor_smooth(cor(task1_2),verbose = F,tol=10^-6)
      task1_3 <- cor_smooth(cor(task1_3),verbose = F,tol=10^-6)
      task1_4 <- cor_smooth(cor(task1_4),verbose = F,tol=10^-6)
      task2_1 <- cor_smooth(cor(task2_1),verbose = F,tol=10^-6)
      task2_2 <- cor_smooth(cor(task2_2),verbose = F,tol=10^-6)
      task2_3 <- cor_smooth(cor(task2_3),verbose = F,tol=10^-6)
      task2_4 <- cor_smooth(cor(task2_4),verbose = F,tol=10^-6)
      task3_1 <- cor_smooth(cor(task3_1),verbose = F,tol=10^-6)
      task3_2 <- cor_smooth(cor(task3_2),verbose = F,tol=10^-6)
      task3_3 <- cor_smooth(cor(task3_3),verbose = F,tol=10^-6)
      task3_4 <- cor_smooth(cor(task3_4),verbose = F,tol=10^-6)
      
      task_time[id,] <- rand_task
      
      task1_1_list[[id]] <- task1_1
      task1_2_list[[id]] <- task1_2
      task1_3_list[[id]] <- task1_3
      task1_4_list[[id]] <- task1_4
      task2_1_list[[id]] <- task2_1
      task2_2_list[[id]] <- task2_2
      task2_3_list[[id]] <- task2_3
      task2_4_list[[id]] <- task2_4
      task3_1_list[[id]] <- task3_1
      task3_2_list[[id]] <- task3_2
      task3_3_list[[id]] <- task3_3
      task3_4_list[[id]] <- task3_4
      
      task1_1_deg[id,] <- UT(task1_1)
      task1_2_deg[id,] <- UT(task1_2)
      task1_3_deg[id,] <- UT(task1_3)
      task1_4_deg[id,] <- UT(task1_4)
      task2_1_deg[id,] <- UT(task2_1)
      task2_2_deg[id,] <- UT(task2_2)
      task2_3_deg[id,] <- UT(task2_3)
      task2_4_deg[id,] <- UT(task2_4)
      task3_1_deg[id,] <- UT(task3_1)
      task3_2_deg[id,] <- UT(task3_2)
      task3_3_deg[id,] <- UT(task3_3)
      task3_4_deg[id,] <- UT(task3_4)
      
      task1_1_top005[id,] <- task1_1_deg[id,] > quantile(task1_1_deg[id,],1-.005)
      task1_2_top005[id,] <- task1_2_deg[id,] > quantile(task1_2_deg[id,],1-.005)
      task1_3_top005[id,] <- task1_3_deg[id,] > quantile(task1_3_deg[id,],1-.005)
      task1_4_top005[id,] <- task1_4_deg[id,] > quantile(task1_4_deg[id,],1-.005)
      task2_1_top005[id,] <- task2_1_deg[id,] > quantile(task2_1_deg[id,],1-.005)
      task2_2_top005[id,] <- task2_2_deg[id,] > quantile(task2_2_deg[id,],1-.005)
      task2_3_top005[id,] <- task2_3_deg[id,] > quantile(task2_3_deg[id,],1-.005)
      task2_4_top005[id,] <- task2_4_deg[id,] > quantile(task2_4_deg[id,],1-.005)
      task3_1_top005[id,] <- task3_1_deg[id,] > quantile(task3_1_deg[id,],1-.005)
      task3_2_top005[id,] <- task3_2_deg[id,] > quantile(task3_2_deg[id,],1-.005)
      task3_3_top005[id,] <- task3_3_deg[id,] > quantile(task3_3_deg[id,],1-.005)
      task3_4_top005[id,] <- task3_4_deg[id,] > quantile(task3_4_deg[id,],1-.005)
            
      task1_1_top20[id,] <- task1_1_deg[id,] > quantile(task1_1_deg[id,],1-.2)
      task1_2_top20[id,] <- task1_2_deg[id,] > quantile(task1_2_deg[id,],1-.2)
      task1_3_top20[id,] <- task1_3_deg[id,] > quantile(task1_3_deg[id,],1-.2)
      task1_4_top20[id,] <- task1_4_deg[id,] > quantile(task1_4_deg[id,],1-.2)
      task2_1_top20[id,] <- task2_1_deg[id,] > quantile(task2_1_deg[id,],1-.2)
      task2_2_top20[id,] <- task2_2_deg[id,] > quantile(task2_2_deg[id,],1-.2)
      task2_3_top20[id,] <- task2_3_deg[id,] > quantile(task2_3_deg[id,],1-.2)
      task2_4_top20[id,] <- task2_4_deg[id,] > quantile(task2_4_deg[id,],1-.2)
      task3_1_top20[id,] <- task3_1_deg[id,] > quantile(task3_1_deg[id,],1-.2)
      task3_2_top20[id,] <- task3_2_deg[id,] > quantile(task3_2_deg[id,],1-.2)
      task3_3_top20[id,] <- task3_3_deg[id,] > quantile(task3_3_deg[id,],1-.2)
      task3_4_top20[id,] <- task3_4_deg[id,] > quantile(task3_4_deg[id,],1-.2)
    }
    for(i in 1:(dim(loopdat)[1]-1)) {
      for(j in (i+1):dim(loopdat)[1]) {
        if(loopdat[i,"Task"] == loopdat[j,"Task"]) {
          df[counter,c("IDA","Task_IDA","Group_IDA")] <- loopdat[i,c("ID","Task","Group")]
          df[counter,c("IDB","Task_IDB","Group_IDB")] <- loopdat[j,c("ID","Task","Group")]
          df[counter,"Time_IDA"] <- task_time[df$IDA[counter],df$Task_IDA[counter]]
          df[counter,"Time_IDB"] <- task_time[df$IDB[counter],df$Task_IDB[counter]]
          df[counter,names(binary_vars)] <-  binary_vars[df$IDA[counter],] != binary_vars[df$IDB[counter],]
          #df[counter,"Sex"] <- 1*(binary_vars[df$IDA[counter],"Sex"]=="Male") - 1*(binary_vars[df$IDB[counter],"Sex"]=="Male")
          #df[counter,"Trt"] <- 1*(binary_vars[df$IDA[counter],"Trt"]=="Treatment") - 1*(binary_vars[df$IDB[counter],"Trt"]=="Treatment")
          df[counter,names(contin_vars)] <- abs(contin_vars[df$IDA[counter],] -  contin_vars[df$IDB[counter],])
          df[counter,"Weight"] <- wt
          df[counter,"KS"] <- ks.test(
            get(paste0("task",df$Task_IDA[counter],"_",df$Group_IDA[counter],"_deg"))[df$IDA[counter],],
            get(paste0("task",df$Task_IDB[counter],"_",df$Group_IDB[counter],"_deg"))[df$IDB[counter],]
          )$statistic
          df[counter,"JD_top005"] <- binary_coefs(
            get(paste0("task",df$Task_IDA[counter],"_",df$Group_IDA[counter],"_top005"))[df$IDA[counter],],
            get(paste0("task",df$Task_IDB[counter],"_",df$Group_IDB[counter],"_top005"))[df$IDB[counter],])
          df[counter,"JD_top20"] <- binary_coefs(
            get(paste0("task",df$Task_IDA[counter],"_",df$Group_IDA[counter],"_top20"))[df$IDA[counter],],
            get(paste0("task",df$Task_IDB[counter],"_",df$Group_IDB[counter],"_top20"))[df$IDB[counter],])
          df[counter,"EUC"] <- dist(rbind(
            get(paste0("task",df$Task_IDA[counter],"_",df$Group_IDA[counter],"_deg"))[df$IDA[counter],],
            get(paste0("task",df$Task_IDB[counter],"_",df$Group_IDB[counter],"_deg"))[df$IDB[counter],]),
            method = "euclidean")
          # df[counter,"AIRM"] <- pdDist(
          #    get(paste0("task",df$Task_IDA[counter],"_",df$Group_IDA[counter],"_list"))[[df$IDA[counter]]], 
          #    get(paste0("task",df$Task_IDB[counter],"_",df$Group_IDB[counter],"_list"))[[df$IDB[counter]]], 
          #    metric = "Riemannian")
          df[counter,"LERM"] <- pdDist(
            get(paste0("task",df$Task_IDA[counter],"_",df$Group_IDA[counter],"_list"))[[df$IDA[counter]]], 
            get(paste0("task",df$Task_IDB[counter],"_",df$Group_IDB[counter],"_list"))[[df$IDB[counter]]], 
            metric = "logEuclidean")
          counter = counter +1
          }
      }
    }
  }
  output_file_name <- paste0("aim2mixed_sim_",sprintf("%06d",seed),".Rds")
  output_file_path <- file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets/",output_file_name)
  saveRDS(df,output_file_path)
}
aim2simdat(task_ID)
