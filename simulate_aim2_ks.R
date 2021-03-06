library(tidyverse)
## Wrap with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#task_ID <- 0

################################################
##### ----- Step 1: simulate dataset ----- #####
################################################

jaccard_index <- function(x,y) {
  # calculates Jaccard Index between
  # x and y (binary vectors of equal length)
  J <- sum(x==y)/length(x)
  return(J)
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

# For these coefficients we follow the notation of
#  Sneath and Sokal (1973): let "U" equal the number of
#  disagreements (either "s 1 99-" 0 " or "0"- 1 "), "a" equal
#  the number of shared presences, "d" equal the number
#  of shared absences. "N" is the number of characters. For
#  a particular pair of individuals a character can then be
#  referred to as an "a-type" character, for example, if both
#  individuals have the I state.

binary_coefs <- function(x,y) {
  # calculates weighted Jaccard similarity between
  # x and y (positive vectors of equal length)
  if(length(x)!=length(y)) {
    warning("vectors must be of the same length")
  } else if(any(!(c(x,y) %in%  c(0,1)))) {
    warning("all elements must be 0 or 1")
  } else {
    U <- sum(x!=y)
    a <- sum(x==1 & y==1)
    d <- sum(x==0 & y==0)
    N <- length(x)
    Jaccard <- (a+d)/N
    vec <- c(Jaccard)
    names(vec) <- c("Jaccard")
#    Russel_Rao <- a/N
#    idkJaccard <- a/(a+U)
#    Kulzinski <- 2*a/(2*a+U+2*d)
#    BaroniUrbani_Bauer <- (a + sqrt(a*d))/(a + U + sqrt(a*d))
#    vec <- c(Jaccard,Russel_Rao,idkJaccard,Kulzinski,BaroniUrbani_Bauer)
#    names(vec) <- c("Jaccard","Russel_Rao","idkJaccard","Kulzinski","BaroniUrbani_Bauer")
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
#    jaccard_wt <- jaccard_wt(x,y)
    ks <- ks.test(x,y)$statistic
#    maxim <- dist(rbind(x,y),method = "maximum")
#    manhat <- dist(rbind(x,y),method = "manhattan")
    minkow0.5 <- dist(rbind(x,y),method = "minkowski",p=0.5)
    canber <- dist(rbind(x,y),method = "canberra")
    minkow1.4 <- dist(rbind(x,y),method = "minkowski",p=1.4)
    euc <- dist(rbind(x,y),method = "euclidean")
    minkow3 <- dist(rbind(x,y),method = "minkowski",p=3)
#    vec <- c(jaccard_wt,ks,euc,maxim,manhat,canber,minkow0.5,minkow1.4,minkow3)
#    names(vec) <- c("Jaccard_wt","KS","Euclidean","Maximum","Manhattan","Canberra","Minkowski_0.5", "Minkowski_1.4", "Minkowski_3")
    vec <- c(ks,minkow0.5,canber,minkow1.4,euc,minkow3)
    names(vec) <- c("KS","Minkowski_0.5","Canberra","Minkowski_1.4","Euclidean","Minkowski_3")
    return(t(vec))
  }
}

aim2simdat <- function(seed) {
  
  set.seed(seed)
  
  ##### -- SUBJECT INFO -- #####
  np <- 100 # 20 total participants (make this an even number)
  nn <- 268 # number of nodes in each subject's network
  
  nn_r1 <- 15 # number of nodes in region 1 (signal nodes regardless of covariates)
  nn_r2 <- 15 # number of nodes in region 2 (signal nodes based on covariates - see amod)
  # region 0 (noise) is everything not in region 1 or 2
  
  dat <- tibble(
    ID = as.factor(1:np),
    Sex = ifelse(rbinom(np,1,.5)==1,"Male","Female"),
    Trt = sample(rep(c("Treatment","Placebo"),np/2)),
    IQ = round(rnorm(np,mean=100,sd=15)),
    Age = round(rnorm(np,mean=100,sd=15))
  )
  
  # initialize a data frame to keep track of all pairwise values
  df <- as_tibble(matrix(dat$ID[1],nrow = choose(np,2), ncol = 2))
  colnames(df) <- c("ID1","ID2")
  df$ID1 <- factor(df$ID1,levels=unique(dat$ID))
  df$ID2 <- factor(df$ID2,levels=unique(dat$ID))
  
  # pick out binary and continuous variables - 
  binary_vars <- cbind(select_if(dat, is.character),select_if(dat, is.factor))
  binary_vars <- binary_vars[,!(names(binary_vars) %in% "ID")]
  contin_vars <- select_if(dat,is.numeric)
  contin_vars <- contin_vars[,!(names(contin_vars) %in% "ID")]
  
  #counter for current indice of vector JM
  counter <- 1
  
  for(i in 1:(np-1)) {
    for(j in (i+1):np) {
      df[counter,"ID1"] <- dat$ID[i]
      df[counter,"ID2"] <- dat$ID[j]
      df[counter,names(binary_vars)] <-  1*(binary_vars[i,] != binary_vars[j,]) 
      df[counter,names(contin_vars)] <- abs(contin_vars[i,] - contin_vars[j,])
      counter <- counter + 1
    }
  }
  
  amod <- pmax(4/3,4 + (dat$IQ-100)*.2 + (dat$Trt=="Treatment")*6)
  amod_ks <- .2*(amod - mean(amod))/sd(amod)
  
  #connmats <- list()
  sim4_deg <- matrix(nrow = np,ncol = nn,dimnames = list(1:np,1:nn))  
  sim4_top20 <- matrix(nrow = np,ncol = nn,dimnames = list(1:np,1:nn))
  
  for(wt in seq(100,0,-10)) {
#  for(wt in 100) {
    for(id in 1:np) {
      sim4_deg[id,] <- rnorm(nn,mean=wt/100*amod_ks[id]+100,sd=1)
      sim4_top20[id,] <- sim4_deg[id,] > quantile(sim4_deg[id,],1-.2)
    }
    
    counter <- 1
    for(i in 1:(np-1)) {
      for(j in (i+1):np) {
        binary_indices <- binary_coefs(sim4_top20[i,],sim4_top20[j,])
        cont_indices <- cont_coefs(sim4_deg[i,],sim4_deg[j,])
        df[counter,paste0("Sim4_wt",sprintf("%03d",wt),"_top20_",colnames(binary_indices))] <- binary_indices
        df[counter,paste0("Sim4_wt",sprintf("%03d",wt),"_corr_",colnames(cont_indices))] <- cont_indices
        counter <- counter + 1
      }
    }
  }
  
  output_file_name <- paste0("sim4_",sprintf("%06d",seed),".Rds")
  output_file_path <- file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets/",output_file_name)
  saveRDS(df,output_file_path)
}

#neededfiles <- paste0("sim4_",sprintf("%06d",1:10^5),".Rds")
#donefiles <- list.files("/pine/scr/c/h/chalmer/aim2/simulations/datasets/", all.files=FALSE, full.names=FALSE)
#xxx <- neededfiles[!(neededfiles %in% donefiles)]
#xxx <- as.numeric(gsub('^....|....$', '', xxx))
#print(xxx)
#
#for(task_ID in xxx) aim2simdat(task_ID)

#if(file.exists(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets/","simsleft.Rds"))) {
#  iddat <- readRDS(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets/","simsleft.Rds"))
#  simids <- iddat[iddat$task==task_ID,]$id
#  for(simid in simids) {
#    aim2simdat(simid)
#  }
#} else {
#  for(i in 1:100) {
#  simid <- task_ID*100+i
#  aim2simdat(simid)
#}
#}

for(i in 1:100) {
simid <- task_ID*100+i
aim2simdat(simid)
}