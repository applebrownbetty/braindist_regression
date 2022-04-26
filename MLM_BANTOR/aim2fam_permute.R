.libPaths("/nas/longleaf/home/chalmer/bin/R3.6.0")
library(tidyverse)
library(lme4)
library(lmerTest)
library(MDMR)

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


# https://stackoverflow.com/questions/13762224/how-to-sort-files-list-by-date
details = file.info(list.files(path="/pine/scr/c/h/chalmer/aim2/simulations/datasets/",pattern="*.Rds", full.names = TRUE))
details = details[with(details, order(as.POSIXct(mtime))), ]
files = rownames(details)
files = sub(".Rds","",sub("/pine/scr/c/h/chalmer/aim2/simulations/datasets//dyn_fam_", "", files))

task_ID <- as.numeric(files[as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))])
#task_ID <- 1
#task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(task_ID)

simIDs <- task_ID
signals <- seq(0,100,20)
distances <- c("Jaccard","Euclidean","logEuclidean")
sims <- 1:6


dat_results <- data.frame(ID = rep(simIDs,each=7*length(signals)*length(distances)*length(sims)))

output_file_name <- paste0("results_aim2fam_",sprintf("%06d",simIDs[1]),".Rds")
output_file_path <- file.path("/pine/scr/c/h/chalmer/aim2/simulations/results/",output_file_name)

  ##### -- SUBJECT INFO -- ##### load individual level covariates for MDMR method
  #########################################
  ##### ------ DEMOGRAPHIC DATA ----- #####
  #########################################
  set.seed(task_ID)
  np <- 100 # minimum number of participants
  
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
  

counter = 1
for(simID in simIDs) {
  df <- readRDS(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets/",paste0("dyn_fam_",sprintf("%06d",simID),".Rds")))
  
#  df$KS <- log(df$KS)
#  names(df)[names(df)=="KS"] <- "logKS"
  
  lvls <- unique(c(df$IDA,df$IDB))
  df$IDA <- factor(df$IDA,levels = lvls)
  df$IDB <- factor(df$IDB,levels = lvls)

  lvls <- unique(c(df$FamIDA,df$FamIDB))
  df$FamIDA <- factor(df$FamIDA,levels = lvls)
  df$FamIDB <- factor(df$FamIDB,levels = lvls)
  
  df$TwinFam0IDA <- ifelse(grepl("Twin",df$TwinFamIDA, fixed = T),df$TwinFamIDA,0)
  df$TwinFam0IDB <- ifelse(grepl("Twin",df$TwinFamIDB, fixed = T),df$TwinFamIDB,0)
  tmp <- table(dat$FamID)
  famsize <- function(x) return(tmp[x])
  df$TwinFam0IDA <- ifelse(famsize(as.character(df$FamIDA)) <= 2,0,df$TwinFam0IDA)
  df$TwinFam0IDB <- ifelse(famsize(as.character(df$FamIDB)) <= 2,0,df$TwinFam0IDB)
  lvls <- unique(c(df$TwinFam0IDA,df$TwinFam0IDB))
  df$TwinFam0IDA <- factor(df$TwinFam0IDA,levels = lvls)
  df$TwinFam0IDB <- factor(df$TwinFam0IDB,levels = lvls)
  
  lvls <- unique(c(df$TwinFamIDA,df$TwinFamIDB))
  df$TwinFamIDA <- factor(df$TwinFamIDA,levels = lvls)
  df$TwinFamIDB <- factor(df$TwinFamIDB,levels = lvls)
  
  df$FamID_FamID <- factor(paste0(df$FamIDA,"_",df$FamIDB))
  df$TwinFamID_TwinFamID <- factor(paste0(df$TwinFamIDA,"_",df$TwinFamIDB))
  df$TwinFam0ID_TwinFam0ID <- factor(paste0(
    ifelse(df$TwinFam0IDA==0,as.character(df$TwinFam0IDB),as.character(df$TwinFam0IDA)),
    "_",
    ifelse(df$TwinFam0IDA==0,as.character(df$TwinFam0IDA),as.character(df$TwinFam0IDB))))
  
  for(sim in sims) {
  
  ##### -- SUBJECT INFO -- ##### load individual level covariates for MDMR method
  #########################################
  ##### ------ DEMOGRAPHIC DATA ----- #####
  #########################################
    if(sim==1) {
        cov_corr_family <- 0
        cov_corr_twin <- 0
        
        corr_family <- 0
        corr_twin <- 0
    }
    if(sim==2) {
        cov_corr_family <- 0
        cov_corr_twin <- 0
        
        corr_family <- 0.5
        corr_twin <- 0.8
    }
    if(sim==3) {
        cov_corr_family <- 0.5
        cov_corr_twin <- 0.8
        
        corr_family <- 0
        corr_twin <- 0
    }
    if(sim==4) {
        cov_corr_family <- 0.5
        cov_corr_twin <- 0.8
        
        corr_family <- 0.5
        corr_twin <- 0.8
    }
    if(sim==5) {
        cov_corr_family <- 1
        cov_corr_twin <- 1
        
        corr_family <- 0.5
        corr_twin <- 0.8
    }
    if(sim==6) {
        cov_corr_family <- 0.5
        cov_corr_twin <- 1
        
        corr_family <- 0.5
        corr_twin <- 0.8
    }    
    set.seed(task_ID)
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
  
    for(signal in signals) {
      dff <- df %>% filter(Weight==signal, Sim == sim)
      dfff <- dff %>% filter(FamIDA != FamIDB)
      
      IDmat <- model.matrix(~IDA-1,data=dfff)+model.matrix(~IDB-1,data=dfff)
      colnames(IDmat) <- gsub("^.{0,3}", "", colnames(IDmat))
      IDmatfull <- model.matrix(~IDA-1,data=dff)+model.matrix(~IDB-1,data=dff)
      colnames(IDmatfull) <- gsub("^.{0,3}", "", colnames(IDmatfull))
      # FamIDmat <- model.matrix(~FamIDA-1,data=dfff)+model.matrix(~FamIDB-1,data=dfff)
      # colnames(FamIDmat) <- gsub("^.{0,6}", "", colnames(FamIDmat))
      # TwinFamIDmat <- model.matrix(~TwinFamIDA-1,data=dfff)+model.matrix(~TwinFamIDB-1,data=dfff)
      # colnames(TwinFamIDmat) <- gsub("^.{0,10}", "", colnames(TwinFamIDmat))
      # TwinFam0IDmat <- model.matrix(~TwinFam0IDA-1,data=dfff)+model.matrix(~TwinFam0IDB-1,data=dfff)
      # colnames(TwinFam0IDmat) <- gsub("^.{0,11}", "", colnames(TwinFam0IDmat))

      for(distance in distances) {
        tempvars <- c("Age","IQ","Sex","Trt")
        
        # F-Test
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"F-Test"))
        dat_results[counter,paste0(rep(tempvars,times=4),rep(c("_est","_se","_t","_p"),each=length(tempvars)))] <- as.vector(summary(lm(get(distance) ~ Age + IQ + Sex + Trt,data=dfff))$coefficients[tempvars,])
        counter <- counter + 1
        
        # F-Test with ILE
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"F-Test with ILE"))
        dat_results[counter,paste0(rep(tempvars,times=4),rep(c("_est","_se","_t","_p"),each=length(tempvars)))] <- as.vector(summary(lm(get(distance) ~ Age + IQ + Sex + Trt - 1 + IDmat,data=dfff))$coefficients[tempvars,])
        counter <- counter + 1
        
        # Random Intercepts
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"Random Intercepts"))
        dat_results[counter,paste0(rep(tempvars,times=5),rep(c("_est","_se","_df","_t","_p"),each=length(tempvars)))] <- as.vector(summary(lmer(get(distance) ~ Age + IQ + Sex + Trt - 1 + IDmat + (1|FamID_FamID) + (1|TwinFam0ID_TwinFam0ID), data=dfff))$coefficients[tempvars,])
        counter <- counter + 1
        
        # Random Intercepts and Twin Slopes
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"R.I. and Twin Slopes"))
        dat_results[counter,paste0(rep(tempvars,times=5),rep(c("_est","_se","_df","_t","_p"),each=length(tempvars)))] <- as.vector(summary(lmer(get(distance) ~ Age + IQ + Sex + Trt - 1 + IDmat + (1|FamID_FamID) + (Age + IQ + Sex + Trt|TwinFam0ID_TwinFam0ID), data=dfff))$coefficients[tempvars,])
        counter <- counter + 1
        
        # Random Intercepts and Slopes
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"R.I. and Slopes"))
        dat_results[counter,paste0(rep(tempvars,times=5),rep(c("_est","_se","_df","_t","_p"),each=length(tempvars)))] <- as.vector(summary(lmer(get(distance) ~ Age + IQ + Sex + Trt - 1 + IDmat + (Age + IQ + Sex + Trt|FamID_FamID) + (Age + IQ + Sex + Trt|TwinFam0ID_TwinFam0ID), data=dfff))$coefficients[tempvars,])
        counter <- counter + 1
        
        D <- as.dist(VDM(dff[,distance]))
        
        mdmr.res <- mdmr(X = dat[,c("Age","IQ","Sex","Trt")], D = D, perm.p = TRUE,nperm=5000)
        x <- mdmr.res$pv
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"MDMR Permutation"))
        dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- x[c("Age","IQ","SexMale","TrtTreatment"),"perm.pvals"]
        counter = counter+1
        
        mixed.res <- mixed.mdmr(~ Age+IQ+Sex+Trt + (1|TwinFamID),data = dat, D = D)
        x <- mixed.res$pv
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"MDMR Mixed"))
        dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- x[c("Age","IQ","Sex","Trt")]
        counter = counter+1
        
#        mixed.res <- mixed.mdmr(~ Age+IQ+Sex+Trt + (Age+IQ+Sex+Trt|FamID),data = dat, D = D)
#        x <- mixed.res$pv
#        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),sim,signal,distance,"MDMR Mixed R.I.S."))
#        dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- x[c("Age","IQ","Sex","Trt")]
#        counter = counter+1 
      }
    }
  }   
  print(simID)
}
saveRDS(dat_results,output_file_path)