.libPaths("/nas/longleaf/home/chalmer/bin/R3.6.0")
library(tidyverse)
library(MDMR)

task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#task_ID = 1


output_file_name <- paste0("perm_mdmr_",sprintf("%06d",task_ID),".Rds")
output_file_path <- file.path("/pine/scr/c/h/chalmer/aim2/simulations/results/",output_file_name)

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

simIDs <- (task_ID*10+1):(task_ID*10 + 10)
sims <- 1:3
signals <- seq(0,100,10)
distances <- c("Jaccard","KS","Minkowski_0.5","Canberra","Minkowski_1.4","Euclidean","Minkowski_3")

dat_results <- data.frame(ID = rep(simIDs,each=2*(length(sims)+1)*length(signals)*length(distances)))

counter = 1

for(simID in simIDs) {
  df <- readRDS(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets",paste0("sim_",sprintf("%06d",simID),".Rds")))
  set.seed(simID)
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
  
  for(sim in sims) {
    for(signal in signals) {
      for(distance in distances) {
        DM <- paste0("Sim",sim,"_wt",sprintf("%03d",signal),ifelse(distance == "Jaccard","_top20_","_corr_"),distance)
        y <- as.matrix(df[,DM])
        if(distance=="KS") y <- log(y)

        D <- VDM(y)

        invisible(capture.output(suppressWarnings(mdmr.res <- mdmr(X = dat[,2:5], D = D,perm.p = FALSE))))
        x <- mdmr.res$pv
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",simID),sim,signal,distance,"mdmr_analytic"))
        dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- x[c("Age","IQ","SexMale","TrtTreatment"),"analytic.pvals"]
        counter = counter+1
        
        invisible(capture.output(suppressWarnings(mdmr.res <- mdmr(X = dat[,2:5], D = D,perm.p = TRUE,nperm=10000))))
        x <- mdmr.res$pv
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",simID),sim,signal,distance,"mdmr_perm"))
        dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- x[c("Age","IQ","SexMale","TrtTreatment"),"perm.pvals"]
        counter = counter+1
      }
    }
  }

  df <- readRDS(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets",paste0("sim4_",sprintf("%06d",simID),".Rds")))
  sim <- 4
  for(signal in signals) {
    for(distance in distances) {
      DM <- paste0("Sim",sim,"_wt",sprintf("%03d",signal),ifelse(distance == "Jaccard","_top20_","_corr_"),distance)
      y <- as.matrix(df[,DM])
      if(distance=="KS") y <- log(y)

      D <- VDM(y)

      invisible(capture.output(suppressWarnings(mdmr.res <- mdmr(X = dat[,2:5], D = D,perm.p = FALSE))))
      x <- mdmr.res$pv
      dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",simID),sim,signal,distance,"mdmr_analytic"))
      dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- x[c("Age","IQ","SexMale","TrtTreatment"),"analytic.pvals"]
      counter = counter+1
      
      invisible(capture.output(suppressWarnings(mdmr.res <- mdmr(X = dat[,2:5], D = D,perm.p = TRUE,nperm=10000))))
      x <- mdmr.res$pv
      dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",simID),sim,signal,distance,"mdmr_perm"))
      dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- x[c("Age","IQ","SexMale","TrtTreatment"),"perm.pvals"]
      counter = counter+1
    }
  }
}

saveRDS(dat_results,output_file_path)