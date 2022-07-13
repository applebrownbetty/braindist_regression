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
files = sub(".Rds","",sub("/pine/scr/c/h/chalmer/aim2/simulations/datasets//aim2mixed_sim_", "", files))

task_ID <- as.numeric(files[as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))])

print(task_ID)
simIDs <- task_ID
signals <- seq(0,100,10)
distances <- c("JD_top005","JD_top20","logKS","EUC","PCD","LERM")

dat_results <- data.frame(ID = rep(simIDs,each=5*length(signals)*length(distances)))

output_file_name <- paste0("aim2mixed_results_",simIDs[1],".Rds")
output_file_path <- file.path("/pine/scr/c/h/chalmer/aim2/simulations/results/",output_file_name)

##### -- SUBJECT INFO -- ##### load individual level covariates for MDMR method
set.seed(task_ID)
np <- 30 # 30 total participants (make this an even number)
dat <- tibble(
  ID = as.factor(1:np),
  Sex = ifelse(rbinom(np,1,.5)==1,"Male","Female"),
  Trt = ifelse(rbinom(np,1,.5)==1,"Treatment","Placebo"),
  IQ = round(rnorm(np,mean=100,sd=15)),
  Age = round(rnorm(np,mean=100,sd=15))
)
dat <- dat[rep(seq_len(nrow(dat)), each = 4), ]

counter = 1
for(simID in simIDs) {
  df <- readRDS(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets/",paste0("aim2mixed_sim_",sprintf("%06d",simID),".Rds")))
  df[,c("Sex","Trt")] <- 1* df[,c("Sex","Trt")] # fix from sim
  df[,c("Age","IQ")] <- abs(df[,c("Age","IQ")]) # fix from sim
  df$KS <- log(df$KS)
  names(df)[names(df)=="KS"] <- "logKS"
  df$Task <- factor(df$Task_IDA)
  lvls <- unique(c(df$IDA,df$IDB))
  df$IDA <- factor(df$IDA,levels = lvls)
  df$IDB <- factor(df$IDB,levels = lvls)
  df$ID_Task_A <- paste0("ID",df$IDA,"_Task",df$Task_IDA)
  df$ID_Task_B <- paste0("ID",df$IDB,"_Task",df$Task_IDB)
  lvls <- unique(c(df$ID_Task_A,df$ID_Task_B))
  df$ID_Task_A <- factor(df$ID_Task_A, levels = lvls)
  df$ID_Task_B <- factor(df$ID_Task_B, levels = lvls)
  df$ID_Group_A <- paste0("ID",df$IDA,"_Group",df$Group_IDA)
  df$ID_Group_B <- paste0("ID",df$IDB,"_Group",df$Group_IDB)
  lvls <- unique(c(df$ID_Group_A,df$ID_Group_B))
  df$ID_Group_A <- factor(df$ID_Group_A, levels = lvls)
  df$ID_Group_B <- factor(df$ID_Group_B, levels = lvls)
  df$ID_Task_Group_A <- paste0("ID",df$IDA,"_Task",df$Task_IDA,"_Group",df$Group_IDA)
  df$ID_Task_Group_B <- paste0("ID",df$IDB,"_Task",df$Task_IDB,"_Group",df$Group_IDB)
  lvls <- unique(c(df$ID_Task_Group_A,df$ID_Task_Group_B))
  df$ID_Task_Group_A <- factor(df$ID_Task_Group_A, levels = lvls)
  df$ID_Task_Group_B <- factor(df$ID_Task_Group_B, levels = lvls)
  for(signal in signals) {
    dff <- df %>% filter(Weight==signal)
    dfff <- dff %>% filter(IDA != IDB)
    IDmat <- model.matrix(~IDA-1,data=dfff)+model.matrix(~IDB-1,data=dfff)
    colnames(IDmat) <- paste0("ID",gsub("^.{0,3}", "", colnames(IDmat)))
    IDTaskmat <- model.matrix(~ID_Task_A-1,data=dfff)+model.matrix(~ID_Task_B-1,data=dfff)
    colnames(IDTaskmat) <- gsub("^.{0,9}", "", colnames(IDTaskmat))
    IDGroupmat <- model.matrix(~ID_Group_A-1,data=dfff)+model.matrix(~ID_Group_B-1,data=dfff)
    colnames(IDGroupmat) <- gsub("^.{0,10}", "", colnames(IDGroupmat))
    IDTaskGroupmat <- model.matrix(~ID_Task_Group_A-1,data=dfff)+model.matrix(~ID_Task_Group_B-1,data=dfff)
    colnames(IDTaskGroupmat) <- gsub("^.{0,15}", "", colnames(IDTaskGroupmat))
    dfff$ID_ID <- factor(paste0("ID",dfff$IDA,"_ID",dfff$IDB))
    dfff$ID_ID_Group <- factor(paste0("ID",dfff$IDA,"_ID",dfff$IDB,ifelse(dfff$Group_IDA == dfff$Group_IDB,"_GroupSame","_GroupDiff")))
    for(distance in distances) {
    
      tempvars <- c("Age:Task1","Age:Task2","Age:Task3","IQ:Task1","IQ:Task2","IQ:Task3","Sex:Task1","Sex:Task2","Sex:Task3","Trt:Task1","Trt:Task2","Trt:Task3")
      
      # F-Test
      dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),simID,signal,distance,"F-Test"))
      dat_results[counter,paste0(rep(tempvars,times=4),rep(c("_est","_se","_t","_p"),each=length(tempvars)))] <- as.vector(summary(lm(get(distance) ~ (Age + IQ + Sex + Trt):Task + Task - 1,data=dfff))$coefficients[tempvars,])
      counter <- counter + 1
      
      # F-Test with ILE
      dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),simID,signal,distance,"F-Test with ILE"))
      dat_results[counter,paste0(rep(tempvars,times=4),rep(c("_est","_se","_t","_p"),each=length(tempvars)))] <- as.vector(summary(lm(get(distance) ~ (Age + IQ + Sex + Trt):Task - 1 + IDTaskGroupmat,data=dfff))$coefficients[tempvars,])
      counter <- counter + 1
      
      # Mixed
      dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),simID,signal,distance,"Mixed Effects"))
      dat_results[counter,paste0(rep(tempvars,times=5),rep(c("_est","_se","_df","_t","_p"),each=length(tempvars)))] <- as.vector(summary(lmer(get(distance) ~ (Age + IQ + Sex + Trt):Task - 1 + IDTaskGroupmat + (Task|ID_ID), data=dfff))$coefficients[tempvars,])
      counter <- counter + 1
      
      for(task in 1:3) {
        dfff_mdmr <- dff %>% filter(Task == task)
        D <- as.dist(VDM(dfff_mdmr[,distance]))
        
        mdmr.res <- mdmr(X = dat[,c("Age","IQ","Sex","Trt")], D = D, perm.p = TRUE,nperm=5000)
        x <- mdmr.res$pv
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),simID,signal,distance,"mdmr_perm"))
        dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),":Task",task,"_p")] <- x[c("Age","IQ","SexMale","TrtTreatment"),"perm.pvals"]

        mixed.res <- mixed.mdmr(~ Age+IQ+Sex+Trt + (1 | ID),data = dat, D = D)
        x <- mixed.res$pv
        dat_results[counter+1,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",task_ID),simID,signal,distance,"mdmr_mixed"))
        dat_results[counter+1,paste0(c("Age","IQ","Sex","Trt"),":Task",task,"_p")] <- x[c("Age","IQ","Sex","Trt")]
      }
      counter = counter+2
    }
  }   
  print(simID)
}
saveRDS(dat_results,output_file_path)