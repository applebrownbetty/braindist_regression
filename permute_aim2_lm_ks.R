.libPaths("/nas/longleaf/home/chalmer/bin/R3.6.0")

task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#task_ID <- 0

df <- readRDS(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets",paste0("sim4_",sprintf("%06d",25),".Rds")))
IDmat <- model.matrix(~ID1-1,data=df)+model.matrix(~ID2-1,data=df)
colnames(IDmat) <- paste0("ID_",1:100)

simIDs <- (task_ID*100+1):(task_ID*100 + 100)
sims <- 4
signals <- seq(0,100,10)
distances <- c("Jaccard","KS","Minkowski_0.5","Canberra","Minkowski_1.4","Euclidean","Minkowski_3")

dat_results <- data.frame(ID = rep(simIDs,each=length(sims)*length(signals)*length(distances)))

output_file_name <- paste0("perm_lm_ks_",simIDs[1],"_",simIDs[length(simIDs)],".Rds")
output_file_path <- file.path("/pine/scr/c/h/chalmer/aim2/simulations/permutations/",output_file_name)

counter = 1
for(simID in simIDs) {
  df <- readRDS(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets",paste0("sim4_",sprintf("%06d",simID),".Rds")))
  for(sim in sims) {
    for(signal in signals) {
    #print(counter)
      for(distance in distances) {
        DM <- paste0("Sim",sim,"_wt",sprintf("%03d",signal),ifelse(distance == "Jaccard","_top20_","_corr_"),distance)
        y <- df[,DM]
        names(y) <- "y"
        dftmp <- cbind(y,df[,c("Age","IQ","Sex","Trt")],IDmat)
        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",simID),sim,signal,distance,"parameric_ID"))
        dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- t(summary(lm(y ~ .-1,data=dftmp))$coefficients[c("Age","IQ","Sex","Trt"),"Pr(>|t|)"])
        counter <- counter + 1
      }
    }
  }
}
saveRDS(dat_results,output_file_path)
