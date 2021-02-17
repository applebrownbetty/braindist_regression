.libPaths("/nas/longleaf/home/chalmer/bin/R3.6.0")
library(lme4)
library(lmerTest)

task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#task_ID <- 1

simIDs <- (task_ID*10+1):(task_ID*10 + 10)
sims <- 1:3
signals <- seq(0,100,10)
distances <- c("Jaccard","KS","Minkowski_0.5","Canberra","Minkowski_1.4","Euclidean","Minkowski_3")

dat_results <- data.frame(ID = rep(simIDs,each=2*length(sims)*length(signals)*length(distances)))

output_file_name <- paste0("perm_lme_",simIDs[1],"_",simIDs[length(simIDs)],".Rds")
output_file_path <- file.path("/pine/scr/c/h/chalmer/aim2/simulations/permutations/",output_file_name)

counter = 1
for(simID in simIDs) {
  df <- readRDS(file.path("/pine/scr/c/h/chalmer/aim2/simulations/datasets",paste0("sim_",sprintf("%06d",simID),".Rds")))
  df$ID11 <- df$ID1
  df$ID22 <- df$ID2
  df$ID11[seq(1,4950,2)] <- df$ID2[seq(1,4950,2)]
  df$ID22[seq(1,4950,2)] <- df$ID1[seq(1,4950,2)]
  for(sim in sims) {
    for(signal in signals) {
    #print(counter)
      for(distance in distances) {
        DM <- paste0("Sim",sim,"_wt",sprintf("%03d",signal),ifelse(distance == "Jaccard","_top20_","_corr_"),distance)

        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",simID),sim,signal,distance,"lme_unbal"))
        dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- t(summary(lmer(get(DM) ~ Age + IQ + Sex + Trt - 1 + (1|ID1) + (1|ID2), data=df))$coefficients[,"Pr(>|t|)"])
        counter <- counter + 1

        dat_results[counter,c("ID","Sim","Signal_Percent","Distance","Method")] <- t(c(sprintf("%06d",simID),sim,signal,distance,"lme_bal"))
        dat_results[counter,paste0(c("Age","IQ","Sex","Trt"),"_p")] <- t(summary(lmer(get(DM) ~ Age + IQ + Sex + Trt - 1 + (1|ID11) + (1|ID22), data=df))$coefficients[,"Pr(>|t|)"])
        counter <- counter + 1
      }
    }
  }
}
saveRDS(dat_results,output_file_path)
