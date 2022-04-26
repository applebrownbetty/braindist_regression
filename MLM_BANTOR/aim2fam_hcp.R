.libPaths("/nas/longleaf/home/chalmer/bin/R4.1.0")
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)

mod_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

HCP_fdf <- readRDS("HCP_famdistdf.Rds") %>% filter(FamIDa != FamIDb)
IDs <- unique(c(HCP_fdf$IDa,HCP_fdf$IDb))
HCP_dem <- as.data.frame(read_excel("Demographics_Restricted.xlsx")) %>% filter(Subject %in% IDs)
rownames(HCP_dem) <- HCP_dem$Subject
HCP_dem_unrestricted <- as.data.frame(read_csv("Demographics_Unrestricted.csv")) %>% filter(Subject %in% IDs)
rownames(HCP_dem_unrestricted) <- HCP_dem_unrestricted$Subject

df <- HCP_fdf

  lvls <- unique(c(df$IDa,df$IDb)) 
  df$IDa <- factor(df$IDa,levels = lvls)
  df$IDb <- factor(df$IDb,levels = lvls)
  IDmat <- model.matrix(~IDa-1,data=df)+model.matrix(~IDb-1,data=df)
    colnames(IDmat) <- gsub("^.{0,3}", "", colnames(IDmat))

  df$Fam0ID_Fam0ID <- factor(paste0(
    ifelse(substr(df$Fam0IDa, 1, 5) < substr(df$Fam0IDb, 1, 5),df$Fam0IDa,df$Fam0IDb),
    "_",
    ifelse(substr(df$Fam0IDa, 1, 5) < substr(df$Fam0IDb, 1, 5),df$Fam0IDb,df$Fam0IDa)))
  df$Twin0ID_Twin0ID <- factor(paste0(
    ifelse(substr(df$Twin0IDa, 1, 5) < substr(df$Twin0IDb, 1, 5),df$Twin0IDa,df$Twin0IDb),
    "_",
    ifelse(substr(df$Twin0IDa, 1, 5) < substr(df$Twin0IDb, 1, 5),df$Twin0IDb,df$Twin0IDa)))
    
if(mod_ID==1) {
mod1 <- lmer(JAC_RS ~ FluidIntl + Gender + Age + Education + BMI + Race + Ethnicity + Handedness + Income + DSM4_Alc_Abuse + SmokeStatus - 1 + IDmat + (Age + BMI + Income + SmokeStatus|Fam0ID_Fam0ID) + (1|Twin0ID_Twin0ID), data=df)
write.csv(summary(mod1)$coefficients,file="hcp_results_RS_spec_JAC.csv")
}
if(mod_ID==2) {
mod1 <- lmer(EUC_RS ~ FluidIntl + Gender + Age + Education + BMI + Race + Ethnicity + Handedness + Income + DSM4_Alc_Abuse + SmokeStatus - 1 + IDmat + (Age + BMI + Income + SmokeStatus|Fam0ID_Fam0ID) + (1|Twin0ID_Twin0ID), data=df)
write.csv(summary(mod1)$coefficients,file="hcp_results_RS_spec_EUC.csv")
}
if(mod_ID==3) {
mod1 <- lmer(LERM_RS ~ FluidIntl + Gender + Age + Education + BMI + Race + Ethnicity + Handedness + Income + DSM4_Alc_Abuse + SmokeStatus - 1 + IDmat + (Age + BMI + Income + SmokeStatus|Fam0ID_Fam0ID) + (1|Twin0ID_Twin0ID), data=df)
write.csv(summary(mod1)$coefficients,file="hcp_results_RS_spec_LERM.csv")
}
if(mod_ID==4) {
mod1 <- lmer(JAC_WM ~ FluidIntl + Gender + Age + Education + BMI + Race + Ethnicity + Handedness + Income + DSM4_Alc_Abuse + SmokeStatus - 1 + IDmat + (Age + BMI + Income + SmokeStatus|Fam0ID_Fam0ID) + (1|Twin0ID_Twin0ID), data=df)
write.csv(summary(mod1)$coefficients,file="hcp_results_WM_spec_JAC.csv")
}
if(mod_ID==5) {
mod1 <- lmer(EUC_WM ~ FluidIntl + Gender + Age + Education + BMI + Race + Ethnicity + Handedness + Income + DSM4_Alc_Abuse + SmokeStatus - 1 + IDmat + (Age + BMI + Income + SmokeStatus|Fam0ID_Fam0ID) + (1|Twin0ID_Twin0ID), data=df)
write.csv(summary(mod1)$coefficients,file="hcp_results_WM_spec_EUC.csv")
}
if(mod_ID==6) {
mod1 <- lmer(LERM_WM ~ FluidIntl + Gender + Age + Education + BMI + Race + Ethnicity + Handedness + Income + DSM4_Alc_Abuse + SmokeStatus - 1 + IDmat + (Age + BMI + Income + SmokeStatus|Fam0ID_Fam0ID) + (1|Twin0ID_Twin0ID), data=df)
write.csv(summary(mod1)$coefficients,file="hcp_results_WM_spec_LERM.csv")
}
