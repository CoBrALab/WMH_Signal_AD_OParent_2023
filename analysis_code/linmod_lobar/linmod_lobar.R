# Generate graphs of WMH characteristics by age/group

library(readxl)
library(effects)
library(RMINC)
library(SpatioTemporal)

ADB_subset = read.csv("../ADB_subset_imputed_RF_new.csv")
ADB_subset = ADB_subset[,-grep("NAWM", colnames(ADB_subset))]

cat_vars = c(1, 2, 4, 6, 16, 17, 18, 19, 20)
colnames(ADB_subset)[cat_vars]
colnames(ADB_subset)[-cat_vars]
ADB_subset[,cat_vars] = lapply(ADB_subset[,cat_vars], as.factor)
ADB_subset[,-cat_vars] = lapply(ADB_subset[,-cat_vars], as.numeric)

# Make graphs
library(ggplot2)
library(effects)
library(viridis)
library(gridExtra)
library(ggpubr)
library(multcomp)

# To assign y_category later
volume_labels = c("volume_WMH", "volume_WMH_PV","volume_WMH_deep","volume_WMH_SWM","volume_NAWM","volume_NAWM_PV","volume_NAWM_deep","volume_NAWM_SWM",
                  "ICV", "TBV", "TBV_ICV_ratio")
qT1_labels = c("qT1_median_WMH", "qT1_median_WMH_PV","qT1_median_WMH_deep","qT1_median_WMH_SWM","qT1_median_NAWM","qT1_median_NAWM_PV","qT1_median_NAWM_deep","qT1_median_NAWM_SWM",
              "qT1_stddev_WMH", "qT1_stddev_WMH_PV","qT1_stddev_WMH_deep","qT1_stddev_WMH_SWM","qT1_stddev_NAWM","qT1_stddev_NAWM_PV","qT1_stddev_NAWM_deep","qT1_stddev_NAWM_SWM")
qT2star_labels = c("qT2star_median_WMH", "qT2star_median_WMH_PV","qT2star_median_WMH_deep","qT2star_median_WMH_SWM","qT2star_median_NAWM","qT2star_median_NAWM_PV","qT2star_median_NAWM_deep","qT2star_median_NAWM_SWM",
                  "qT2star_stddev_WMH", "qT2star_stddev_WMH_PV","qT2star_stddev_WMH_deep","qT2star_stddev_WMH_SWM","qT2star_stddev_NAWM","qT2star_stddev_NAWM_PV","qT2star_stddev_NAWM_deep","qT2star_stddev_NAWM_SWM")
QSM_labels = c("QSM_median_WMH", "QSM_median_WMH_PV","QSM_median_WMH_deep","QSM_median_WMH_SWM","QSM_median_NAWM","QSM_median_NAWM_PV","QSM_median_NAWM_deep","QSM_median_NAWM_SWM",
               "QSM_stddev_WMH", "QSM_stddev_WMH_PV","QSM_stddev_WMH_deep","QSM_stddev_WMH_SWM","QSM_stddev_NAWM","QSM_stddev_NAWM_PV","QSM_stddev_NAWM_deep","QSM_stddev_NAWM_SWM")
t1t2ratio_labels = c("t1t2ratio_median_WMH", "t1t2ratio_median_WMH_PV","t1t2ratio_median_WMH_deep","t1t2ratio_median_WMH_SWM","t1t2ratio_median_NAWM","t1t2ratio_median_NAWM_PV","t1t2ratio_median_NAWM_deep","t1t2ratio_median_NAWM_SWM",
                     "t1t2ratio_stddev_WMH", "t1t2ratio_stddev_WMH_PV","t1t2ratio_stddev_WMH_deep","t1t2ratio_stddev_WMH_SWM","t1t2ratio_stddev_NAWM","t1t2ratio_stddev_NAWM_PV","t1t2ratio_stddev_NAWM_deep","t1t2ratio_stddev_NAWM_SWM")
T1w_labels = c("T1w_median_WMH", "T1w_median_WMH_PV","T1w_median_WMH_deep","T1w_median_WMH_SWM","T1w_median_NAWM","T1w_median_NAWM_PV","T1w_median_NAWM_deep","T1w_median_NAWM_SWM",
               "T1w_stddev_WMH", "T1w_stddev_WMH_PV","T1w_stddev_WMH_deep","T1w_stddev_WMH_SWM","T1w_stddev_NAWM","T1w_stddev_NAWM_PV","T1w_stddev_NAWM_deep","T1w_stddev_NAWM_SWM")
T2w_labels = c("T2w_median_WMH", "T2w_median_WMH_PV","T2w_median_WMH_deep","T2w_median_WMH_SWM","T2w_median_NAWM","T2w_median_NAWM_PV","T2w_median_NAWM_deep","T2w_median_NAWM_SWM",
               "T2w_stddev_WMH", "T2w_stddev_WMH_PV","T2w_stddev_WMH_deep","T2w_stddev_WMH_SWM","T2w_stddev_NAWM","T2w_stddev_NAWM_PV","T2w_stddev_NAWM_deep","T2w_stddev_NAWM_SWM")
FLAIR_labels = c("FLAIR_median_WMH", "FLAIR_median_WMH_PV","FLAIR_median_WMH_deep","FLAIR_median_WMH_SWM","FLAIR_median_NAWM","FLAIR_median_NAWM_PV","FLAIR_median_NAWM_deep","FLAIR_median_NAWM_SWM",
                 "FLAIR_stddev_WMH", "FLAIR_stddev_WMH_PV","FLAIR_stddev_WMH_deep","FLAIR_stddev_WMH_SWM","FLAIR_stddev_NAWM","FLAIR_stddev_NAWM_PV","FLAIR_stddev_NAWM_deep","FLAIR_stddev_NAWM_SWM")

# ratio_labels = c("T1w_WMH_NAWM_ratio", "T1w_WMH_PV_NAWM_PV_ratio","T1w_WMH_deep_NAWM_deep_ratio","T1w_WMH_SWM_NAWM_SWM_ratio",
#                  "T2w_WMH_NAWM_ratio", "T2w_WMH_PV_NAWM_PV_ratio","T2w_WMH_deep_NAWM_deep_ratio","T2w_WMH_SWM_NAWM_SWM_ratio",
#                  "FLAIR_WMH_NAWM_ratio", "FLAIR_WMH_PV_NAWM_PV_ratio","FLAIR_WMH_deep_NAWM_deep_ratio","FLAIR_WMH_SWM_NAWM_SWM_ratio")

# To assign x_category later
demo = c("Cohort", "Group", "Age", "Sex", "BMI", "APOE4_status", "Years_school")
cog = c("MoCA_total", "RBANS_immediate_memory", "RBANS_visuospatial_memory", 
        "RBANS_language", "RBANS_attention", "RBANS_delayed_memory",
        "RBANS_total", "AD8")
ct_weights = c("CT_comp1", "CT_comp2", "CT_comp3", "CT_comp4", "CT_comp5", "CT_comp6", 
               "CT_comp7", "CT_comp8", "CT_comp9", "CT_comp10", "CT_comp11", "CT_comp12")
risk_factors = c("Concussion", "TIA", "High_BP", "High_chol", "Diabetes", "Alcohol", "Smoking")

x_vars = c(seq(2, 28),64,65,66) # can be continous and categorical
ADB_subset_x = ADB_subset[,x_vars]
y_vars = seq(29, 63) # have to be continous variables
ADB_subset_y = ADB_subset[,y_vars]
ADB_subset_y = ADB_subset_y[,order(colnames(ADB_subset_y))] # Order ADB_subset_y columns by alphabetical order
ADB_subset = cbind(ADB_subset$ID, ADB_subset_x, ADB_subset_y)

color_var = ADB_subset$Age # Has to be continous variable

lm_values = as.data.frame(matrix(ncol=9))
colnames(lm_values) = c("x", "y", "std_beta", "95CI_min", "95CI_max", "Fstat", "pval", "x_category", "y_category")

x_vars = seq(2, 31)
x_vars = x_vars[!x_vars %in% c(3,4,7)] # remove age, sex and Age from x_vars
y_vars= seq(32, ncol(ADB_subset))

# Color with age
for (xvar in x_vars){
  for (yvar in y_vars){

    # Categorical x variables (covarying for age and sex)
    if (is.factor(ADB_subset[,xvar])==TRUE){
      results_folder = paste0("./x.", colnames(ADB_subset)[xvar],".color.age/")
      #dir.create(file.path(results_folder), showWarnings = FALSE)
      cat(paste0("\n-----------\ny = ", colnames(ADB_subset)[yvar], "\nx = ", colnames(ADB_subset)[xvar], "\ncolor = Age\n-----------\n"))
      
      # Run linear model
      lin_mod = summary(aov(ADB_subset[,yvar] ~ ADB_subset[,xvar] + ADB_subset$Age + ADB_subset$Sex + ADB_subset$Years_school))
      fstat = lin_mod[[1]][["F value"]][1]
      pval = lin_mod[[1]][["Pr(>F)"]][1]
      
      # Assign y_category between volume, T1, T2star, t1t2ratio, QSM, ratios
      if (colnames(ADB_subset)[yvar] %in% volume_labels) {y_categ = "Volume"}
      else if (colnames(ADB_subset)[yvar] %in% qT1_labels) {y_categ = "qT1"}
      else if (colnames(ADB_subset)[yvar] %in% qT2star_labels) {y_categ = "qT2star"}
      else if (colnames(ADB_subset)[yvar] %in% QSM_labels) {y_categ = "QSM"}
      else if (colnames(ADB_subset)[yvar] %in% t1t2ratio_labels) {y_categ = "T1w/T2w"}
      #else if (colnames(ADB_subset)[yvar] %in% ratio_labels) {y_categ = "Ratio"}
      else if (colnames(ADB_subset)[yvar] %in% T1w_labels) {y_categ = "T1w"}
      else if (colnames(ADB_subset)[yvar] %in% T2w_labels) {y_categ = "T2w"}
      else if (colnames(ADB_subset)[yvar] %in% FLAIR_labels) {y_categ = "FLAIR"}
      
      # Assign x_category between demo and cog
      if (colnames(ADB_subset)[xvar] %in% demo) {x_categ = "Demo"}
      else if (colnames(ADB_subset)[xvar] %in% cog) {x_categ = "Cognition"}
      else if (colnames(ADB_subset)[xvar] %in% ct_weights) {x_categ = "Atrophy"}
      else if (colnames(ADB_subset)[xvar] %in% risk_factors) {x_categ = "Risk Factors"}
      
      lm_values[nrow(lm_values)+1,] = c(colnames(ADB_subset)[xvar], colnames(ADB_subset)[yvar], "", "", "", fstat, pval, x_categ, y_categ)
      
      # ggplot(ADB_subset, aes(as.factor(x=ADB_subset[,xvar]), y=ADB_subset[,yvar])) +
      #   geom_jitter(aes(color=color_var),height=0) +
      #   stat_boxplot(fill=NA) +
      #   labs(y=colnames(ADB_subset)[yvar]) +
      #   ggtitle(paste0("F value = ", round(fstat,3), "\np = ", round(pval,7))) +
      #   scale_color_viridis("Age") +
      #   scale_fill_viridis() +
      #   theme(text=element_text(size=20), plot.title=element_text(size=15)) +
      #   scale_x_discrete(name=colnames(ADB_subset)[xvar]) +
      #   #scale_x_discrete(name="Group", breaks=c("0", "1", "2", "3", "4"), labels=c("HC", "FAMHX", "SCI", "MCI", "AD")) +
      #   ggsave(paste0(results_folder, "x.", colnames(ADB_subset)[xvar], ".y.", colnames(ADB_subset)[yvar], ".color.age.jpg"))
      
    }
    
    # Continuous x variables (residualizing for age and sex)
    else {
      results_folder = paste0("./x.", colnames(ADB_subset)[xvar],".color.age/")
      #dir.create(file.path(results_folder), showWarnings = FALSE)
      cat(paste0("\n-----------\ny = ", colnames(ADB_subset)[yvar], "\nx = ", colnames(ADB_subset)[xvar], "\ncolor = Age\n-----------\n"))
      
      # residualize x variable for age and sex
      lin_mod = lm(scale(ADB_subset[,yvar]) ~ scale(ADB_subset[,xvar]) + scale(ADB_subset$Age) + ADB_subset$Sex + scale(ADB_subset$Years_school))
      
      std_beta = summary(lin_mod)$coefficients[2,1]
      pval = summary(lin_mod)$coefficients[2,4]
      r2 = summary(lin_mod)$r.squared
      fstat = as.numeric(summary(lin_mod)$fstatistic[1])
      
      conf_int = confint(lin_mod, "scale(ADB_subset[, xvar])", level=0.95)
      
      # Assign y_category between volume, T1, T2star, t1t2ratio, QSM, ratios
      if (colnames(ADB_subset)[yvar] %in% volume_labels) {y_categ = "Volume"}
      else if (colnames(ADB_subset)[yvar] %in% qT1_labels) {y_categ = "qT1"}
      else if (colnames(ADB_subset)[yvar] %in% qT2star_labels) {y_categ = "qT2star"}
      else if (colnames(ADB_subset)[yvar] %in% QSM_labels) {y_categ = "QSM"}
      else if (colnames(ADB_subset)[yvar] %in% t1t2ratio_labels) {y_categ = "T1w/T2w"}
      #else if (colnames(ADB_subset)[yvar] %in% ratio_labels) {y_categ = "Ratio"}
      else if (colnames(ADB_subset)[yvar] %in% T1w_labels) {y_categ = "T1w"}
      else if (colnames(ADB_subset)[yvar] %in% T2w_labels) {y_categ = "T2w"}
      else if (colnames(ADB_subset)[yvar] %in% FLAIR_labels) {y_categ = "FLAIR"}
      
      # Assign x_category between demo and cog
      if (colnames(ADB_subset)[xvar] %in% demo) {x_categ = "Demo"}
      else if (colnames(ADB_subset)[xvar] %in% cog) {x_categ = "Cognition"}
      else if (colnames(ADB_subset)[xvar] %in% ct_weights) {x_categ = "Atrophy"}
      else if (colnames(ADB_subset)[xvar] %in% risk_factors) {x_categ = "Risk Factors"}
      
      lm_values[nrow(lm_values)+1,] = c(colnames(ADB_subset)[xvar], colnames(ADB_subset)[yvar], std_beta, conf_int[1,1], conf_int[1,2], fstat, pval, x_categ, y_categ)
      
      # ggplot(ADB_subset, aes(x=x, y=y)) +
      #   geom_point(aes(color=color_var)) +
      #   geom_smooth(aes(y=y, x=x), method=lm) +
      #   labs(x=colnames(ADB_subset)[xvar], y=colnames(ADB_subset)[yvar]) +
      #   ggtitle(paste0("Std beta = ", round(std_beta, 3), "\np = ", round(pval, 7), "\nR2 = ", round(r2, 3))) +
      #   scale_color_viridis("Age") +
      #   #scale_color_viridis(name="Group", breaks=c("0", "1", "2", "3", "4"), labels=c("HC", "FAMHX", "SCI", "MCI", "AD"), discrete=TRUE) +
      #   scale_fill_viridis() +
      #   theme(text=element_text(size=20), plot.title=element_text(size=15)) +
      #   ggsave(paste0(results_folder, "x.", colnames(ADB_subset)[xvar], ".y.", colnames(ADB_subset)[yvar], ".color.age.jpg"))
    }
  }
  
  # Combine global results
  # dir.create(file.path("./temp"), showWarnings = FALSE)
  # system(paste0("convert ", results_folder,"x.", colnames(ADB_subset)[xvar], ".y.WMH.*_median.color.age.jpg ", results_folder,"x.", colnames(ADB_subset)[xvar], ".y.WMH.volume.color.age.jpg ", "+append ./temp/WMH.jpg"))
  # system(paste0("convert ", results_folder,"x.", colnames(ADB_subset)[xvar], ".y.NAWM.*_median.color.age.jpg ", results_folder,"x.", colnames(ADB_subset)[xvar], ".y.NAWM.volume.color.age.jpg ", "+append ./temp/NAWM.jpg"))
  # system(paste0("convert ", results_folder,"x.", colnames(ADB_subset)[xvar], ".y.T*w_WMH_NAWM_ratio.*.jpg ",  results_folder,"x.", colnames(ADB_subset)[xvar], ".y.TBV_ICV_ratio.*.jpg +append ./temp/ratio.jpg"))
  # system(paste0("convert ./temp/WMH.jpg ./temp/NAWM.jpg ./temp/ratio.jpg -append ", "./x.", colnames(ADB_subset)[xvar],".y.all_global.color.age.jpg"))
}

# FDR correction
lm_values = lm_values[-1,]
lm_values$pval_fdr = p.adjust(lm_values$pval, method="fdr")
lm_values$pval = as.numeric(lm_values$pval)

# row with only significant pvalues, std_beta
lm_values$pval_sig = NA
lm_values$pval_fdr_sig = NA
lm_values$std_beta_sig = NA
lm_values$std_beta_fdr_sig = NA

for (i in 1:nrow(lm_values)){
  if (lm_values$pval[i] <= 0.01) {
    lm_values$pval_sig[i] = as.numeric(lm_values$pval[i])
    lm_values$std_beta_sig[i] = as.numeric(lm_values$std_beta[i])
  }
  if (lm_values$pval_fdr[i] <= 0.10) {
    lm_values$pval_fdr_sig[i] = as.numeric(lm_values$pval_fdr[i])
    lm_values$std_beta_fdr_sig[i] = as.numeric(lm_values$std_beta[i])
  }
}
sig_vars = sig_vars[-1,]

# Post-hoc pair-wise group differences
library(multcomp)
library(tidyverse)
library(broom)

levels(ADB_subset$Group) = c("HC", "FAMHX", "MCI", "AD")
levels(ADB_subset$Sex) = c("Male", "Female")

results_folder = paste0("./group_diff/")
dir.create(file.path(results_folder), showWarnings = FALSE)

for (i in 1:nrow(lm_values)){
  if (lm_values$x[i] == 'Group' && lm_values$pval[i] <= 0.01){
    print(lm_values$y[i])
    
    group_anova = aov(ADB_subset[, grep(paste0("^",lm_values$y[i],"$"), names(ADB_subset))] ~ ADB_subset$Group + ADB_subset$Age + ADB_subset$Sex)
    ph_group = glht(group_anova, linfct=mcp('ADB_subset$Group'='Tukey'))
    print(summary(ph_group))
    
    conf_int = tidy(confint(ph_group))
    
    ggplot(conf_int, aes(contrast, y=estimate, ymin=conf.low, ymax=conf.high)) + 
      geom_hline(yintercept=0, linetype="11", color="grey60") + 
      geom_errorbar(width=0.1) + 
      geom_point() + 
      coord_flip() + 
      theme_classic() + 
      ggtitle(lm_values$y[i]) +
      ggsave(paste0(results_folder,"group_", lm_values$y[i],".png"))
    
    # png(file=paste0(results_folder,"group_", lm_values$y[i],".png"))
    # plot(confint(ph_group))
    # dev.off()
  }
}

# Heat map of associations
library(viridis)
library(grid)
library(scales)

lm_values$x_category = as.factor(lm_values$x_category)
lm_values$y_category = as.factor(lm_values$y_category)
lm_values$x = as.factor(lm_values$x)
lm_values$y = as.factor(lm_values$y)

lm_values = lm_values[order(lm_values$y_category, lm_values$x_category),]

# Graphs for pvals
ggplot(lm_values, aes(x=x, y=y, fill=pval_sig)) +
  geom_tile(color = "white", lwd = 0.1, linetype = 1) +
  scale_fill_gradientn(colors = viridis_pal(direction=-1)(9), limits=c(0, 0.01), breaks=c(0, 0.01), name="P-value") +
  facet_grid(rows = vars(y_category), cols = vars(x_category), scales="free", space="free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        legend.text = element_text(size=12), legend.title = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        strip.text = element_text(size=12)) + 
  ggsave("./all_lm_pval.jpg", height=10, width=13)

ggplot(lm_values, aes(x=x, y=y, fill=pval_fdr_sig)) +
  geom_tile(color = "white", lwd = 0.1, linetype = 1) +
  scale_fill_gradientn(colors = viridis_pal(direction=-1)(9), limits=c(0, 0.10), breaks=c(0, 0.10), name="Q-value") +
  facet_grid(rows = vars(y_category), cols = vars(x_category), scales="free", space="free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        legend.text = element_text(size=12), legend.title = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        strip.text = element_text(size=12)) + 
  ggsave("./all_lm_pval_fdr.jpg", height=10, width=13)

#graphs for std_beta
ggplot(lm_values, aes(x=x, y=y, fill=as.numeric(std_beta))) +
  geom_tile(color = "white", lwd = 0.1, linetype = 1) +
  scale_fill_gradientn(colours=c("#0000FFFF", "#FFFFFFFF", "#FF0000FF")) +
  facet_grid(rows = vars(y_category), cols = vars(x_category), scales="free", space="free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        strip.text = element_text(size=12)) + 
  ggsave("./all_lm_std_beta.jpg", height=10, width=13)

ggplot(lm_values, aes(x=x, y=y, fill=std_beta_sig)) +
  geom_tile(color = "white", lwd = 0.1, linetype = 1) +
  scale_fill_gradientn(colours=c("#0000FFFF", "#FFFFFFFF", "#FF0000FF")) +
  facet_grid(rows = vars(y_category), cols = vars(x_category), scales="free", space="free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        strip.text = element_text(size=12)) +
  ggsave("./all_lm_std_beta_sig.jpg", height=10, width=13)

ggplot(lm_values, aes(x=x, y=y, fill=std_beta_fdr_sig)) +
  geom_tile(color = "white", lwd = 0.1, linetype = 1) +
  scale_fill_gradientn(colours=c("#0000FFFF", "#FFFFFFFF", "#FF0000FF")) +
  facet_grid(rows = vars(y_category), cols = vars(x_category), scales="free", space="free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        strip.text = element_text(size=12)) +
  ggsave("./all_lm_std_beta_sig_fdr.jpg", height=10, width=13)

# Write lm values to csv
write.csv(lm_values, "./lm_values.csv", col.names=TRUE, row.names=FALSE)
