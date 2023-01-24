# Linear relationships between WMH signal measures and cortical atrophy, 
# cognition, demographics, group differences, and cardiovascular risk factors
# covarying for age, sex, and education, AND corresponding WMH volume
# Eg. for qT1 in deep WMH, covarying for deep WMH

# Parcellation: Periventricular / deep / superficial white matter

ADB_subset = read.csv("../../df_pvdeepswm.csv")
ADB_subset = ADB_subset[,-grep("NAWM", colnames(ADB_subset))]
ADB_subset = ADB_subset[-c(60,61,62,63)]

# Assign column types
cat_vars = c(1, 2, 4, 6, 16, 17, 18, 19, 20)
colnames(ADB_subset)[cat_vars]
colnames(ADB_subset)[-cat_vars]
ADB_subset[,cat_vars] = lapply(ADB_subset[,cat_vars], as.factor)
ADB_subset[,-cat_vars] = lapply(ADB_subset[,-cat_vars], as.numeric)

# Make graphs
library(ggplot2)
library(viridis)

# To assign x_category later
demo = c("Group", "BMI", "APOE4_status")
cog = c("MoCA_total", "RBANS_immediate_memory", "RBANS_visuospatial_memory", 
        "RBANS_language", "RBANS_attention", "RBANS_delayed_memory",
        "RBANS_total", "AD8")
ct_weights = c("CT_comp1", "CT_comp2", "CT_comp3", "CT_comp4", "CT_comp5", "CT_comp6", 
               "CT_comp7", "CT_comp8")
risk_factors = c("High_BP", "High_chol", "Diabetes", "Alcohol", "Smoking")

# Order dataset and assign x and y variables
x_vars = c(seq(2, 28), 57, 58, 59) # can be continous and categorical
ADB_subset_x = ADB_subset[,x_vars]
y_vars = seq(29, 56) # have to be continous variables
ADB_subset_y = ADB_subset[,y_vars]
ADB_subset_y = ADB_subset_y[,order(colnames(ADB_subset_y))] # Order ADB_subset_y columns by alphabetical order
ADB_subset = cbind(ADB_subset$anon_ID, ADB_subset_x, ADB_subset_y)

# Initialize results storage df
lm_values = as.data.frame(matrix(ncol=9, nrow=0))
colnames(lm_values) = c("x", "y", "std_beta", "95CI_min", "95CI_max", "Fstat", "pval", "x_category", "y_category")

# Reassign x and y variables
x_vars = seq(2, 31)
x_vars = x_vars[!x_vars %in% c(3, 4, 7)] # remove age, sex and edu
y_vars= seq(32, ncol(ADB_subset))
y_vars = y_vars[!y_vars %in% seq(56, 59)] # Remove WMH volumes

# Relate every x variable to every y variable
for (xvar in x_vars){
  for (yvar in y_vars){
    
    # If categorical x variable, do ANOVA
    if (is.factor(ADB_subset[,xvar])==TRUE){
      cat(paste0("\n-----------\ny = ", colnames(ADB_subset)[yvar], "\nx = ", colnames(ADB_subset)[xvar], "\ncolor = Age\n-----------\n"))
      
      # Assign y_category
      if (grepl("T1_", colnames(ADB_subset)[yvar])) {y_categ = "T1"}
      else if (grepl("T2star_", colnames(ADB_subset)[yvar])) {y_categ = "T2star"}
      else if (grepl("t1t2ratio_", colnames(ADB_subset)[yvar])) {y_categ = "T1w/T2w"}
      else if (grepl("T1w_", colnames(ADB_subset)[yvar])) {y_categ = "T1w"}
      else if (grepl("T2w_", colnames(ADB_subset)[yvar])) {y_categ = "T2w"}
      else if (grepl("FLAIR_", colnames(ADB_subset)[yvar])) {y_categ = "FLAIR"}
      else {y_categ = "Volume"}
      
      # Assign x_category
      if (colnames(ADB_subset)[xvar] %in% demo) {x_categ = "Demo"}
      else if (colnames(ADB_subset)[xvar] %in% cog) {x_categ = "Cognition"}
      else if (colnames(ADB_subset)[xvar] %in% ct_weights) {x_categ = "Atrophy"}
      else if (colnames(ADB_subset)[xvar] %in% risk_factors) {x_categ = "Risk Factors"}
      
      # Run linear model with corresponding WMH volume as covariate
      
      # If volume, don't add WMH volume as covariate
      if (y_categ == "Volume"){lin_mod = summary(aov(ADB_subset[,yvar] ~ ADB_subset[,xvar] + ADB_subset$Age + ADB_subset$Sex))}
      # If signal measure, add corresponding WMH volume as covariate
      else {
        if (grepl("PV", colnames(ADB_subset)[yvar], fixed=TRUE)==TRUE){
          lin_mod = summary(aov(ADB_subset[,yvar] ~ ADB_subset[,xvar] + ADB_subset$volume_WMH_PV + ADB_subset$Age + ADB_subset$Sex))
          cat("PV\n")
        }
        else if (grepl("deep", colnames(ADB_subset)[yvar], fixed=TRUE)==TRUE){
          lin_mod = summary(aov(ADB_subset[,yvar] ~ ADB_subset[,xvar] + ADB_subset$volume_WMH_deep + ADB_subset$Age + ADB_subset$Sex))
          cat("deep\n")
        }
        else if (grepl("SWM", colnames(ADB_subset)[yvar], fixed=TRUE)==TRUE){
          lin_mod = summary(aov(ADB_subset[,yvar] ~ ADB_subset[,xvar] + ADB_subset$volume_WMH_SWM + ADB_subset$Age + ADB_subset$Sex))
          cat("SWM\n")
        }
        else {
          lin_mod = summary(aov(ADB_subset[,yvar] ~ ADB_subset[,xvar] + ADB_subset$volume_WMH + ADB_subset$Age + ADB_subset$Sex))
          cat("global\n")
        }
      }
      
      fstat = lin_mod[[1]][["F value"]][1]
      pval = lin_mod[[1]][["Pr(>F)"]][1]
      
      # Store results
      lm_values[nrow(lm_values)+1,] = c(colnames(ADB_subset)[xvar], colnames(ADB_subset)[yvar], "", "", "", fstat, pval, x_categ, y_categ)
      }
    
    # If categorical x variable, do regression
    else {
      cat(paste0("\n-----------\ny = ", colnames(ADB_subset)[yvar], "\nx = ", colnames(ADB_subset)[xvar], "\ncolor = Age\n-----------\n"))
      
      # Assign y_category
      if (grepl("T1_", colnames(ADB_subset)[yvar])) {y_categ = "T1"}
      else if (grepl("T2star_", colnames(ADB_subset)[yvar])) {y_categ = "T2star"}
      else if (grepl("t1t2ratio_", colnames(ADB_subset)[yvar])) {y_categ = "T1w/T2w"}
      else if (grepl("T1w_", colnames(ADB_subset)[yvar])) {y_categ = "T1w"}
      else if (grepl("T2w_", colnames(ADB_subset)[yvar])) {y_categ = "T2w"}
      else if (grepl("FLAIR_", colnames(ADB_subset)[yvar])) {y_categ = "FLAIR"}
      else {y_categ = "Volume"}
      
      # Assign x_category
      if (colnames(ADB_subset)[xvar] %in% demo) {x_categ = "Demo"}
      else if (colnames(ADB_subset)[xvar] %in% cog) {x_categ = "Cognition"}
      else if (colnames(ADB_subset)[xvar] %in% ct_weights) {x_categ = "Atrophy"}
      else if (colnames(ADB_subset)[xvar] %in% risk_factors) {x_categ = "Risk Factors"}
      
      # If volume, don't add WMH volume as covariate
      if (y_categ == "Volume"){
        lin_mod = lm(scale(ADB_subset[,yvar]) ~ scale(ADB_subset[,xvar]) + scale(ADB_subset$Age) + ADB_subset$Sex + scale(ADB_subset$Years_school))
      }
      # If signal measure, add corresponding WMH volume as covariate
      else {
        if (grepl("PV", colnames(ADB_subset)[yvar], fixed=TRUE)==TRUE){
          lin_mod = lm(scale(ADB_subset[,yvar]) ~ scale(ADB_subset[,xvar]) + scale(ADB_subset$Age) + ADB_subset$Sex + scale(ADB_subset$Years_school) + scale(ADB_subset$volume_WMH_PV))
          cat("PV\n")
        }
        else if (grepl("deep", colnames(ADB_subset)[yvar], fixed=TRUE)==TRUE){
          lin_mod = lm(scale(ADB_subset[,yvar]) ~ scale(ADB_subset[,xvar]) + scale(ADB_subset$Age) + ADB_subset$Sex + scale(ADB_subset$Years_school) + scale(ADB_subset$volume_WMH_deep))
          cat("deep\n")
        }
        else if (grepl("SWM", colnames(ADB_subset)[yvar], fixed=TRUE)==TRUE){
          lin_mod = lm(scale(ADB_subset[,yvar]) ~ scale(ADB_subset[,xvar]) + scale(ADB_subset$Age) + ADB_subset$Sex + scale(ADB_subset$Years_school) + scale(ADB_subset$volume_WMH_SWM))
          cat("SWM\n")
        }
        else {
          lin_mod = lm(scale(ADB_subset[,yvar]) ~ scale(ADB_subset[,xvar]) + scale(ADB_subset$Age) + ADB_subset$Sex + scale(ADB_subset$Years_school) + scale(ADB_subset$volume_WMH))
          cat("global\n")
        }
      }
      
      std_beta = summary(lin_mod)$coefficients[2,1]
      pval = summary(lin_mod)$coefficients[2,4]
      r2 = summary(lin_mod)$r.squared
      fstat = as.numeric(summary(lin_mod)$fstatistic[1])
      conf_int = confint(lin_mod, "scale(ADB_subset[, xvar])", level=0.95)
      
      # Store results
      lm_values[nrow(lm_values)+1,] = c(colnames(ADB_subset)[xvar], colnames(ADB_subset)[yvar], std_beta, conf_int[1,1], conf_int[1,2], fstat, pval, x_categ, y_categ)
    }
  }
}

# FDR correction
lm_values$pval_fdr = p.adjust(lm_values$pval, method="fdr")
lm_values$pval = as.numeric(lm_values$pval)

# Add columns with only significant pvalues, std_beta
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

# Post-hoc pair-wise group differences
library(multcomp)
library(tidyverse)
library(broom)

levels(ADB_subset$Group) = c("HC", "FAMHX", "MCI", "AD")
levels(ADB_subset$Sex) = c("Male", "Female")

results_folder = paste0("./group_diff/")
dir.create(file.path(results_folder), showWarnings = FALSE)

# For each significant group ANOVA
for (i in 1:nrow(lm_values)){
  if (lm_values$x[i] == 'Group' && lm_values$pval[i] <= 0.01){
    print(lm_values$y[i])
    
    # Tukey contrasts
    group_anova = aov(scale(ADB_subset[, grep(paste0("^",lm_values$y[i],"$"), names(ADB_subset))]) ~ ADB_subset$Group + scale(ADB_subset$Age) + ADB_subset$Sex +scale(ADB_subset$Years_school))
    ph_group = glht(group_anova, linfct=mcp('ADB_subset$Group'='Tukey'))
    print(summary(ph_group))
    
    conf_int = tidy(confint(ph_group))
    
    # Plots
    ggplot(conf_int, aes(contrast, y=estimate, ymin=conf.low, ymax=conf.high)) + 
      geom_hline(yintercept=0, linetype="11", color="grey60") + 
      geom_errorbar(width=0.1) + 
      geom_point() + 
      coord_flip() + 
      theme_classic() + 
      theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14),
            axis.title.x = element_text(size=14), axis.title.y = element_blank(),
            plot.title = element_text(size=16)) + 
      ggtitle(lm_values$y[i]) +
      ggsave(paste0(results_folder,"group_", lm_values$y[i],".png"), height=5, width=6)
  }
}

# Heat map of associations
library(grid)

dir.create(file.path("visualization"), showWarnings = FALSE)

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
  ggsave("./visualization/all_lm_pval.jpg", height=10, width=13)

ggplot(lm_values, aes(x=x, y=y, fill=pval_fdr_sig)) +
  geom_tile(color = "white", lwd = 0.1, linetype = 1) +
  scale_fill_gradientn(colors = viridis_pal(direction=-1)(9), limits=c(0, 0.10), breaks=c(0, 0.10), name="Q-value") +
  facet_grid(rows = vars(y_category), cols = vars(x_category), scales="free", space="free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        legend.text = element_text(size=12), legend.title = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        strip.text = element_text(size=12)) + 
  ggsave("./visualization/all_lm_pval_fdr.jpg", height=10, width=13)

# Write lm values to csv
dir.create(file.path("results"), showWarnings = FALSE)
write.csv(lm_values, "./results/lm_values.csv", row.names=FALSE)
