# Generate graphs of WMH characteristics by age/group
library(readxl)
library(effects)
library(RMINC)
library(SpatioTemporal)


ADB_subset = read.csv("../ADB_subset_imputed_RF_new.csv")
#ADB_subset = ADB_subset[,-grep("NAWM", colnames(ADB_subset))]

# Assign column types
cat_vars = c(1, 2, 4, 6, 16, 17, 18, 19, 20)
colnames(ADB_subset)[cat_vars]
colnames(ADB_subset)[-cat_vars]
ADB_subset[,cat_vars] = lapply(ADB_subset[,cat_vars], as.factor)
ADB_subset[,-cat_vars] = lapply(ADB_subset[,-cat_vars], as.numeric)

# Remove raw volumes
ADB_subset = ADB_subset[,!grepl("raw", names(ADB_subset))]

# Divide dataset in NAWM and WMH
ADB_subset_WMH = ADB_subset[,grepl("WMH", names(ADB_subset))]
ADB_subset_WMH$WM_type = "WMH"
colnames(ADB_subset_WMH) = gsub("_WMH", "", colnames(ADB_subset_WMH))
ADB_subset_WMH = cbind(ADB_subset[,1:22], ADB_subset_WMH)

ADB_subset_NAWM = ADB_subset[,grepl("NAWM", names(ADB_subset))]
ADB_subset_NAWM$WM_type = "NAWM"
ADB_subset_NAWM = cbind(ADB_subset[,1:22], ADB_subset_NAWM)
names(ADB_subset_NAWM) = names(ADB_subset_WMH)

ADB_subset_WMH_NAWM = rbind(ADB_subset_WMH, ADB_subset_NAWM)
ADB_subset_WMH_NAWM$WM_type = as.factor(ADB_subset_WMH_NAWM$WM_type)
new_volume_WMH = rbind(as.data.frame(ADB_subset$volume_WMH), as.data.frame(ADB_subset$volume_WMH))
names(new_volume_WMH) = c("WMH_volume")
ADB_subset_WMH_NAWM = cbind(ADB_subset_WMH_NAWM, new_volume_WMH)
ADB_subset_WMH_NAWM$WMH_volume = as.numeric(ADB_subset_WMH_NAWM$WMH_volume)

# Make graphs
library(ggplot2)
library(effects)
library(viridis)
library(gridExtra)
library(ggpubr)
library(multcomp)

# To assign y_category later
volume_labels = c("volume", "volume_PV","volume_deep","volume_SWM","volume","volume_PV","volume_deep","volume_SWM",
                  "ICV", "TBV", "TBV_ICV_ratio")
qT1_labels = c("qT1_median", "qT1_median_PV","qT1_median_deep","qT1_median_SWM","qT1_median","qT1_median_PV","qT1_median_deep","qT1_median_SWM",
              "qT1_stddev", "qT1_stddev_PV","qT1_stddev_deep","qT1_stddev_SWM","qT1_stddev","qT1_stddev_PV","qT1_stddev_deep","qT1_stddev_SWM")
qT2star_labels = c("qT2star_median", "qT2star_median_PV","qT2star_median_deep","qT2star_median_SWM","qT2star_median","qT2star_median_PV","qT2star_median_deep","qT2star_median_SWM",
                  "qT2star_stddev", "qT2star_stddev_PV","qT2star_stddev_deep","qT2star_stddev_SWM","qT2star_stddev","qT2star_stddev_PV","qT2star_stddev_deep","qT2star_stddev_SWM")
QSM_labels = c("QSM_median", "QSM_median_PV","QSM_median_deep","QSM_median_SWM","QSM_median","QSM_median_PV","QSM_median_deep","QSM_median_SWM",
                "QSM_stddev", "QSM_stddev_PV","QSM_stddev_deep","QSM_stddev_SWM","QSM_stddev","QSM_stddev_PV","QSM_stddev_deep","QSM_stddev_SWM")
t1t2ratio_labels = c("t1t2ratio_median", "t1t2ratio_median_PV","t1t2ratio_median_deep","t1t2ratio_median_SWM","t1t2ratio_median","t1t2ratio_median_PV","t1t2ratio_median_deep","t1t2ratio_median_SWM",
                      "t1t2ratio_stddev", "t1t2ratio_stddev_PV","t1t2ratio_stddev_deep","t1t2ratio_stddev_SWM","t1t2ratio_stddev","t1t2ratio_stddev_PV","t1t2ratio_stddev_deep","t1t2ratio_stddev_SWM")
T1w_labels = c("T1w_median", "T1w_median_PV","T1w_median_deep","T1w_median_SWM","T1w_median","T1w_median_PV","T1w_median_deep","T1w_median_SWM",
               "T1w_stddev", "T1w_stddev_PV","T1w_stddev_deep","T1w_stddev_SWM","T1w_stddev","T1w_stddev_PV","T1w_stddev_deep","T1w_stddev_SWM")
T2w_labels = c("T2w_median", "T2w_median_PV","T2w_median_deep","T2w_median_SWM","T2w_median","T2w_median_PV","T2w_median_deep","T2w_median_SWM",
               "T2w_stddev", "T2w_stddev_PV","T2w_stddev_deep","T2w_stddev_SWM","T2w_stddev","T2w_stddev_PV","T2w_stddev_deep","T2w_stddev_SWM")
FLAIR_labels = c("FLAIR_median", "FLAIR_median_PV","FLAIR_median_deep","FLAIR_median_SWM","FLAIR_median","FLAIR_median_PV","FLAIR_median_deep","FLAIR_median_SWM",
                 "FLAIR_stddev", "FLAIR_stddev_PV","FLAIR_stddev_deep","FLAIR_stddev_SWM","FLAIR_stddev","FLAIR_stddev_PV","FLAIR_stddev_deep","FLAIR_stddev_SWM")

# ratio_labels = c("T1w_WMH_NAWM_ratio", "T1w_WMH_PV_NAWM_PV_ratio","T1w_WMH_deep_NAWM_deep_ratio","T1w_WMH_SWM_NAWM_SWM_ratio",
#                  "T2w_WMH_NAWM_ratio", "T2w_WMH_PV_NAWM_PV_ratio","T2w_WMH_deep_NAWM_deep_ratio","T2w_WMH_SWM_NAWM_SWM_ratio",
#                  "FLAIR_WMH_NAWM_ratio", "FLAIR_WMH_PV_NAWM_PV_ratio","FLAIR_WMH_deep_NAWM_deep_ratio","FLAIR_WMH_SWM_NAWM_SWM_ratio")
#
# # To assign x_category later
# demo = c("study_gr_final", "age_testing_actual", "gender_identity", "bmi_final", "APOE4_status", "demo_years_school")
# cog = c("moca_total", "rbans_immediate_memory_index_finalscore", "rbans_visuospatial_memory_index_finalscore",
#         "rbans_language_index_finalscore", "rbans_attention_index_finalscore", "rbans_delayed_memory_index_finalscore",
#         "rbans_total_scale_finalscore", "ad8_total")

x_vars = c(3, 52) # can be continous and categorical
y_vars = seq(24, 50) # have to be continous variables
color_var = ADB_subset_WMH_NAWM$WM_type
lm_values = as.data.frame(matrix(ncol=9))
colnames(lm_values) = c("x", "y", "beta_WMH", "pval_WMH", "beta_NAWM", "pval_NAWM", "beta_interaction", "pval_interaction", "y_category")

# Color with WM_type
for (xvar in x_vars){
  for (yvar in y_vars){
    results_folder = paste0("./x.",colnames(ADB_subset_WMH_NAWM)[xvar],".color.WM_type/")
    dir.create(file.path(results_folder), showWarnings = FALSE)
    cat(paste0("\n-----------\ny = ", colnames(ADB_subset_WMH_NAWM)[yvar], "\nx = ",colnames(ADB_subset_WMH_NAWM)[xvar],"\ncolor = WM_type\n-----------\n"))
    
    lin_mod_WMH = lm(ADB_subset_WMH_NAWM[,yvar] ~ ADB_subset_WMH_NAWM[,xvar], subset = ADB_subset_WMH_NAWM$WM_type == "WMH")
    beta_WMH = summary(lin_mod_WMH)$coefficients[2,1]
    pval_WMH = summary(lin_mod_WMH)$coefficients[2,4]
    r2_WMH = summary(lin_mod_WMH)$r.squared
    fstat_WMH = as.numeric(summary(lin_mod_WMH)$fstatistic[1])
    
    lin_mod_NAWM = lm(ADB_subset_WMH_NAWM[,yvar] ~ ADB_subset_WMH_NAWM[,xvar], subset = ADB_subset_WMH_NAWM$WM_type == "NAWM")
    beta_NAWM = summary(lin_mod_NAWM)$coefficients[2,1]
    pval_NAWM = summary(lin_mod_NAWM)$coefficients[2,4]
    r2_NAWM = summary(lin_mod_NAWM)$r.squared
    fstat_NAWM = as.numeric(summary(lin_mod_NAWM)$fstatistic[1])
    
    lin_mod_interaction = lm(ADB_subset_WMH_NAWM[,yvar] ~ ADB_subset_WMH_NAWM[,xvar] * ADB_subset_WMH_NAWM$WM_type)
    beta_interaction = summary(lin_mod_interaction)$coefficients[4,1]
    pval_interaction = summary(lin_mod_interaction)$coefficients[4,4]
    
    # Assign y_category between volume, T1, T2star, t1t2ratio, QSM, ratios
    if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% volume_labels) {y_categ = "Volume"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% qT1_labels) {y_categ = "qT1"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% qT2star_labels) {y_categ = "qT2star"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% QSM_labels) {y_categ = "QSM"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% t1t2ratio_labels) {y_categ = "T1w/T2w"}
    #else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% ratio_labels) {y_categ = "Ratio"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% T1w_labels) {y_categ = "T1w"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% T2w_labels) {y_categ = "T2w"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% FLAIR_labels) {y_categ = "FLAIR"}
   
    lm_values[nrow(lm_values)+1,] = c(colnames(ADB_subset_WMH_NAWM)[xvar], colnames(ADB_subset_WMH_NAWM)[yvar], beta_WMH, pval_WMH, beta_NAWM, pval_NAWM, beta_interaction, pval_interaction, y_categ)
    
    ggplot(ADB_subset_WMH_NAWM, aes(x=ADB_subset_WMH_NAWM[,xvar], y=ADB_subset_WMH_NAWM[,yvar], fill=WM_type, color=WM_type)) +
      geom_point() +
      geom_smooth(method=lm) +
      #labs(x=colnames(ADB_subset_WMH_NAWM)[xvar], y=colnames(ADB_subset_WMH_NAWM)[yvar]) +
      #ggtitle(paste0("Interaction pval = ", round(pval_interaction, 7), "\nSlope WMH = ", round(beta_WMH, 7), "\nSlope NAWM = ", round(beta_NAWM, 7))) +
      ggtitle(colnames(ADB_subset_WMH_NAWM)[yvar]) + 
      theme(text=element_text(size=30), plot.title=element_text(size = 30, hjust=0.5),
            axis.title.x=element_blank(), axis.title.y = element_blank(),
            legend.position="none")
      ggsave(paste0(results_folder, "x.", colnames(ADB_subset_WMH_NAWM)[xvar], ".y.", colnames(ADB_subset_WMH_NAWM)[yvar], ".color.WM_type.jpg"))
  }
}

# FDR correction (only for interaction term)
lm_values = lm_values[-1,]
lm_values = subset(lm_values, (lm_values$y_category != "Volume"))
lm_values$pval_interaction_fdr = p.adjust(lm_values$pval_interaction, method="fdr")
lm_values$pval_interaction = as.numeric(lm_values$pval_interaction)

# row with only significant pvalues
lm_values$pval_interaction_sig = NA
lm_values$pval_interaction_fdr_sig = NA
for (i in 1:nrow(lm_values)){
  if (lm_values$pval_interaction[i] <= 0.01) {lm_values$pval_interaction_sig[i] = as.numeric(lm_values$pval_interaction[i])}
  if (lm_values$pval_interaction_fdr[i] <= 0.05) {lm_values$pval_interaction_fdr_sig[i] = as.numeric(lm_values$pval_interaction_fdr[i])}
}

# Heat map of associations
library(viridis)
library(grid)


# lm_values$x_category = as.factor(lm_values$x_category)
lm_values$y_category = as.factor(lm_values$y_category)
lm_values$x = factor(lm_values$x, levels=c("Age", "volume_WMH_true"))
lm_values$y = as.factor(lm_values$y)
#levels(lm_values$y_category) = c("T1w", "T2w", "T1w/T2w", "FLAIR", "qT1", "qT2star")

lm_values = lm_values[order(lm_values$y_category),]

ggplot(lm_values, aes(x=x, y=y, fill=as.numeric(pval_interaction_sig))) +
  geom_tile(color = "white", lwd = 0.1, linetype = 1) +
  scale_fill_gradientn(colors = viridis_pal(direction = -1)(9), limits=c(0, 0.01), breaks=c(0, 0.01), name="P-value") +
  facet_grid(rows = vars(y_category), scales="free", space="free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
    legend.text = element_text(size=12), legend.title = element_text(size=12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    strip.text = element_text(size=12)) +
  scale_x_discrete(labels=c("Age", "WMH volume"))
  ggsave("./all_lm_pval.jpg", height=8, width=6)

ggplot(lm_values, aes(x=x, y=y, fill=as.numeric(pval_interaction_fdr_sig))) +
  geom_tile(color = "white", lwd = 0.1, linetype = 1) +
  scale_fill_gradientn(colors = viridis_pal(direction = -1)(9), limits=c(0, 0.05), breaks=c(0, 0.05), name="Q-value") +
  facet_grid(rows = vars(y_category), scales="free", space="free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
        legend.text = element_text(size=12), legend.title = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        strip.text = element_text(size=12)) +
  scale_x_discrete(labels=c("Age", "WMH volume"))
  ggsave("./all_lm_pval_fdr.jpg", height=8, width=6)

# Concatenate graphs by region
names = c("median_global.color", "median_PV.color", "median_deep.color", "median_SWM.color")
names_1 = c("global", "PV", "deep", "SWM")

for (i in 1:length(names)){
  system(paste0("montage -density 300 -tile 3x0 -geometry +10+10 -border 10 ./x.Age.color.WM_type/*", names[i], "* graphs_age_", names_1[i],".png"))
  print(paste0("montage -density 300 -tile 3x0 -geometry +10+10 -border 10 ./x.Age.color.WM_type/*", names[i], "* graphs_age_", names_1[i],".png"))
  
  # system(paste0("montage -density 300 -tile 3x0 -geometry +10+10 -border 10 ./x.volume_WMH_true.color.WM_type/*", names[i], "* graphs_volWMH_", names_1[i],".png"))
  # print(paste0("montage -density 300 -tile 3x0 -geometry +10+10 -border 10 ./x.volume_WMH_true.color.WM_type/*", names[i], "* graphs_volWMH_", names_1[i],".png"))
  
  system(paste0("montage -density 300 -tile 3x0 -geometry +10+10 -border 10 ./x.WMH_volume.color.WM_type/*", names[i], "* graphs_volWMH_", names_1[i],".png"))
  
}


# Add WMH_volume correcting for age
# for (yvar in y_vars){
#   results_folder = paste0("./x.volume_WMH_cov_age.color.WM_type/")
#   dir.create(file.path(results_folder), showWarnings = FALSE)
#   cat(paste0("\n-----------\ny = ", colnames(ADB_subset_WMH_NAWM)[yvar], "\nx = volume_WMH_cov_age\ncolor = WM_type\n-----------\n"))
#   
#   lin_mod_WMH = lm(ADB_subset_WMH_NAWM[,yvar] ~ ADB_subset_WMH_NAWM$volume_WMH_true + ADB_subset_WMH_NAWM$Age, subset = ADB_subset_WMH_NAWM$WM_type == "WMH")
#   beta_WMH = summary(lin_mod_WMH)$coefficients[2,1]
#   pval_WMH = summary(lin_mod_WMH)$coefficients[2,4]
#   r2_WMH = summary(lin_mod_WMH)$r.squared
#   fstat_WMH = as.numeric(summary(lin_mod_WMH)$fstatistic[1])
#   
#   lin_mod_NAWM = lm(ADB_subset_WMH_NAWM[,yvar] ~ ADB_subset_WMH_NAWM$volume_WMH_true + ADB_subset_WMH_NAWM$Age, subset = ADB_subset_WMH_NAWM$WM_type == "NAWM")
#   beta_NAWM = summary(lin_mod_NAWM)$coefficients[2,1]
#   pval_NAWM = summary(lin_mod_NAWM)$coefficients[2,4]
#   r2_NAWM = summary(lin_mod_NAWM)$r.squared
#   fstat_NAWM = as.numeric(summary(lin_mod_NAWM)$fstatistic[1])
#   
#   lin_mod_interaction = lm(ADB_subset_WMH_NAWM[,yvar] ~ ADB_subset_WMH_NAWM$volume_WMH_true * ADB_subset_WMH_NAWM$WM_type + ADB_subset_WMH_NAWM$Age)
#   beta_interaction = summary(lin_mod_interaction)$coefficients[4,1]
#   pval_interaction = summary(lin_mod_interaction)$coefficients[4,4]
#   
#   # Assign y_category between volume, T1, T2star, t1t2ratio, QSM, ratios
#   if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% volume_labels) {y_categ = "Volume"}
#   else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% qT1_labels) {y_categ = "qT1"}
#   else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% qT2star_labels) {y_categ = "qT2star"}
#   else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% QSM_labels) {y_categ = "QSM"}
#   else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% t1t2ratio_labels) {y_categ = "T1w/T2w"}
#   #else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% ratio_labels) {y_categ = "Ratio"}
#   else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% T1w_labels) {y_categ = "T1w"}
#   else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% T2w_labels) {y_categ = "T2w"}
#   else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% FLAIR_labels) {y_categ = "FLAIR"}
#   #
#   # # Assign x_category between demo and cog
#   # if (colnames(ADB_subset_WMH_NAWM)[xvar] %in% demo) {x_categ = "Demographics"}
#   # else if (colnames(ADB_subset_WMH_NAWM)[xvar] %in% cog) {x_categ = "Cognition"}
#   #
#   
#   lm_values[nrow(lm_values)+1,] = c("volume_WMH_cov_age", colnames(ADB_subset_WMH_NAWM)[yvar], beta_WMH, pval_WMH, beta_NAWM, pval_NAWM, beta_interaction, pval_interaction, y_categ)
#   
#   ggplot(ADB_subset_WMH_NAWM, aes(x=ADB_subset_WMH_NAWM$volume_WMH_true, y=ADB_subset_WMH_NAWM[,yvar], fill=WM_type, color=WM_type)) +
#     geom_point() +
#     geom_smooth(method=lm) +
#     labs(x=colnames(ADB_subset_WMH_NAWM)[xvar], y=colnames(ADB_subset_WMH_NAWM)[yvar]) +
#     ggtitle(paste0("Interaction pval = ", round(pval_interaction, 7), "\nSlope WMH = ", round(beta_WMH, 7), "\nSlope NAWM = ", round(beta_NAWM, 7))) +
#     theme(text=element_text(size=20), plot.title=element_text(size=15)) +
#     ggsave(paste0(results_folder, "x.", colnames(ADB_subset_WMH_NAWM)[xvar], ".y.", colnames(ADB_subset_WMH_NAWM)[yvar], ".color.WM_type.jpg"))
# }



