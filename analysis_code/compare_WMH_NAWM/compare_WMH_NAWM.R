
# Analyzing if WMH signal measure slopes are significantly different between WMH and NAWM
# relative to age and WMH volume, which we interpret as indicators of 
# time and vascular burden, respectively. We used linear models with 
# interaction terms. 

# If the interaction is significant, our interpretation is that 
# the WMH signal measure is sensitive to tissue deterioration within the lesion above
# the deterioration of the global white matter

# Parcellation: Periventricular / deep / superficial white matter

ADB_subset = read.csv("../../df_pvdeepswm.csv")

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
library(viridis)

# To assign y_category later
qT1_labels = c("qT1_median", "qT1_median_PV","qT1_median_deep","qT1_median_SWM","qT1_median","qT1_median_PV","qT1_median_deep","qT1_median_SWM")
qT2star_labels = c("qT2star_median", "qT2star_median_PV","qT2star_median_deep","qT2star_median_SWM","qT2star_median","qT2star_median_PV","qT2star_median_deep","qT2star_median_SWM")
t1t2ratio_labels = c("t1t2ratio_median", "t1t2ratio_median_PV","t1t2ratio_median_deep","t1t2ratio_median_SWM","t1t2ratio_median","t1t2ratio_median_PV","t1t2ratio_median_deep","t1t2ratio_median_SWM")
T1w_labels = c("T1w_median", "T1w_median_PV","T1w_median_deep","T1w_median_SWM","T1w_median","T1w_median_PV","T1w_median_deep","T1w_median_SWM")
T2w_labels = c("T2w_median", "T2w_median_PV","T2w_median_deep","T2w_median_SWM","T2w_median","T2w_median_PV","T2w_median_deep","T2w_median_SWM")
FLAIR_labels = c("FLAIR_median", "FLAIR_median_PV","FLAIR_median_deep","FLAIR_median_SWM","FLAIR_median","FLAIR_median_PV","FLAIR_median_deep","FLAIR_median_SWM")

# Select x variables and y variables
x_vars = c(3, 52) # Age and WMH volume
colnames(ADB_subset_WMH_NAWM)[x_vars]
y_vars = c(seq(24, 29), seq(31,36), seq(38,43), seq(45,50)) # WMH signal measures
colnames(ADB_subset_WMH_NAWM)[y_vars]
color_var = ADB_subset_WMH_NAWM$WM_type

# Initialize results storage
lm_values = as.data.frame(matrix(ncol=5, nrow=0))
colnames(lm_values) = c("x", "y", "beta_interaction", "pval_interaction", "y_category")

# Run linear models and plot results
for (xvar in x_vars){
  for (yvar in y_vars){
    # Create results folder
    results_folder = "./visualization/"
    dir.create(file.path(results_folder), showWarnings = FALSE)
    cat(paste0("\n-----------\ny = ", colnames(ADB_subset_WMH_NAWM)[yvar], "\nx = ",colnames(ADB_subset_WMH_NAWM)[xvar],"\ncolor = WM_type\n-----------\n"))
    
    # Run linear models and store results
    lin_mod_interaction = lm(ADB_subset_WMH_NAWM[,yvar] ~ ADB_subset_WMH_NAWM[,xvar] * ADB_subset_WMH_NAWM$WM_type)
    beta_interaction = summary(lin_mod_interaction)$coefficients[4,1]
    pval_interaction = summary(lin_mod_interaction)$coefficients[4,4]
    
    # Assign y_category
    if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% qT1_labels) {y_categ = "qT1"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% qT2star_labels) {y_categ = "qT2star"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% t1t2ratio_labels) {y_categ = "T1w/T2w"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% T1w_labels) {y_categ = "T1w"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% T2w_labels) {y_categ = "T2w"}
    else if (colnames(ADB_subset_WMH_NAWM)[yvar] %in% FLAIR_labels) {y_categ = "FLAIR"}
   
   # Store results
    lm_values[nrow(lm_values)+1,] = c(colnames(ADB_subset_WMH_NAWM)[xvar], colnames(ADB_subset_WMH_NAWM)[yvar], beta_interaction, pval_interaction, y_categ)
    
    # Make plots
    ggplot(ADB_subset_WMH_NAWM, aes(x=ADB_subset_WMH_NAWM[,xvar], y=ADB_subset_WMH_NAWM[,yvar], fill=WM_type, color=WM_type)) +
      geom_point() +
      geom_smooth(method=lm) +
      ggtitle(colnames(ADB_subset_WMH_NAWM)[yvar]) + 
      theme(text=element_text(size=30), plot.title=element_text(size = 30, hjust=0.5),
            axis.title.x=element_blank(), axis.title.y = element_blank(),
            legend.position="none")
      ggsave(paste0(results_folder, "x.", colnames(ADB_subset_WMH_NAWM)[xvar], ".y.", colnames(ADB_subset_WMH_NAWM)[yvar], ".color.WM_type.jpg"))
  }
}

# FDR correction (only for interaction term)
lm_values$pval_interaction_fdr = p.adjust(lm_values$pval_interaction, method="fdr")
lm_values$pval_interaction = as.numeric(lm_values$pval_interaction)

# Add column with only significant pvalues
lm_values$pval_interaction_sig = NA
lm_values$pval_interaction_fdr_sig = NA
for (i in 1:nrow(lm_values)){
  if (lm_values$pval_interaction[i] <= 0.01) {lm_values$pval_interaction_sig[i] = as.numeric(lm_values$pval_interaction[i])}
  if (lm_values$pval_interaction_fdr[i] <= 0.05) {lm_values$pval_interaction_fdr_sig[i] = as.numeric(lm_values$pval_interaction_fdr[i])}
}

# Save results to csv
dir.create(file.path("./results"), showWarnings = FALSE)
write.csv(lm_values, "./results/lm_values.csv", row.names = FALSE)

# Heat maps of associations
library(grid)

lm_values$y_category = as.factor(lm_values$y_category)
lm_values$x = factor(lm_values$x, levels=c("Age", "volume_WMH_true"))
lm_values$y = as.factor(lm_values$y)

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
  ggsave("./visualization/all_lm_pval.jpg", height=8, width=6)

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
  ggsave("./visualization/all_lm_pval_fdr.jpg", height=8, width=6)

  