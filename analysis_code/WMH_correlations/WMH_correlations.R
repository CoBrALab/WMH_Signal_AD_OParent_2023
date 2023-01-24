# Correlations of WMH signal and volume measures between themselves
# First, in a between-measure within-region fashion
# Second, in a between-region within measure fashion

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

# Make subsets by region and by measure
ADB_subset_parc = list()

ADB_subset_parc[[1]] = subset(ADB_subset, select=c("volume_WMH", "qT2star_median_WMH", "qT1_median_WMH","T1w_median_WMH", "T2w_median_WMH",
                                                   "t1t2ratio_median_WMH", "FLAIR_median_WMH"))
colnames(ADB_subset_parc[[1]]) = c("Volume", "qT2star", "qT1", "T1w", "T2w", "T1w/T2w", "FLAIR")
ADB_subset_parc[[2]] = subset(ADB_subset, select=c("volume_WMH_PV", "qT2star_median_WMH_PV", "qT1_median_WMH_PV", "T1w_median_WMH_PV", 
                                                   "T2w_median_WMH_PV", "t1t2ratio_median_WMH_PV", "FLAIR_median_WMH_PV"))
colnames(ADB_subset_parc[[2]]) = c("Volume", "qT2star", "qT1", "T1w", "T2w", "T1w/T2w", "FLAIR")
ADB_subset_parc[[3]] = subset(ADB_subset, select=c("volume_WMH_deep", "qT2star_median_WMH_deep", "qT1_median_WMH_deep", "T1w_median_WMH_deep", 
                                                   "T2w_median_WMH_deep", "t1t2ratio_median_WMH_deep", "FLAIR_median_WMH_deep"))
colnames(ADB_subset_parc[[3]]) = c("Volume", "qT2star", "qT1", "T1w", "T2w", "T1w/T2w", "FLAIR")
ADB_subset_parc[[4]] = subset(ADB_subset, select=c("volume_WMH_SWM", "qT2star_median_WMH_SWM", "qT1_median_WMH_SWM", "T1w_median_WMH_SWM", 
                                                   "T2w_median_WMH_SWM", "t1t2ratio_median_WMH_SWM", "FLAIR_median_WMH_SWM"))
colnames(ADB_subset_parc[[4]]) = c("Volume", "qT2star", "qT1", "T1w", "T2w", "T1w/T2w", "FLAIR")
ADB_subset_parc[[5]] = ADB_subset[,grepl("volume_", names(ADB_subset))]
colnames(ADB_subset_parc[[5]]) = c("Global", "PV", "Deep", "SWM")
ADB_subset_parc[[6]] = ADB_subset[,grepl("qT1_", names(ADB_subset))]
colnames(ADB_subset_parc[[6]]) = c("Global", "PV", "Deep", "SWM")
ADB_subset_parc[[7]] = ADB_subset[,grepl("qT2star_", names(ADB_subset))]
colnames(ADB_subset_parc[[7]]) = c("Global", "PV", "Deep", "SWM")
ADB_subset_parc[[8]] = ADB_subset[,grepl("t1t2ratio_", names(ADB_subset))]
colnames(ADB_subset_parc[[8]]) = c("Global", "PV", "Deep", "SWM")
ADB_subset_parc[[9]] = ADB_subset[,grepl("T1w_", names(ADB_subset))]
colnames(ADB_subset_parc[[9]]) = c("Global", "PV", "Deep", "SWM")
ADB_subset_parc[[10]] = ADB_subset[,grepl("T2w_", names(ADB_subset))]
colnames(ADB_subset_parc[[10]]) = c("Global", "PV", "Deep", "SWM")
ADB_subset_parc[[11]] = ADB_subset[,grepl("FLAIR_", names(ADB_subset))]
colnames(ADB_subset_parc[[11]]) = c("Global", "PV", "Deep", "SWM")

names_parc = c("global_medians", "PV_medians", "deep_medians", "SWM_medians",
               "volume", "T1", "T2star", "t1t2ratio", "T1w", "T2w", "FLAIR")

# Correlation matrices

# Make color palette
mypalette = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                   "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                   "#4393C3", "#2166AC", "#053061")))

library("Hmisc")
library(corrplot)

dir.create("./visualization", showWarnings = FALSE)
dir.create("./results", showWarnings = FALSE)

# Between-measure within-region associations
for (i in seq(1, 4)){
  print(names_parc[i])

  # Compute correlation matrix
  corr_matrix = rcorr(as.matrix(ADB_subset_parc[[i]]), type="pearson")
  
  # Plot
  png(height=1200, width=1800, pointsize=37, file=paste0("./visualization/bmeasure_", names_parc[i],".png"))
  corrplot(as.matrix(corr_matrix$r), type="upper", order="original", method="circle", addCoef.col = "black", addgrid.col = TRUE,
           tl.srt = 45, tl.col = "black", col=mypalette(200), diag=FALSE, tl.cex=1.5, cl.length = 5, cl.cex = 1.5, number.cex=1.3,
           p.mat=corr_matrix$P, sig.level=0.01, insig="blank")
  dev.off()
  
  # Save results
  write.csv(corr_matrix$r, paste0("./results/bmeasure_",names_parc[i],"_r.csv"), row.names = TRUE)
  write.csv(corr_matrix$P, paste0("./results/bmeasure_",names_parc[i],"_p.csv"), row.names = TRUE)
}

# Between-region within-measure associations
for (i in seq(5, 11)){
  print(names_parc[i])
  
  # Compute correlation matrix
  corr_matrix = rcorr(as.matrix(ADB_subset_parc[[i]]), type="pearson")
  
  # Plot
  png(height=1200, width=1800, pointsize=37, file=paste0("./visualization/bregion_", names_parc[i],".png"))
  corrplot(as.matrix(corr_matrix$r), type="upper", order="original", method="circle", addCoef.col = "black", addgrid.col = TRUE,
           tl.srt = 45, tl.col = "black", col=mypalette(200), diag=FALSE, tl.cex=2.5, cl.length = 5, cl.cex = 2, number.cex=2.5,
           cl.ratio = 0.3, p.mat=as.matrix(corr_matrix$P), sig.level=0.01, insig="blank")
  dev.off()
  
  # Save results
  write.csv(corr_matrix$r, paste0("./results/bregion_",names_parc[i],"_r.csv"), row.names = TRUE)
  write.csv(corr_matrix$P, paste0("./results/bregion_",names_parc[i],"_p.csv"), row.names = TRUE)
}
