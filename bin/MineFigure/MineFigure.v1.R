
library(ggplot2)
library(RColorBrewer)


################################# all gene methylation rate ################################# 

MethylRate_data_curation <- function (data) {
  data <- as.data.frame(data)
  data_r_methyl <- data.frame(r_methyl = data$S1.r_methyl, sample = 'S1', color = 'S1')
  data_r_methyl <- rbind(data_r_methyl, data.frame(r_methyl = data$S2.r_methyl, sample = 'S2', color = "S2"))
  data_r_methyl <- rbind(data_r_methyl, data.frame(r_methyl = data$PIG_3.r_methyl, sample = 'PIG_3', color = "PIG_3"))
  data_r_methyl <- rbind(data_r_methyl, data.frame(r_methyl = data$PIG_4.r_methyl, sample = 'PIG_4', color = "PIG_4"))
  data_r_methyl <- rbind(data_r_methyl, data.frame(r_methyl = data$Liver_1.r_methyl, sample = 'Liver_1', color = "Liver_1"))
  data_r_methyl <- rbind(data_r_methyl, data.frame(r_methyl = data$Liver_2.r_methyl, sample = 'Liver_2', color = "Liver_2"))
  data_r_methyl <- rbind(data_r_methyl, data.frame(r_methyl = data$Liver_3.r_methyl, sample = 'Liver_3', color = "Liver_3"))
  data_r_methyl <- rbind(data_r_methyl, data.frame(r_methyl = data$Liver_26.r_methyl, sample = 'Liver_26', color = "Liver_26"))
  return(data_r_methyl)
}

MethylRate_data_draw_boxplot <- function (data) {
  
  figure <- ggplot(data = data, aes(x = sample, y = r_methyl, group = sample, col = color)) + geom_boxplot() +
            scale_x_discrete(name = "Samples") +
            scale_y_continuous(name = "Methylation Rate in all genes (%)", breaks = seq(0, 100, 10), limits=c(0, 100)) +
            theme_bw() +
            scale_color_manual(name="Samples", values=brewer.pal(8, "Set1") )
  
  return(figure)
}

mf_home <- "/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD190324/MineFigure"

context_type <- "CpG"
context_type <- "CHG"
context_type <- "CHH"

boundary_type <- "extend"
boundary_type <- "origin"
boundary_type <- "up"
boundary_type <- "down"

strand_type <- "OA" #all
strand_type <- "OB" #bottom
strand_type <- "OT" #top

data <- read.table(paste(mf_home, "AllGeneMethylRate", paste("AllGeneMethylRate", context_type, strand_type, "_methyl", boundary_type, sep = "."), sep = "/"), header = TRUE)
data <- MethylRate_data_curation(data)
(figure <- MethylRate_data_draw_boxplot(data))
ggsave(paste(mf_home, "AllGeneMethylRate", paste("AllGeneMethylRate", context_type, strand_type, boundary_type, "jpg", sep = "."), sep = "/"), dpi = 300, width = 7, height = 5)


#################################  all gene methylation rate by bin ################################# 

MethylRateByBin_data_curation <- function(data){
  bin_order <- c(1:70)
  data_mtx <- data.frame(region = bin_order, r_methyl = data$S1.r_methyl, sample = "S1", color = "S1")
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$S2.r_methyl, sample = "S2", color = "S2"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$PIG_3.r_methyl, sample = "PIG_3", color = "PIG_3"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$PIG_4.r_methyl, sample = "PIG_4", color = "PIG_4"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$Liver_1.r_methyl, sample = "Liver_1", color = "Liver_1"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$Liver_2.r_methyl, sample = "Liver_2", color = "Liver_2"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$Liver_3.r_methyl, sample = "Liver_3", color = "Liver_3"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$Liver_26.r_methyl, sample = "Liver_26", color = "Liver_26"))
  return(data_mtx)
}

MethylRateByBin_data_draw_lineplot <- function (data_mtx) {
  figure <- ggplot(data = data_mtx, aes(x = region, y = r_methyl, group = sample, col = color)) + geom_line(size=0.6) + 
            scale_y_continuous(limits=c(0, 100)) + 
            xlab('Up2k-GeneBody-Down2k') + ylab('Methylation Rate (%)') + 
            theme_bw() +
            scale_color_manual(name="Samples", values=brewer.pal(8, "Set1"))
  return(figure)
}

mf_home <- "/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD190324/MineFigure"

context_type <- "CpG"
context_type <- "CHG"
context_type <- "CHH"

strand_type <- "OA" #all
strand_type <- "OB" #bottom
strand_type <- "OT" #top

data <- read.table(paste(mf_home, "AllGeneMethylRateByBin", paste("AllGeneMethylRateByBin", context_type, strand_type, "_methylByBinSum", sep = "."), sep = "/"), header = TRUE)
data <- MethylRateByBin_data_curation(data)
(figure <- MethylRateByBin_data_draw_lineplot(data))
ggsave(paste(mf_home, "AllGeneMethylRateByBin", paste("AllGeneMethylRateByBin", context_type, strand_type, "_methylByBinSum.jpg", sep = "."), sep = "/"), dpi = 300, width = 10, height = 3)


################################# all gene DMC ################################# 

parse_DMC_ByBinSum_for_d_bin_dmc <- function(data) {
  data <- as.data.frame(data)
  data_d_bin_dmc <- data.frame(region_type = c(1:length(data$region_type)), cont_case = "Liver_cont_vs_Liver_case", d_bin_dmc = data$Liver_cont_vs_Liver_case.d_bin_dmc)
  data_d_bin_dmc <- rbind(data_d_bin_dmc, data.frame(region_type = c(1:length(data$region_type)), cont_case = "Liver_1_vs_Liver_2", d_bin_dmc = data$Liver_1_vs_Liver_2.d_bin_dmc))
  data_d_bin_dmc <- rbind(data_d_bin_dmc, data.frame(region_type = c(1:length(data$region_type)), cont_case = "Liver_1_vs_Liver_3", d_bin_dmc = data$Liver_1_vs_Liver_3.d_bin_dmc))
  data_d_bin_dmc <- rbind(data_d_bin_dmc, data.frame(region_type = c(1:length(data$region_type)), cont_case = "Liver_26_vs_Liver_2", d_bin_dmc = data$Liver_26_vs_Liver_2.d_bin_dmc))
  data_d_bin_dmc <- rbind(data_d_bin_dmc, data.frame(region_type = c(1:length(data$region_type)), cont_case = "Liver_26_vs_Liver_3", d_bin_dmc = data$Liver_26_vs_Liver_3.d_bin_dmc))

  data_plot <- ggplot(data_d_bin_dmc, aes(x=region_type, y=d_bin_dmc, col=cont_case)) + geom_line() +
    xlab('Up2k-GeneBody-Down2k') + ylab('Density of DMC in BIN') + 
    theme_bw() +
    scale_color_manual(name="Comparison pair", values=brewer.pal(6, 'Set1'))
  data_plot
  return (data_plot)
}

parse_DMC_ByBinSum_for_avg_diffdelta <- function(data) {
  data <- as.data.frame(data)
  data_avg_diffdelta <- data.frame(region_type = c(1:length(data$region_type)), cont_case = "Liver_cont_vs_Liver_case", avg_diffdelta = data$Liver_cont_vs_Liver_case.avg_diffdelta)
  data_avg_diffdelta <- rbind(data_avg_diffdelta, data.frame(region_type = c(1:length(data$region_type)), cont_case = "Liver_1_vs_Liver_2", avg_diffdelta = data$Liver_1_vs_Liver_2.avg_diffdelta))
  data_avg_diffdelta <- rbind(data_avg_diffdelta, data.frame(region_type = c(1:length(data$region_type)), cont_case = "Liver_1_vs_Liver_3", avg_diffdelta = data$Liver_1_vs_Liver_3.avg_diffdelta))
  data_avg_diffdelta <- rbind(data_avg_diffdelta, data.frame(region_type = c(1:length(data$region_type)), cont_case = "Liver_26_vs_Liver_2", avg_diffdelta = data$Liver_26_vs_Liver_2.avg_diffdelta))
  data_avg_diffdelta <- rbind(data_avg_diffdelta, data.frame(region_type = c(1:length(data$region_type)), cont_case = "Liver_26_vs_Liver_3", avg_diffdelta = data$Liver_26_vs_Liver_3.avg_diffdelta))
  
  data_plot <- ggplot(data_avg_diffdelta, aes(x=region_type, y=avg_diffdelta, col=cont_case)) + geom_line() +
    xlab('Up2k-GeneBody-Down2k') + ylab('Avg.Delta of DMC in BIN') + 
    theme_bw() +
    scale_color_manual(name="Comparison pair", values=brewer.pal(6, 'Set1'))
  data_plot
  return (data_plot)
}

mf_home <- "/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD190324/MineFigure"

context_type <- "CpG"
context_type <- "CHG"
context_type <- "CHH"

strand_type <- "OA" #all
strand_type <- "OB" #bottom
strand_type <- "OT" #top

data <- read.table(paste(mf_home, "AllGeneMethylRateDMC", paste("AllGeneMethylRateDMC", context_type, strand_type, "DMC.ByBinSum", sep = "."), sep = "/"), header = TRUE)

data_plot <- parse_DMC_ByBinSum_for_d_bin_dmc(data)
ggsave(paste(mf_home, "AllGeneMethylRateDMC", paste("AllGeneMethylRateDMC", context_type, strand_type, "DMC.ByBinSum", "density", "jpg", sep = "."), sep = "/"), dpi = 300, width = 10, height = 3)

data_plot <- parse_DMC_ByBinSum_for_avg_diffdelta(data)
ggsave(paste(mf_home, "AllGeneMethylRateDMC", paste("AllGeneMethylRateDMC", context_type, strand_type, "DMC.ByBinSum", "diffdelta", "jpg", sep = "."), sep = "/"), dpi = 300, width = 10, height = 3)

## Interested 
data <- read.table(paste(mf_home, "AllGeneMethylRateDMC_v2", paste("AllGeneMethylRateDMC", context_type, strand_type, "DMC.ByBinSum.interested.WBS", sep = "."), sep = "/"), header = TRUE)

data_plot <- parse_DMC_ByBinSum_for_d_bin_dmc(data)
ggsave(paste(mf_home, "AllGeneMethylRateDMC_v2", paste("AllGeneMethylRateDMC", context_type, strand_type, "DMC.ByBinSum.interested.WBS", "density", "jpg", sep = "."), sep = "/"), dpi = 300, width = 10, height = 3)

data_plot <- parse_DMC_ByBinSum_for_avg_diffdelta(data)
ggsave(paste(mf_home, "AllGeneMethylRateDMC_v2", paste("AllGeneMethylRateDMC", context_type, strand_type, "DMC.ByBinSum.interested.WBS", "diffdelta", "jpg", sep = "."), sep = "/"), dpi = 300, width = 10, height = 3)



#################################  all gene methylation rate by bin (with Interest) ################################# 

MethylRateByBinInterested_data_draw_lineplot <- function (data_mtx) {
  figure <- ggplot(data = data_mtx, aes(x = region, y = r_methyl, group = sample, col = color)) + geom_line(size=0.6) + 
    scale_y_continuous(limits=c(0, 100)) + 
    xlab('Up2k-GeneBody-Down2k') + ylab('Methylation Rate (%)') + 
    theme_bw() +
    scale_color_manual(name="Samples", values=brewer.pal(8, "Set1"))
  return(figure)
}

MethylRateByBinInterested_data_curation <- function(data){
  bin_order <- c(1:70)
  data_mtx <- data.frame(region = bin_order, r_methyl = data$S1.r_methyl, sample = "S1(WBS)", color = "S1(WBS)")
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$S2.r_methyl, sample = "S2(WBS)", color = "S2(WBS)"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$PIG_3.r_methyl, sample = "PIG_3(WBS)", color = "PIG_3(WBS)"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$PIG_4.r_methyl, sample = "PIG_4(WBS)", color = "PIG_4(WBS)"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$Liver_1.r_methyl, sample = "Liver_1(WBS)", color = "Liver_1(WBS)"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$Liver_2.r_methyl, sample = "Liver_2(WBS)", color = "Liver_2(WBS)"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$Liver_3.r_methyl, sample = "Liver_3(WBS)", color = "Liver_3(WBS)"))
  data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$Liver_26.r_methyl, sample = "Liver_26(WBS)", color = "Liver_26(WBS)"))
  return(data_mtx)
}

mf_home <- "/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD190324/MineFigure"

context_type <- "CpG"
context_type <- "CHG"
context_type <- "CHH"

strand_type <- "OA" #all
strand_type <- "OB" #bottom
strand_type <- "OT" #top

interested_name <- "WBS"

bin_order <- c(1:70)
data_all <- read.table(paste(mf_home, "AllGeneMethylRateByBin", paste("AllGeneMethylRateByBin", context_type, strand_type, "_methylByBinSum", sep = "."), sep = "/"), header = TRUE)
data_interested <- read.table(paste(mf_home, "AllGeneMethylRateByBin", paste("AllGeneMethylRateByBin", context_type, strand_type, "_methylByBin", "Interested", interested_name, "Sum", sep = "."), sep = "/"), header = TRUE)

###
sample_name <- "S1"
data_mtx <- data.frame(region = bin_order, r_methyl = data_all$"S1.r_methyl", sample = sample_name, color = sample_name)
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data_interested$"S1.r_methyl", sample = paste(sample_name, "(", interested_name, ")", sep=""), color  = paste(sample_name, "(", interested_name, ")", sep="")))
sample_name <- "S2"
data_mtx <- data.frame(region = bin_order, r_methyl = data_all$"S2.r_methyl", sample = sample_name, color = sample_name)
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data_interested$"S2.r_methyl", sample = paste(sample_name, "(", interested_name, ")", sep=""), color  = paste(sample_name, "(", interested_name, ")", sep="")))
sample_name <- "PIG_3"
data_mtx <- data.frame(region = bin_order, r_methyl = data_all$"PIG_3.r_methyl", sample = sample_name, color = sample_name)
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data_interested$"PIG_3.r_methyl", sample = paste(sample_name, "(", interested_name, ")", sep=""), color  = paste(sample_name, "(", interested_name, ")", sep="")))
sample_name <- "PIG_4"
data_mtx <- data.frame(region = bin_order, r_methyl = data_all$"PIG_4.r_methyl", sample = sample_name, color = sample_name)
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data_interested$"PIG_4.r_methyl", sample = paste(sample_name, "(", interested_name, ")", sep=""), color  = paste(sample_name, "(", interested_name, ")", sep="")))
sample_name <- "Liver_1"
data_mtx <- data.frame(region = bin_order, r_methyl = data_all$"Liver_1.r_methyl", sample = sample_name, color = sample_name)
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data_interested$"Liver_1.r_methyl", sample = paste(sample_name, "(", interested_name, ")", sep=""), color  = paste(sample_name, "(", interested_name, ")", sep="")))
sample_name <- "Liver_2"
data_mtx <- data.frame(region = bin_order, r_methyl = data_all$"Liver_2.r_methyl", sample = sample_name, color = sample_name)
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data_interested$"Liver_2.r_methyl", sample = paste(sample_name, "(", interested_name, ")", sep=""), color  = paste(sample_name, "(", interested_name, ")", sep="")))
sample_name <- "Liver_3"
data_mtx <- data.frame(region = bin_order, r_methyl = data_all$"Liver_3.r_methyl", sample = sample_name, color = sample_name)
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data_interested$"Liver_3.r_methyl", sample = paste(sample_name, "(", interested_name, ")", sep=""), color  = paste(sample_name, "(", interested_name, ")", sep="")))
sample_name <- "Liver_26"
data_mtx <- data.frame(region = bin_order, r_methyl = data_all$"Liver_26.r_methyl", sample = sample_name, color = sample_name)
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data_interested$"Liver_26.r_methyl", sample = paste(sample_name, "(", interested_name, ")", sep=""), color  = paste(sample_name, "(", interested_name, ")", sep="")))

(figure <- MethylRateByBinInterested_data_draw_lineplot(data_mtx))
ggsave(paste(mf_home, "AllGeneMethylRateByBin", paste("AllGeneMethylRateByBin", context_type, strand_type, "_methylByBin", "Interested", interested_name, "Sum", sample_name, 'jpg', sep = "."), sep = "/"), dpi = 300, width = 10, height = 3)

### 
data_mtx <- MethylRateByBinInterested_data_curation(data_interested)
(figure <- MethylRateByBinInterested_data_draw_lineplot(data_mtx))
ggsave(paste(mf_home, "AllGeneMethylRateByBin", paste("AllGeneMethylRateByBin", context_type, strand_type, "_methylByBin", "Interested", interested_name, "Sum", 'jpg', sep = "."), sep = "/"), dpi = 300, width = 10, height = 3)











