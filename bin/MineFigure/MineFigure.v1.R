
library(ggplot2)
library(RColorBrewer)


# all gene methylation rate
data_curation <- function (data) {
  data <- as.data.frame(data)
  data_r_methyl <- data.frame(r_methyl = data$S1.r_methyl, sample = 'S1', color = 'S1')
  data_r_methyl <- rbind(data_r_methyl, data.frame(r_methyl = data$S2.r_methyl, sample = 'S2', color = "S2"))
  data_r_methyl <- rbind(data_r_methyl, data.frame(r_methyl = data$PIG_3.r_methyl, sample = 'PIG_3', color = "PIG_3"))
  data_r_methyl <- rbind(data_r_methyl, data.frame(r_methyl = data$PIG_4.r_methyl, sample = 'PIG_4', color = "PIG_4"))
  return(data_r_methyl)
}

data_draw_boxplot <- function (data) {
  
  figure <- ggplot(data = data, aes(x = sample, y = r_methyl, group = sample, col = color)) + geom_boxplot() +
            scale_x_discrete(name = "Samples") +
            scale_y_continuous(name = "Methylation Rate in all genes (%)", breaks = seq(0, 100, 10), limits=c(0, 100)) +
            theme_bw() +
            scale_color_manual(name="Samples", values=c("black","#E41A1C","#377EB8","#4DAF4A") )
  
  return(figure)
}

data <- read.table("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRate/AllGeneMethylRate.CpG.OA._methyl.extend", header = TRUE)
data <- data_curation(data)
figure <- data_draw_boxplot(data)
ggsave("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRate/AllGeneMethylRate.CpG.jpg", dpi = 300, width = 7, height = 5)

data <- read.table("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRate/AllGeneMethylRate.CHG.OA._methyl.extend", header = TRUE)
data <- data_curation(data)
figure <- data_draw_boxplot(data)
ggsave("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRate/AllGeneMethylRate.CHG.jpg", dpi = 300, width = 7, height = 5)

data <- read.table("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRate/AllGeneMethylRate.CHH.OA._methyl.extend", header = TRUE)
data <- data_curation(data)
figure <- data_draw_boxplot(data)
ggsave("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRate/AllGeneMethylRate.CHH.jpg", dpi = 300, width = 7, height = 5)


# all gene methylation rate by bin

data <- read.table("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRateByBin/AllGeneMethylRateByBin.CpG.OA._methylByBinSum", header = TRUE)
bin_order <- c(1:70)

data_mtx <- data.frame(region = bin_order, r_methyl = data$S1.r_methyl, sample = "S1", color = "S1")
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$S2.r_methyl, sample = "S2", color = "S2"))
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$PIG_3.r_methyl, sample = "PIG_3", color = "PIG_3"))
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$PIG_4.r_methyl, sample = "PIG_4", color = "PIG_4"))

ggplot(data = data_mtx, aes(x = region, y = r_methyl, group = sample, col = color)) + geom_line(size=0.6) + 
  scale_y_continuous(limits=c(0, 100)) + 
  xlab('Up2k-GeneBody-Down2k') + ylab('Methylation Rate (%)') + 
  theme_bw() +
  scale_color_manual(name="Samples", values=c("black","#E41A1C","#377EB8","#4DAF4A"))

ggsave("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRateByBin/MethylationByBin.CpG.jpg", dpi = 300, width = 10, height = 3)


# all gene DMC
all_gene_methylation_rate_DMC <- function(data) {
  data <- as.data.frame(data)
  data_d_bin_dmc <- data.frame(region_type = c(1:length(data$region_type)), cont_case = "S1_vs_S2", d_bin_dmc = data$S1_vs_S2.d_bin_dmc)
  data_d_bin_dmc <- rbind(data_d_bin_dmc, data.frame(region_type = c(1:length(data$region_type)), cont_case = "S1_vs_PIG_3", d_bin_dmc = data$S1_vs_PIG_3.d_bin_dmc))
  data_d_bin_dmc <- rbind(data_d_bin_dmc, data.frame(region_type = c(1:length(data$region_type)), cont_case = "S1_vs_PIG_4", d_bin_dmc = data$S1_vs_PIG_4.d_bin_dmc))
  data_d_bin_dmc <- rbind(data_d_bin_dmc, data.frame(region_type = c(1:length(data$region_type)), cont_case = "S2_vs_PIG_3", d_bin_dmc = data$S2_vs_PIG_3.d_bin_dmc))
  data_d_bin_dmc <- rbind(data_d_bin_dmc, data.frame(region_type = c(1:length(data$region_type)), cont_case = "S2_vs_PIG_4", d_bin_dmc = data$S2_vs_PIG_4.d_bin_dmc))
  data_d_bin_dmc <- rbind(data_d_bin_dmc, data.frame(region_type = c(1:length(data$region_type)), cont_case = "PIG_3_vs_PIG_4", d_bin_dmc = data$PIG_3_vs_PIG_4.d_bin_dmc))
  
  data_plot <- ggplot(data_d_bin_dmc, aes(x=region_type, y=d_bin_dmc, col=cont_case)) + geom_line() +
    xlab('Up2k-GeneBody-Down2k') + ylab('Density of DMC in BIN') + 
    theme_bw() +
    scale_color_manual(name="Comparison pair", values=brewer.pal(6, 'Set1'))
  
  return (data_plot)
}

data <- read.table("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRateDMC/AllGeneMethylRateDMC.CHH.OA.DMC.ByBinSum", header = TRUE)
all_gene_methylation_rate_DMC(data)
ggsave("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRateDMC/DMCByBin.CHH.jpg", dpi = 300, width = 10, height = 3)


# all gene methylation rate by bin (with Interest)

data <- read.table("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRateByBin/AllGeneMethylRateByBin.CpG.OA._methylByBinSum", header = TRUE)
bin_order <- c(1:70)
data_mtx <- data.frame(region = bin_order, r_methyl = data$S1.r_methyl, sample = "S1", color = "S1")
data_mtx <- data.frame(region = bin_order, r_methyl = data$S2.r_methyl, sample = "S2", color = "S2")
data_mtx <- data.frame(region = bin_order, r_methyl = data$PIG_3.r_methyl, sample = "PIG_3", color = "PIG_3")
data_mtx <- data.frame(region = bin_order, r_methyl = data$PIG_4.r_methyl, sample = "PIG_4", color = "PIG_4")

data <- read.table("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRateByBin/AllGeneMethylRateByBin.CpG.OA._methylByBin.Interested.WBS.Sum", header = TRUE)
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$S1.r_methyl, sample = "S1(WBS)", color = "S1(WBS)"))
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$S2.r_methyl, sample = "S2(WBS)", color = "S2(WBS)"))
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$PIG_3.r_methyl, sample = "PIG_3(WBS)", color = "PIG_3(WBS)"))
data_mtx <- rbind(data_mtx, data.frame(region = bin_order, r_methyl = data$PIG_4.r_methyl, sample = "PIG_4(WBS)", color = "PIG_4(WBS)"))

ggplot(data = data_mtx, aes(x = region, y = r_methyl, group = sample, col = color)) + geom_line(size=0.6) + 
  scale_y_continuous(limits=c(0, 100)) + 
  xlab('Up2k-GeneBody-Down2k') + ylab('Methylation Rate (%)') + 
  theme_bw() +
  scale_color_manual(name="Samples", values=c("black","#E41A1C","#377EB8","#4DAF4A"))

ggsave("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository/TBD180834-ALLBT-Pig-Bisulfite-Report-20190222/MineFigure/AllGeneMethylRateByBin/MethylationByBin.CpG.PIG_4.WBS.jpg", dpi = 300, width = 10, height = 3)


