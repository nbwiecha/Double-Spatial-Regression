# Simulation result figures 
# for Two-Stage Estimators for Spatial Confounding with Point-Referenced Data
# Code by Nate Wiecha, North Carolina State University

rm(list=ls())
setwd("~/GitHub/Double-Spatial-Regression")
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(cowplot)
source("code//HPC//simulation_functions_hpc.R")

load("outputs//HPC_outputs_revision//sims_deterministic_n1000.Rdata")

load("outputs//HPC_outputs_revision//sims_highdimvar_n1000.Rdata")

load("outputs//HPC_outputs_revision//sims_exp_nonspatial_grid_n1000.Rdata")

load("outputs//HPC_outputs_revision//sims_U_gneiting_n1000.Rdata")

load("outputs//HPC_outputs_revision//sims_U_matern_n1000.Rdata")

load("outputs//HPC_outputs_revision//sims_trans_n1000.Rdata")

# Coverage table for paper
U_matern_sims$out1$output["95% CI Coverage",c("DML_alt", "DML_GP", "gSEM_noint", "LMM", "SpatialPlus")]
U_matern_sims$out2$output["95% CI Coverage",c("DML_alt", "DML_GP", "gSEM_noint", "LMM", "SpatialPlus")]
U_gneiting_sims$out1$output["95% CI Coverage",c("DML_alt", "DML_GP", "gSEM_noint", "LMM", "SpatialPlus")]
U_gneiting_sims$out2$output["95% CI Coverage",c("DML_alt", "DML_GP", "gSEM_noint", "LMM", "SpatialPlus")]

U_matern_sims$out1$output["95% CI Length",c("DML_alt", "DML_GP", "gSEM_noint", "LMM", "SpatialPlus")]
U_matern_sims$out2$output["95% CI Length",c("DML_alt", "DML_GP", "gSEM_noint", "LMM", "SpatialPlus")]
U_gneiting_sims$out1$output["95% CI Length",c("DML_alt", "DML_GP", "gSEM_noint", "LMM", "SpatialPlus")]
U_gneiting_sims$out2$output["95% CI Length",c("DML_alt", "DML_GP", "gSEM_noint", "LMM", "SpatialPlus")]

################################################################################
#                                  Plots                                       #
################################################################################
methods <- c("LMM", "SpatialPlus", "gSEM_noint", "DML_GP",  "DML_alt")
method_names <- c("LMM", "Spatial+", "gSEM", "DSR (theory)",  "DSR")
make_boxplots <- function(output, methods, method_names){
  num_plots <- length(output)/2
  plots_point <- list(length=num_plots)
  # taus <- c(1.5, 2.5, 50)
  subtitles <- c("Rough A", "Smooth A")
  for(i in 1:num_plots){
    # Make point est boxplots
    dat <- output[[i]]$results %>% 
      as.data.frame %>% 
      select(ends_with(".thetahat"))
    colnames(dat) <- str_remove(colnames(dat), ".thetahat")
    dat <- select(dat, all_of(methods)) 
    colnames(dat) <- method_names
    dat <- dat %>%
      pivot_longer(everything()) %>%
      mutate(Method=as.factor(name),
             Bias=value - 0.5)
    subtitle <- subtitles[i]
    plots_point[[i]] <- ggplot(dat, aes(x=Method, y=Bias)) + 
      geom_boxplot() +
      geom_hline(yintercept=0) +
      ggtitle(subtitle) +
      theme_classic(base_size = 25) +
      ylab(bquote(hat(beta) - beta[0])) +
      theme(text = element_text(family = "serif"),
            legend.key.size=unit(1,'cm'),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      ylim(-1.1, 1.1)
  }
  
  return(plots_point)
}



plots_U_matern <- make_boxplots(U_matern_sims, methods, method_names)

plots_U_matern_grid <- ggarrange(plotlist=plots_U_matern)

U_matern.title <- ggdraw() + 
  draw_label(
    bquote(hat(beta)-beta[0]~with~rough~confounder),
    fontfamily = 'serif',
    x = 0,
    hjust = 0,
    size=30
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

# arrange title and grid
U_matern.plot <- plot_grid(
  U_matern.title, plots_U_matern_grid,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

plots_U_gneiting <- make_boxplots(U_gneiting_sims, methods, method_names)

plots_U_gneiting_grid <- ggarrange(plotlist=plots_U_gneiting)

U_gneiting.title <- ggdraw() + 
  draw_label(
    bquote(hat(beta)-beta[0]~with~smooth~confounder),
    fontfamily = 'serif',
    x = 0,
    hjust = 0,
    size=30
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

# arrange title and grid
U_gneiting.plot <- plot_grid(
  U_gneiting.title, plots_U_gneiting_grid,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

png("outputs//U_matern_plot.png", width = 1000, height = 800)
U_matern.plot
dev.off()

png("outputs//U_gneiting_plot.png", width = 1000, height = 800)
U_gneiting.plot
dev.off()
