#Load in necessary libraries:
library(tidyverse)
library(reshape2)
library(cowplot)

#Functions for plotting:
reorder_demes <- function(Qmatrix, focal_indivs, K, Kmax, all_demes) {
   focal_demes <- apply(subset(Qmatrix, Individual %in% focal_indivs, select=all_demes), 1, which.max);
   names(focal_demes) <- subset(Qmatrix, Individual %in% focal_indivs)$Individual;
   focal_indivs_indistinct <- FALSE;
   if (length(unique(focal_demes)) < length(focal_demes)) {
      focal_indivs_indistinct <- TRUE;
   }
   pops <- c();
   pops[1] <- paste0("Deme", focal_demes[focal_indivs[1]]);
   for (i in seq(2, Kmax-1)) {
      if (K < i || focal_indivs_indistinct) {
         pops[i] <- sample(setdiff(all_demes, pops), 1);
      } else {
         pops[i] <- paste0("Deme", focal_demes[focal_indivs[i]]);
      }
   }
   pops[Kmax] <- setdiff(all_demes, pops);
   pops
}

prepare_maverick_df <- function(MavericK, arm, focal_indivs, Kmax) {
   all_demes <- paste0("Deme", 1:Kmax);
   MavericK_df <- data.frame(Window=c(), K=c(), Individual=c());
   for (h in seq(1, Kmax)) {
      MavericK_df[,paste0("Deme", h)] <- c();
   }
   #Reorder demes:
   for (i in seq(1, length(MavericK[[arm]]))) {
      K <- MavericK[[arm]][[i]][1,"K"];
      ordered_pops <- reorder_demes(MavericK[[arm]][[i]], focal_indivs, K, Kmax, all_demes);
      Qdf <- MavericK[[arm]][[i]][,c("Window", "K", "Individual")];
      for (j in seq(1, Kmax)) {
         Qdf[,paste0("Deme", j)] <- MavericK[[arm]][[i]][,ordered_pops[j]];
      }
      MavericK_df <- rbind(MavericK_df, Qdf);
   }
   #Windows are in units of 100 kb, so rescale to Mb:
   MavericK_df$Window.Start <- (MavericK_df$Window-1)/10;
   #Melt the DF for plotting:
   MavericK_df_melted <- melt(MavericK_df, id.vars=c("Window", "Window.Start", "K", "Individual"), variable.name="Deme", value.name="Ancestry.Fraction")
   MavericK_df_melted
}

plot_maverick <- function(melted_df, arm, num_indivs, breakpoints) {
   if (nrow(breakpoints) > 0) {
      ggplot(melted_df, aes(x=Window.Start, y=Ancestry.Fraction, fill=Deme, group=Individual)) +
         geom_col() +
         geom_hline(yintercept=seq(0,num_indivs)) +
         geom_segment(data=breakpoints, aes(x=x, xend=xend, y=y, yend=yend, colour=Inversion), size=2, inherit.aes=FALSE) +
         xlab("Window Start (Mb)") +
         ylab("Ancestry Fraction") +
         ggtitle(paste0("MavericK for 100kb windows of chr", arm)) +
         scale_fill_brewer(palette="Paired") +
         scale_colour_brewer(palette="Set2") +
         guides(size=FALSE, fill=FALSE, colour=guide_legend(override.aes=list(size=4))) +
         theme_bw() +
         theme(plot.title=element_text(size=10),
               axis.text=element_text(size=10),
               axis.title=element_text(size=10, face="bold"),
               axis.title.x=element_text(vjust=-4),
               axis.title.y=element_text(vjust=4),
               legend.title=element_blank(),
               legend.text=element_text(size=10),
               legend.text.align=0,
               legend.position="top",
               legend.spacing.x=unit(0.3, "cm"),
               legend.box.margin=margin(t=-0.0, r=0,
                                        b=-0.0, l=0, unit='cm'),
               legend.margin=margin(t=-0.3, r=0,
                                    b=-0.3, l=0, unit="cm"),
               strip.text=element_text(size=10),
               plot.margin=margin(t=0.2, r=0.5, b=1.0, l=1.0, unit="cm")) +
         scale_y_continuous(breaks=seq(0.5,num_indivs-0.5), labels=rev(levels(melted_df$Individual)), expand=c(0.01,0.001))
   } else {
      ggplot(melted_df, aes(x=Window.Start, y=Ancestry.Fraction, fill=Deme, group=Individual)) +
         geom_col() +
         geom_hline(yintercept=seq(0,num_indivs)) +
         xlab("Window Start (Mb)") +
         ylab("Ancestry Fraction") +
         ggtitle(paste0("MavericK for 100kb windows of chr", arm)) +
         scale_fill_brewer(palette="Paired") +
         guides(fill=FALSE) +
         theme_bw() +
         theme(plot.title=element_text(size=10),
               axis.text=element_text(size=10),
               axis.title=element_text(size=10, face="bold"),
               axis.title.x=element_text(vjust=-4),
               axis.title.y=element_text(vjust=4),
               legend.title=element_text(size=10, face="bold"),
               legend.text=element_text(size=10),
               legend.text.align=0,
               strip.text=element_text(size=10),
               plot.margin=margin(t=0.2, r=0.5, b=1.0, l=1.0, unit="cm")) +
         scale_y_continuous(breaks=seq(0.5,num_indivs-0.5), labels=rev(levels(melted_df$Individual)), expand=c(0.001,0.001))
   }
}

#Load in the MavericK results:
load('All_Dyak_MavericK_results.Rdata')

#Prep the data.frames for plotting:
inbred_lines <- c("NY73PB", "Tai18E2", "R9b", "2.8.1", "ST", "IC261", "Jess", "CY01A", "CY02B5", "CY04B", "CY08A", "CY13A", "CY17C", "CY20A", "CY21B3", "CY22B", "CY28", "NY42", "NY48", "NY56", "NY62", "NY65", "NY66-2", "NY73", "NY81", "NY85", "NY141")
SD_lines <- c("SD1", "SD2", "SD3", "SD4", "SD5", "SD6", "SD7", "SD8", "SD9", "SD10", "SD11", "SD12", "SD13", "SD14", "SD15", "SD16", "SD17", "SD18", "SD19", "SD20", "SD21", "SD22", "sd25", "sd26", "sd27", "sd28", "sd29", "sd30", "sd31", "sd32")
Breakpoints <- list("X"=data.frame(), "2L"=data.frame(x=c(8.564383), xend=c(13.799295), y=c(27.5), yend=c(27.5), Inversion=c("In(2L)")), "2R"=data.frame(x=c(11.406171, 11.802385, 15.500001), xend=c(20.875023, 18.508443, 19.200001), y=c(27.5, 28.5, 29.5), yend=c(27.5, 28.5, 29.5), Inversion=c("2Rj", "2Rn", "2Rk")), "3L"=data.frame(), "3R"=data.frame())
focal_indivs_2L <- c("Tai18E2", "NY73PB")
focal_indivs_2R <- c("NY73PB", "Tai18E2")
Kmax <- 3

MavericK_dfs <- list()
MavericK_dfs[["X"]] <- prepare_maverick_df(MavericK, "X", focal_indivs_2L, Kmax)
MavericK_dfs[["2L"]] <- prepare_maverick_df(MavericK, "2L", focal_indivs_2L, Kmax)
MavericK_dfs[["2R"]] <- prepare_maverick_df(MavericK, "2R", focal_indivs_2R, Kmax)
MavericK_dfs[["3L"]] <- prepare_maverick_df(MavericK, "3L", focal_indivs_2L, Kmax)
MavericK_dfs[["3R"]] <- prepare_maverick_df(MavericK, "3R", focal_indivs_2L, Kmax)

#Do the actual plotting:
MavericK_plots <- list()
MavericK_plots[["X"]] <- MavericK_dfs[["X"]] %>% filter(Individual %in% inbred_lines) %>% mutate(Individual=factor(Individual, levels=inbred_lines)) %>% plot_maverick(., "X", length(inbred_lines), Breakpoints[["X"]])
MavericK_plots[["2L"]] <- MavericK_dfs[["2L"]] %>% filter(Individual %in% inbred_lines) %>% mutate(Individual=factor(Individual, levels=inbred_lines)) %>% plot_maverick(., "2L", length(inbred_lines), Breakpoints[["2L"]])
MavericK_plots[["2R"]] <- MavericK_dfs[["2R"]] %>% filter(Individual %in% inbred_lines) %>% mutate(Individual=factor(Individual, levels=inbred_lines)) %>% plot_maverick(., "2R", length(inbred_lines), Breakpoints[["2R"]])
MavericK_plots[["3L"]] <- MavericK_dfs[["3L"]] %>% filter(Individual %in% inbred_lines) %>% mutate(Individual=factor(Individual, levels=inbred_lines)) %>% plot_maverick(., "3L", length(inbred_lines), Breakpoints[["3L"]])
MavericK_plots[["3R"]] <- MavericK_dfs[["3R"]] %>% filter(Individual %in% inbred_lines) %>% mutate(Individual=factor(Individual, levels=inbred_lines)) %>% plot_maverick(., "3R", length(inbred_lines), Breakpoints[["3R"]])
#Supplemental figure of the X:
ggsave('FigSX_MavericK_Dyak_inbred_chrX.pdf', plot=MavericK_plots[["X"]], width=12.0, height=16.0, units="cm", dpi=500)

#Figure 4 of chromosome 2:
Fig4_grid <- plot_grid(MavericK_plots[["2L"]], MavericK_plots[["2R"]], align="h", axis="bt", ncol=2)
save_plot('Fig4_MavericK_Dyak_inbred_chr2.pdf', Fig4_grid, ncol=2, base_height=16.0/2.54, base_width=12.0/2.54, dpi=500)

#Supplemental figure of chromosome 3:
FigSY_grid <- plot_grid(MavericK_plots[["3L"]], MavericK_plots[["3R"]], align="h", axis="bt", ncol=2)
save_plot('FigSY_MavericK_Dyak_inbred_chr3.pdf', FigSY_grid, ncol=2, base_height=16.0/2.54, base_width=12.0/2.54, dpi=500)
