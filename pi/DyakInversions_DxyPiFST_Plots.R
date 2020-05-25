#Libraries necessary for processing and plotting:
library(tidyverse)
library(cowplot)

#General variables for IO and filtering
diversity_col_classes <- c("character", "integer", "numeric", "numeric")
min_sites_thresh <- 0.004
#Plotting details:
fontsize <- 16
#Genomewide pi (4-fold of course):
pi_genomewide <- list("Dyak"=0.0218465, "Dsan"=0.015436)
#Y axis boundaries:
ymax <- list("pi_between"=0.09,
             "pi_within"=0.05,
             "F_ST"=0.9)
#Scaffold boundaries for use with xlim(), in Mb scale:
#Better way might have been to read in from FAI
bounds <- list("2L"=c(1, 24263191)/1e6, "2R"=c(1, 23843395)/1e6)
#Inversion spans (NY73PB coordinates):
inversions <- list("2L"=data.frame(x=c(8564383)/1e6,
                                   xend=c(13799295)/1e6,
                                   y=c(0/4),
                                   yend=c(0.8/4),
                                   Inversion=c("In(2L)")),
                   "2R"=data.frame(x=c(11406171, 11802385, 15500001)/1e6,
                                   xend=c(20875023, 18508443, 19200001)/1e6,
                                   y=c(0/4, 1.0/4, 2.0/4),
                                   yend=c(0.8/4, 1.8/4, 2.8/4),
                                   Inversion=factor(c("In(2R)j", "In(2R)n", "In(2R)k"),
                                                    levels=c("In(2R)k", "In(2R)n", "In(2R)j"))))

#Read in the 2L data:
#pi in the standard arrangement:
diversity_2L <- read.table(gzfile('Dyak_Std2L_1_In2L_2_pi_S_f0.5_4fold_w100kb.tsv.gz'),
                           header=FALSE,
                           colClasses=diversity_col_classes,
                           col.names=c("Scaffold", "WindowStart", "pi_S", "UsableFraction_pi_S"))
#pi in the inverted arrangement:
diversity_2L <- read.table(gzfile('Dyak_Std2L_1_In2L_2_pi_I_f0.5_4fold_w100kb.tsv.gz'),
                           header=FALSE,
                           colClasses=diversity_col_classes,
                           col.names=c("Scaffold", "WindowStart", "pi_I", "UsableFraction_pi_I")) %>%
   full_join(diversity_2L, by=c("Scaffold", "WindowStart"))
#d_xy between arrangements:
diversity_2L <- read.table(gzfile('Dyak_Std2L_1_In2L_2_D_SI_f0.5_4fold_w100kb.tsv.gz'),
                           header=FALSE,
                           colClasses=diversity_col_classes,
                           col.names=c("Scaffold", "WindowStart", "D_SI", "UsableFraction_D_SI")) %>%
   full_join(diversity_2L, by=c("Scaffold", "WindowStart"))
#total pi amongst arrangements (here, equivalent to total pi in Dyak):
#diversity_2L <- read.table(gzfile('Dyak_Std2L_1_In2L_2_pi_T_SI_f0.5_4fold_w100kb.tsv.gz'),
#                           header=FALSE,
#                           colClasses=diversity_col_classes,
#                           col.names=c("Scaffold", "WindowStart", "pi_T_SI", "UsableFraction_pi_T_SI")) %>%
#   full_join(diversity_2L, by=c("Scaffold", "WindowStart"))
#d_xy between species:
diversity_2L <- read.table(gzfile('Dsan_projections/Dyak_1_Dsan_2_D_YS_f0.5_4fold_w100kb.tsv.gz'),
                           header=FALSE,
                           colClasses=diversity_col_classes,
                           col.names=c("Scaffold", "WindowStart", "D_YS", "UsableFraction_D_YS")) %>%
   full_join(diversity_2L, by=c("Scaffold", "WindowStart"))

#Make Position into Mb:
diversity_2L <- diversity_2L %>%
   mutate(Position=WindowStart/1e6)

#Hudson, Slatkin, Maddison F_ST estimator: 1-H_w/H_b
#Their H_w is the average of within-pop pi
#Their H_b is D_xy between pops
diversity_2L <- diversity_2L %>%
   mutate(F_ST_SI=1-(pi_S+pi_I)/(D_SI*2))

#Filter for usable fraction of sites in a window, and pivot df:
diversity_2L <- diversity_2L %>%
   filter(UsableFraction_pi_S >= min_sites_thresh,
          UsableFraction_pi_I >= min_sites_thresh,
          UsableFraction_D_SI >= min_sites_thresh,
          UsableFraction_D_YS >= min_sites_thresh) %>%
   gather("Statistic", "Value", -Scaffold, -WindowStart, -Position)

#Make the four plots for the column:
#1) Rectangle of the 2L inversion
#2) pi between arrangements
#3) pi within arrangements
#4) F_ST between arrangements
#legends go on top of each plot, make sure to rename statistics to pretty things
inversion_2L <- inversions[["2L"]] %>%
   ggplot(aes(x=x-1.5, xmin=x, xmax=xend, y=y+0.4/4, ymin=y, ymax=yend, fill=Inversion, label=Inversion)) +
      geom_rect() +
      geom_text() +
      scale_x_continuous(limits=bounds[["2L"]],
                         expand=expand_scale(mult=c(0.05, 0.06))) +
      ylim(0, 3/4) +
      scale_fill_brewer(palette="Paired") +
      labs(title="2L") +
      theme_nothing() +
      theme(text=element_text(size=fontsize-2),
            plot.title=element_text(size=fontsize+2,
                                    margin=margin(b=12)))
print(inversion_2L)

pi_between_2L <- diversity_2L %>%
   filter(Scaffold == "2L",
          Statistic %in% c("D_SI", "D_YS")) %>%
   mutate(Statistic=factor(Statistic, levels=c("D_YS", "D_SI"))) %>%
   ggplot(aes(x=Position, y=Value, colour=Statistic)) +
      geom_line() +
      geom_hline(yintercept=pi_genomewide[["Dyak"]], linetype=2, alpha=0.5) +
      scale_colour_brewer(palette="Set2",
                          labels=c(expression(d[~S-Y]),
                                   expression(pi[~INV-STD]))) +
      labs(x="Position (Mb)", y=expression(pi)) +
      guides(fill=FALSE, colour=guide_legend(override.aes=list(size=3))) +
      theme_bw() +
      ylim(0, ymax[["pi_between"]]) +
      scale_x_continuous(limits=bounds[["2L"]],
                         expand=expand_scale(mult=c(0.05, 0.06))) +
      theme(plot.title=element_text(size=fontsize, face="bold"),
            axis.text=element_text(size=fontsize),
            axis.title=element_text(size=fontsize, face="bold"),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.y=element_text(vjust=4),
            #legend.title=element_text(size=fontsize, face="bold"),
            legend.title=element_blank(),
            legend.text=element_text(size=fontsize-2),
            legend.text.align=0,
            legend.position="top",
            legend.spacing.x=unit(0.3, "cm"),
            legend.box.margin=margin(t=-0.1, r=0,
                                     b=-0.2, l=0, unit='cm'),
            strip.text=element_text(size=fontsize),
            plot.margin=margin(t=0.2, r=0.0, b=0.0, l=1.0, unit="cm"))
print(pi_between_2L)

pi_within_2L <- diversity_2L %>%
   filter(Scaffold == "2L",
          Statistic %in% c("pi_S", "pi_I")) %>%
   mutate(Statistic=factor(Statistic, levels=c("pi_I", "pi_S"))) %>%
   ggplot(aes(x=Position, y=Value, colour=Statistic)) +
      geom_line() +
      geom_hline(yintercept=pi_genomewide[["Dyak"]], linetype=2, alpha=0.5) +
      scale_colour_brewer(palette="Set2",
                          labels=c(expression(pi[~INV]),
                                   expression(pi[~STD]))) +
      labs(x="Position (Mb)", y=expression(pi)) +
      guides(fill=FALSE, colour=guide_legend(override.aes=list(size=3))) +
      theme_bw() +
      ylim(0, ymax[["pi_within"]]) +
      scale_x_continuous(limits=bounds[["2L"]],
                         expand=expand_scale(mult=c(0.05, 0.06))) +
      theme(plot.title=element_text(size=fontsize, face="bold"),
            axis.text=element_text(size=fontsize),
            axis.title=element_text(size=fontsize, face="bold"),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.y=element_text(vjust=4),
            #legend.title=element_text(size=fontsize, face="bold"),
            legend.title=element_blank(),
            legend.text=element_text(size=fontsize-2),
            legend.text.align=0,
            legend.position="top",
            legend.spacing.x=unit(0.3, "cm"),
            legend.box.margin=margin(t=-0.1, r=0,
                                     b=-0.3, l=0, unit='cm'),
            strip.text=element_text(size=fontsize),
            plot.margin=margin(t=0.2, r=0.0, b=0.0, l=1.0, unit="cm"))
print(pi_within_2L)

F_ST_2L <- diversity_2L %>%
   filter(Scaffold == "2L",
          Statistic %in% c("F_ST_SI")) %>%
   ggplot(aes(x=Position, y=Value, colour=Statistic)) +
      geom_line() +
      scale_colour_brewer(palette="Set2",
                          labels=c(expression(F[ST]^{~STD-INV}))) +
      labs(x="Position (Mb)", y=expression(F[ST])) +
      guides(fill=FALSE, colour=guide_legend(override.aes=list(size=3))) +
      theme_bw() +
      ylim(0, ymax[["F_ST"]]) +
      scale_x_continuous(limits=bounds[["2L"]],
                         expand=expand_scale(mult=c(0.05, 0.06))) +
      theme(plot.title=element_text(size=fontsize, face="bold"),
            axis.text=element_text(size=fontsize),
            axis.title=element_text(size=fontsize, face="bold"),
            axis.title.x=element_text(vjust=-4),
            axis.title.y=element_text(vjust=4),
            #legend.title=element_text(size=fontsize, face="bold"),
            legend.title=element_blank(),
            legend.text=element_text(size=fontsize-2),
            legend.text.align=0,
            legend.position="top",
            legend.spacing.x=unit(0.3, "cm"),
            legend.box.margin=margin(t=-0.1, r=0,
                                     b=-0.3, l=0, unit='cm'),
            strip.text=element_text(size=fontsize),
            plot.margin=margin(t=0.2, r=0.0, b=1.0, l=1.0, unit="cm"))
print(F_ST_2L)

#2R
#Read in the 2R data:
#pi in the inverted arrangements:
diversity_2R <- read.table(gzfile('Dyak_2Rj_1_2Rn_2_2Rk_3_pi_2Rj_f0.5_4fold_w100kb.tsv.gz'),
                           header=FALSE,
                           colClasses=diversity_col_classes,
                           col.names=c("Scaffold", "WindowStart", "pi_2Rj", "UsableFraction_pi_2Rj"))
diversity_2R <- read.table(gzfile('Dyak_2Rj_1_2Rn_2_2Rk_3_pi_2Rn_f0.5_4fold_w100kb.tsv.gz'),
                           header=FALSE,
                           colClasses=diversity_col_classes,
                           col.names=c("Scaffold", "WindowStart", "pi_2Rn", "UsableFraction_pi_2Rn")) %>%
   full_join(diversity_2R, by=c("Scaffold", "WindowStart"))
diversity_2R <- read.table(gzfile('Dyak_2Rj_1_2Rn_2_2Rk_3_pi_2Rk_f0.5_4fold_w100kb.tsv.gz'),
                           header=FALSE,
                           colClasses=diversity_col_classes,
                           col.names=c("Scaffold", "WindowStart", "pi_2Rk", "UsableFraction_pi_2Rk")) %>%
   full_join(diversity_2R, by=c("Scaffold", "WindowStart"))
#d_xy between arrangements:
diversity_2R <- read.table(gzfile('Dyak_2Rj_1_2Rn_2_2Rk_3_D_jn_f0.5_4fold_w100kb.tsv.gz'),
                           header=FALSE,
                           colClasses=diversity_col_classes,
                           col.names=c("Scaffold", "WindowStart", "D_jn", "UsableFraction_D_jn")) %>%
   full_join(diversity_2R, by=c("Scaffold", "WindowStart"))
diversity_2R <- read.table(gzfile('Dyak_2Rj_1_2Rn_2_2Rk_3_D_jk_f0.5_4fold_w100kb.tsv.gz'),
                           header=FALSE,
                           colClasses=diversity_col_classes,
                           col.names=c("Scaffold", "WindowStart", "D_jk", "UsableFraction_D_jk")) %>%
   full_join(diversity_2R, by=c("Scaffold", "WindowStart"))
#total pi amongst arrangements (here, equivalent to total pi in Dyak):
#diversity_2R <- read.table(gzfile('Dyak_2Rj_1_2Rn_2_2Rk_3_pi_T_jn_f0.5_4fold_w100kb.tsv.gz'),
#                           header=FALSE,
#                           colClasses=diversity_col_classes,
#                           col.names=c("Scaffold", "WindowStart", "pi_T_jn", "UsableFraction_pi_T_jn")) %>%
#   full_join(diversity_2R, by=c("Scaffold", "WindowStart"))
#diversity_2R <- read.table(gzfile('Dyak_2Rj_1_2Rn_2_2Rk_3_pi_T_jk_f0.5_4fold_w100kb.tsv.gz'),
#                           header=FALSE,
#                           colClasses=diversity_col_classes,
#                           col.names=c("Scaffold", "WindowStart", "pi_T_jk", "UsableFraction_pi_T_jk")) %>%
#   full_join(diversity_2R, by=c("Scaffold", "WindowStart"))
#d_xy between species:
diversity_2R <- read.table(gzfile('Dsan_projections/Dyak_1_Dsan_2_D_YS_f0.5_4fold_w100kb.tsv.gz'),
                           header=FALSE,
                           colClasses=diversity_col_classes,
                           col.names=c("Scaffold", "WindowStart", "D_YS", "UsableFraction_D_YS")) %>%
   full_join(diversity_2R, by=c("Scaffold", "WindowStart"))


#Make Position into Mb:
diversity_2R <- diversity_2R %>%
   mutate(Position=WindowStart/1e6)

#Hudson, Slatkin, Maddison F_ST estimator: 1-H_w/H_b
#Their H_w is the average of within-pop pi
#Their H_b is D_xy between pops
diversity_2R <- diversity_2R %>%
   mutate(F_ST_jn=1-(pi_2Rj+pi_2Rn)/(D_jn*2),
          F_ST_jk=1-(pi_2Rj+pi_2Rk)/(D_jk*2))

#Filter for usable fraction of sites in a window, and pivot df:
diversity_2R <- diversity_2R %>%
   filter(UsableFraction_pi_2Rj >= min_sites_thresh,
          UsableFraction_pi_2Rk >= min_sites_thresh,
          UsableFraction_pi_2Rn >= min_sites_thresh,
          UsableFraction_D_jn >= min_sites_thresh,
          UsableFraction_D_jk >= min_sites_thresh,
          UsableFraction_D_YS >= min_sites_thresh) %>%
   gather("Statistic", "Value", -Scaffold, -WindowStart, -Position)


#Make the four plots for the column:
#1) Rectangle of the 2L inversion
#2) pi between arrangements
#3) pi within arrangements
#4) F_ST between arrangements
#legends go on top of each plot, make sure to rename statistics to pretty things
inversion_2R <- inversions[["2R"]] %>%
   ggplot(aes(x=x-1.5, xmin=x, xmax=xend, y=y+0.4/4, ymin=y, ymax=yend, fill=Inversion, label=Inversion)) +
      geom_rect() +
      geom_text() +
      scale_x_continuous(limits=bounds[["2R"]],
                         expand=expand_scale(mult=c(0.05, 0.06))) +
      ylim(0, 3/4) +
      scale_fill_brewer(palette="Paired") +
      labs(title="2R") +
      theme_nothing() +
      theme(text=element_text(size=fontsize-2),
            plot.title=element_text(size=fontsize+2,
                                    margin=margin(b=12)))
print(inversion_2R)

pi_between_2R <- diversity_2R %>%
   filter(Scaffold == "2R",
          Statistic %in% c("D_jn", "D_jk", "D_YS")) %>%
   mutate(Statistic=factor(Statistic, levels=c("D_YS", "D_jn", "D_jk"))) %>%
   ggplot(aes(x=Position, y=Value, colour=Statistic)) +
      geom_line() +
      geom_hline(yintercept=pi_genomewide[["Dyak"]], linetype=2, alpha=0.5) +
      scale_colour_brewer(palette="Set2",
                          labels=c(expression(d[~S-Y]),
                                   expression(pi[~j-n]),
                                   expression(pi[~j-k]))) +
      labs(x="Position (Mb)", y=expression(pi)) +
      guides(fill=FALSE, colour=guide_legend(override.aes=list(size=3))) +
      theme_bw() +
      ylim(0, ymax[["pi_between"]]) +
      scale_x_continuous(limits=bounds[["2R"]],
                         expand=expand_scale(mult=c(0.05, 0.06))) +
      theme(plot.title=element_text(size=fontsize, face="bold"),
            axis.text=element_text(size=fontsize),
            axis.title=element_text(size=fontsize, face="bold"),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            #legend.title=element_text(size=fontsize, face="bold"),
            legend.title=element_blank(),
            legend.text=element_text(size=fontsize-2),
            legend.text.align=0,
            legend.position="top",
            legend.spacing.x=unit(0.3, "cm"),
            legend.box.margin=margin(t=-0.1, r=0,
                                     b=-0.3, l=0, unit='cm'),
            strip.text=element_text(size=fontsize),
            plot.margin=margin(t=0.2, r=0.5, b=0.0, l=0.0, unit="cm"))
print(pi_between_2R)

pi_within_2R <- diversity_2R %>%
   filter(Scaffold == "2R",
          Statistic %in% c("pi_2Rj", "pi_2Rk", "pi_2Rn")) %>%
   mutate(Statistic=factor(Statistic, levels=c("pi_2Rj", "pi_2Rn", "pi_2Rk"))) %>%
   ggplot(aes(x=Position, y=Value, colour=Statistic)) +
      geom_line() +
      geom_hline(yintercept=pi_genomewide[["Dyak"]], linetype=2, alpha=0.5) +
      scale_colour_brewer(palette="Set2",
                          labels=c(expression(pi[~j]),
                                   expression(pi[~n]),
                                   expression(pi[~k]))) +
      labs(x="Position (Mb)", y=expression(pi)) +
      guides(fill=FALSE, colour=guide_legend(override.aes=list(size=3))) +
      theme_bw() +
      ylim(0, ymax[["pi_within"]]) +
      scale_x_continuous(limits=bounds[["2R"]],
                         expand=expand_scale(mult=c(0.05, 0.06))) +
      theme(plot.title=element_text(size=fontsize, face="bold"),
            axis.text=element_text(size=fontsize),
            axis.title=element_text(size=fontsize, face="bold"),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            #legend.title=element_text(size=fontsize, face="bold"),
            legend.title=element_blank(),
            legend.text=element_text(size=fontsize-2),
            legend.text.align=0,
            legend.position="top",
            legend.spacing.x=unit(0.3, "cm"),
            legend.box.margin=margin(t=-0.1, r=0,
                                     b=-0.3, l=0, unit='cm'),
            strip.text=element_text(size=fontsize),
            plot.margin=margin(t=0.2, r=0.5, b=0.0, l=0.0, unit="cm"))
print(pi_within_2R)

F_ST_2R <- diversity_2R %>%
   filter(Scaffold == "2R",
          Statistic %in% c("F_ST_jk", "F_ST_jn")) %>%
   mutate(Statistic=factor(Statistic, levels=c("F_ST_jk", "F_ST_jn"))) %>%
   ggplot(aes(x=Position, y=Value, colour=Statistic)) +
      geom_line() +
      scale_colour_brewer(palette="Set2",
                          labels=c(expression(F[ST]^{~j-k}),
                                   expression(F[ST]^{~j-n}))) +
      labs(x="Position (Mb)", y=expression(F[ST])) +
      guides(fill=FALSE, colour=guide_legend(override.aes=list(size=3))) +
      theme_bw() +
      ylim(0, ymax[["F_ST"]]) +
      scale_x_continuous(limits=bounds[["2R"]],
                         expand=expand_scale(mult=c(0.05, 0.06))) +
      theme(plot.title=element_text(size=fontsize, face="bold"),
            axis.text=element_text(size=fontsize),
            axis.title=element_text(size=fontsize, face="bold"),
            axis.title.x=element_text(vjust=-4),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            #legend.title=element_text(size=fontsize, face="bold"),
            legend.title=element_blank(),
            legend.text=element_text(size=fontsize-2),
            legend.text.align=0,
            legend.position="top",
            legend.spacing.x=unit(0.3, "cm"),
            legend.box.margin=margin(t=-0.1, r=0,
                                     b=-0.3, l=0, unit='cm'),
            strip.text=element_text(size=fontsize),
            plot.margin=margin(t=0.2, r=0.5, b=1.0, l=0.0, unit="cm"))
print(F_ST_2R)

#Now combine all of these into one figure with cowplot:
column_2L_plots <- align_plots(inversion_2L, pi_between_2L, pi_within_2L, F_ST_2L,
                               align="v",
                               axis="lr")
column_2R_plots <- align_plots(inversion_2R, pi_between_2R, pi_within_2R, F_ST_2R,
                               align="v",
                               axis="lr")
Fig2_inversions <- plot_grid(column_2L_plots[[1]], column_2R_plots[[1]],
                             align="h",
                             axis="b",
                             nrow=1)
Fig2_pi_between <- plot_grid(column_2L_plots[[2]], column_2R_plots[[2]],
                             align="h",
                             axis="bt",
                             nrow=1)
Fig2_pi_within <- plot_grid(column_2L_plots[[3]], column_2R_plots[[3]],
                            align="h",
                            axis="bt",
                            nrow=1)
Fig2_F_ST <- plot_grid(column_2L_plots[[4]], column_2R_plots[[4]],
                       align="h",
                       axis="bt",
                       nrow=1)
Fig2_diversity <- plot_grid(Fig2_inversions, Fig2_pi_between, Fig2_pi_within, Fig2_F_ST,
                            ncol=1,
                            rel_heights=c(0.35, 1, 1, 1),
                            labels=c("", "a", "b", "c"))
#Fig2_diversity <- plot_grid(inversion_2L, inversion_2R,
#                            pi_between_2L, pi_between_2R,
#                            pi_within_2L, pi_within_2R,
#                            F_ST_2L, F_ST_2R,
#                            ncol=2,
#                            rel_heights=c(0.25, 1, 1, 1),
#                            labels=c("", "", "a", "", "b", "", "c", ""))
print(Fig2_diversity)
save_plot('Fig2_InversionDiversityDivergence.pdf',
       plot=Fig2_diversity,
       base_width=16.0/2.54,
       base_height=8.0/2.54,
       ncol=2,
       nrow=4,
       dpi=500)

#Dmel plots:
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In2Lt_pi_S_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "pi_2Lt_S", "UsableFraction_pi_2Lt_S"))
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In2Lt_pi_I_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "pi_2Lt_I", "UsableFraction_pi_2Lt_I")) %>%
   full_join(diversity_Dmel, by=c("Scaffold", "WindowStart"))
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In2Lt_D_SI_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "D_2Lt", "UsableFraction_D_2Lt")) %>%
   full_join(diversity_Dmel, by=c("Scaffold", "WindowStart"))
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In2RNS_pi_S_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "pi_2RNS_S", "UsableFraction_pi_2RNS_S")) %>%
   full_join(diversity_Dmel, by=c("Scaffold", "WindowStart"))
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In2RNS_pi_I_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "pi_2RNS_I", "UsableFraction_pi_2RNS_I")) %>%
   full_join(diversity_Dmel, by=c("Scaffold", "WindowStart"))
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In2RNS_D_SI_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "D_2RNS", "UsableFraction_D_2RNS")) %>%
   full_join(diversity_Dmel, by=c("Scaffold", "WindowStart"))
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In3LOk_pi_S_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "pi_3LOk_S", "UsableFraction_pi_3LOk_S")) %>%
   full_join(diversity_Dmel, by=c("Scaffold", "WindowStart"))
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In3LOk_pi_I_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "pi_3LOk_I", "UsableFraction_pi_3LOk_I")) %>%
   full_join(diversity_Dmel, by=c("Scaffold", "WindowStart"))
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In3LOk_D_SI_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "D_3LOk", "UsableFraction_D_3LOk")) %>%
   full_join(diversity_Dmel, by=c("Scaffold", "WindowStart"))
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In3RK_pi_S_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "pi_3RK_S", "UsableFraction_pi_3RK_S")) %>%
   full_join(diversity_Dmel, by=c("Scaffold", "WindowStart"))
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In3RK_pi_I_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "pi_3RK_I", "UsableFraction_pi_3RK_I")) %>%
   full_join(diversity_Dmel, by=c("Scaffold", "WindowStart"))
diversity_Dmel <- read.table(gzfile('Dmel_DPGP3_Zambia_In3RK_D_SI_f0.5_4fold_w100kb.tsv.gz'),
                             header=FALSE,
                             colClasses=diversity_col_classes,
                             col.names=c("Scaffold", "WindowStart", "D_3RK", "UsableFraction_D_3RK")) %>%
   full_join(diversity_Dmel, by=c("Scaffold", "WindowStart"))

#Make Position into Mb:
diversity_Dmel <- diversity_Dmel %>%
   mutate(Position=WindowStart/1e6)

#Hudson, Slatkin, Maddison F_ST estimator: 1-H_w/H_b
#Their H_w is the average of within-pop pi
#Their H_b is D_xy between pops
diversity_Dmel <- diversity_Dmel %>%
   mutate(F_ST_2Lt=1-(pi_2Lt_S+pi_2Lt_I)/(D_2Lt*2),
          F_ST_2RNS=1-(pi_2RNS_S+pi_2RNS_I)/(D_2RNS*2),
          F_ST_3LOk=1-(pi_3LOk_S+pi_3LOk_I)/(D_3LOk*2),
          F_ST_3RK=1-(pi_3RK_S+pi_3RK_I)/(D_3RK*2))

#Filter for usable fraction of sites in a window, and pivot df:
diversity_Dmel_2L <- diversity_Dmel %>%
   filter(UsableFraction_pi_2Lt_S >= min_sites_thresh,
          UsableFraction_pi_2Lt_I >= min_sites_thresh,
          UsableFraction_D_2Lt >= min_sites_thresh) %>%
   gather("Statistic", "Value", -Scaffold, -WindowStart, -Position)

#Filter for usable fraction of sites in a window, and pivot df:
diversity_Dmel_2R <- diversity_Dmel %>%
   filter(UsableFraction_pi_2RNS_S >= min_sites_thresh,
          UsableFraction_pi_2RNS_I >= min_sites_thresh,
          UsableFraction_D_2RNS >= min_sites_thresh) %>%
   gather("Statistic", "Value", -Scaffold, -WindowStart, -Position)

#Filter for usable fraction of sites in a window, and pivot df:
diversity_Dmel_3L <- diversity_Dmel %>%
   filter(UsableFraction_pi_3LOk_S >= min_sites_thresh,
          UsableFraction_pi_3LOk_I >= min_sites_thresh,
          UsableFraction_D_3LOk >= min_sites_thresh) %>%
   gather("Statistic", "Value", -Scaffold, -WindowStart, -Position)

#Filter for usable fraction of sites in a window, and pivot df:
diversity_Dmel_3R <- diversity_Dmel %>%
   filter(UsableFraction_pi_3RK_S >= min_sites_thresh,
          UsableFraction_pi_3RK_I >= min_sites_thresh,
          UsableFraction_D_3RK >= min_sites_thresh) %>%
   gather("Statistic", "Value", -Scaffold, -WindowStart, -Position)

#Make a list of known inversion breakpoints in Dmel:
Breakpoints_Dmel <- list("2Lt"=data.frame(x=2.225744, xend=13.154180, y=0.55, yend=0.60),
                         "2RNS"=data.frame(x=15.391154, xend=20.276334, y=0.55, yend=0.60),
                         "3LOk"=data.frame(),
                         "3RK"=data.frame(x=11.750567, xend=26.140370, y=0.55, yend=0.60))
bounds_Dmel <- list("2L"=c(1,23513712)/1e6,
                    "2R"=c(1,25286936)/1e6,
                    "3L"=c(1,28110227)/1e6,
                    "3R"=c(1,32079331)/1e6)
#Now make the individual F_ST plots so we can combine them in the end:
F_ST_2Lt <- diversity_Dmel_2L %>%
   filter(Scaffold == "2L",
          Statistic %in% c("F_ST_2Lt")) %>%
   ggplot(aes(x=Position, y=Value, colour=Statistic)) +
   geom_line() +
   geom_rect(data=Breakpoints_Dmel[["2Lt"]], aes(xmin=x, xmax=xend, ymin=y, ymax=yend), inherit.aes=FALSE) +
   geom_text(data=Breakpoints_Dmel[["2Lt"]], aes(x=xend+3, y=y+0.025, label="In(2L)t"), inherit.aes=FALSE, size=fontsize/4) +
   scale_colour_brewer(palette="Set2",
                       labels=c(expression(F[ST]^{~In2Lt}))) +
   labs(x="Position (Mb)", y=expression(F[ST])) +
   guides(fill=FALSE, colour=guide_legend(override.aes=list(size=3))) +
   theme_bw() +
   ylim(0, ymax[["F_ST"]]) +
   xlim(bounds_Dmel[["2L"]]) +
   theme(plot.title=element_text(size=fontsize, face="bold"),
         axis.text=element_text(size=fontsize),
         axis.title=element_text(size=fontsize, face="bold"),
         axis.title.x=element_text(vjust=-4),
         axis.title.y=element_text(vjust=4),
         #legend.title=element_text(size=fontsize, face="bold"),
         legend.title=element_blank(),
         legend.text=element_text(size=fontsize-2),
         legend.text.align=0,
         legend.position="top",
         legend.spacing.x=unit(0.3, "cm"),
         legend.box.margin=margin(t=-0.1, r=0,
                                  b=-0.3, l=0, unit='cm'),
         strip.text=element_text(size=fontsize),
         plot.margin=margin(t=0.2, r=0.5, b=1.0, l=1.0, unit="cm"))
print(F_ST_2Lt)

F_ST_2RNS <- diversity_Dmel_2R %>%
   filter(Scaffold == "2R",
          Statistic %in% c("F_ST_2RNS")) %>%
   ggplot(aes(x=Position, y=Value, colour=Statistic)) +
   geom_line() +
   geom_rect(data=Breakpoints_Dmel[["2RNS"]], aes(xmin=x, xmax=xend, ymin=y, ymax=yend), inherit.aes=FALSE) +
   geom_text(data=Breakpoints_Dmel[["2RNS"]], aes(x=x-3, y=y+0.025, label="In(2R)NS"), inherit.aes=FALSE, size=fontsize/4) +
   scale_colour_brewer(palette="Set2",
                       labels=c(expression(F[ST]^{~In2RNS}))) +
   labs(x="Position (Mb)", y=expression(F[ST])) +
   guides(fill=FALSE, colour=guide_legend(override.aes=list(size=3))) +
   theme_bw() +
   ylim(0, ymax[["F_ST"]]) +
   xlim(bounds_Dmel[["2R"]]) +
   theme(plot.title=element_text(size=fontsize, face="bold"),
         axis.text=element_text(size=fontsize),
         axis.title=element_text(size=fontsize, face="bold"),
         axis.title.x=element_text(vjust=-4),
         axis.title.y=element_text(vjust=4),
         #legend.title=element_text(size=fontsize, face="bold"),
         legend.title=element_blank(),
         legend.text=element_text(size=fontsize-2),
         legend.text.align=0,
         legend.position="top",
         legend.spacing.x=unit(0.3, "cm"),
         legend.box.margin=margin(t=-0.1, r=0,
                                  b=-0.3, l=0, unit='cm'),
         strip.text=element_text(size=fontsize),
         plot.margin=margin(t=0.2, r=0.5, b=1.0, l=1.0, unit="cm"))
print(F_ST_2RNS)

F_ST_3LOk <- diversity_Dmel_3L %>%
   filter(Scaffold == "3L",
          Statistic %in% c("F_ST_3LOk")) %>%
   ggplot(aes(x=Position, y=Value, colour=Statistic)) +
   geom_line() +
#   geom_rect(data=Breakpoints_Dmel[["3LOk"]], aes(xmin=x, xmax=xend, ymin=y, ymax=yend), inherit.aes=FALSE) +
#   geom_text(data=Breakpoints_Dmel[["3LOk"]], aes(x=x-3, y=y+0.025, label="In(3L)Ok"), inherit.aes=FALSE, size=fontsize/4) +
   scale_colour_brewer(palette="Set2",
                       labels=c(expression(F[ST]^{~In3LOk}))) +
   labs(x="Position (Mb)", y=expression(F[ST])) +
   guides(fill=FALSE, colour=guide_legend(override.aes=list(size=3))) +
   theme_bw() +
   ylim(0, ymax[["F_ST"]]) +
   xlim(bounds_Dmel[["3L"]]) +
   theme(plot.title=element_text(size=fontsize, face="bold"),
         axis.text=element_text(size=fontsize),
         axis.title=element_text(size=fontsize, face="bold"),
         axis.title.x=element_text(vjust=-4),
         axis.title.y=element_text(vjust=4),
         #legend.title=element_text(size=fontsize, face="bold"),
         legend.title=element_blank(),
         legend.text=element_text(size=fontsize-2),
         legend.text.align=0,
         legend.position="top",
         legend.spacing.x=unit(0.3, "cm"),
         legend.box.margin=margin(t=-0.1, r=0,
                                  b=-0.3, l=0, unit='cm'),
         strip.text=element_text(size=fontsize),
         plot.margin=margin(t=0.2, r=0.5, b=1.0, l=1.0, unit="cm"))
print(F_ST_3LOk)

F_ST_3RK <- diversity_Dmel_3R %>%
   filter(Scaffold == "3R",
          Statistic %in% c("F_ST_3RK")) %>%
   ggplot(aes(x=Position, y=Value, colour=Statistic)) +
   geom_line() +
   geom_rect(data=Breakpoints_Dmel[["3RK"]], aes(xmin=x, xmax=xend, ymin=y, ymax=yend), inherit.aes=FALSE) +
   geom_text(data=Breakpoints_Dmel[["3RK"]], aes(x=x-3, y=y+0.025, label="In(3R)K"), inherit.aes=FALSE, size=fontsize/4) +
   scale_colour_brewer(palette="Set2",
                       labels=c(expression(F[ST]^{~In3RK}))) +
   labs(x="Position (Mb)", y=expression(F[ST])) +
   guides(fill=FALSE, colour=guide_legend(override.aes=list(size=3))) +
   theme_bw() +
   ylim(0, ymax[["F_ST"]]) +
   xlim(bounds_Dmel[["3R"]]) +
   theme(plot.title=element_text(size=fontsize, face="bold"),
         axis.text=element_text(size=fontsize),
         axis.title=element_text(size=fontsize, face="bold"),
         axis.title.x=element_text(vjust=-4),
         axis.title.y=element_text(vjust=4),
         #legend.title=element_text(size=fontsize, face="bold"),
         legend.title=element_blank(),
         legend.text=element_text(size=fontsize-2),
         legend.text.align=0,
         legend.position="top",
         legend.spacing.x=unit(0.3, "cm"),
         legend.box.margin=margin(t=-0.1, r=0,
                                  b=-0.3, l=0, unit='cm'),
         strip.text=element_text(size=fontsize),
         plot.margin=margin(t=0.2, r=0.5, b=1.0, l=1.0, unit="cm"))
print(F_ST_3RK)

FigSZ_Dmel_divergence <- plot_grid(F_ST_2Lt, F_ST_2RNS, F_ST_3LOk, F_ST_3RK,
                                   align="hv", axis="tblr", nrow=2, ncol=2,
                                   labels="auto")
print(FigSZ_Dmel_divergence)
save_plot('FigSZ_DmelDivergence.pdf', FigSZ_Dmel_divergence,
          nrow=2, ncol=2, base_width=16.0/2.54, base_height=12.0/2.54,
          dpi=500)
