#Libraries necessary for processing and plotting:
library(tidyverse)
library(cowplot)

#Column names for the output of abcreg:
abc_2L_colnames <- c("N_anc", "N_STD", "N_INV", "In2L")
abc_2Rjn_colnames <- c("N_anc", "N_n", "N_j", "N_S", "In2Rn", "In2Rj", "Dsan")
abc_2Rk_colnames <- c("N_anc", "N_j", "N_k", "In2Rk")

#Note: We need to downscale 2L and 2Rk times by 10 generations per year
generations_per_year <- 10

#Read in the abcreg outputs:
ABC <- read.table(gzfile('Dyak_2Rjn_4fold_finalgenomes_priorset2_abcreg_t1e-3.0.tangent.post.gz'),
                  header=FALSE,
                  col.names=abc_2Rjn_colnames) %>%
   mutate(Breakpoint="p2Rn") %>%
   gather("Parameter", "Value", -Breakpoint)
ABC <- rbind(ABC, read.table(gzfile('Dyak_2Rjn_4fold_finalgenomes_priorset2_abcreg_t1e-3.1.tangent.post.gz'),
                             header=FALSE,
                             col.names=abc_2Rjn_colnames) %>%
                mutate(Breakpoint="d2Rn") %>%
                gather("Parameter", "Value", -Breakpoint))
ABC <- rbind(ABC, read.table(gzfile('Dyak_2Rjk_4fold_finalgenomes_priorset3_t1e-3.0.tangent.post.gz'),
                             header=FALSE,
                             col.names=abc_2Rk_colnames) %>%
                mutate(Breakpoint="p2Rk", In2Rk=In2Rk/generations_per_year) %>%
                gather("Parameter", "Value", -Breakpoint))
ABC <- rbind(ABC, read.table(gzfile('Dyak_2Rjk_4fold_finalgenomes_priorset3_t1e-3.1.tangent.post.gz'),
                             header=FALSE,
                             col.names=abc_2Rk_colnames) %>%
                mutate(Breakpoint="d2Rk", In2Rk=In2Rk/generations_per_year) %>%
                gather("Parameter", "Value", -Breakpoint))
ABC <- rbind(ABC, read.table(gzfile('Dyak_2L_4fold_finalgenomes_priorset3_t1e-3.0.tangent.post.gz'),
                             header=FALSE,
                             col.names=abc_2L_colnames) %>%
                mutate(Breakpoint="d2L", In2L=In2L/generations_per_year) %>%
                gather("Parameter", "Value", -Breakpoint))
ABC <- rbind(ABC, read.table(gzfile('Dyak_2L_4fold_finalgenomes_priorset3_t1e-3.1.tangent.post.gz'),
                             header=FALSE,
                             col.names=abc_2L_colnames) %>%
                mutate(Breakpoint="p2L", In2L=In2L/generations_per_year) %>%
                gather("Parameter", "Value", -Breakpoint))

#Now we filter on breakpoints, and plot:
ABC %>% filter(Breakpoint %in% c("p2L", "d2Rn", "d2Rk"),
               Parameter %in% c("In2Rj", "In2Rn", "In2Rk", "In2L", "Dsan")) %>%
   mutate(Parameter=factor(Parameter,
                           levels=c("In2L", "In2Rk", "Dsan", "In2Rj", "In2Rn"))) %>%
   ggplot(aes(x=Parameter, y=Value)) +
      geom_violin(aes(fill=Parameter), scale="width") +
      geom_boxplot(alpha=0.5, width=0.07) +
      coord_flip() +
      guides(fill=FALSE) +
      scale_fill_brewer(palette="Set2") +
      scale_x_discrete(labels=c("In(2L)", "2Rk", "Dsan split", "2Rj", "2Rn")) +
      theme_bw() +
      labs(x="Inversion", y="Age (years)") +
      theme(text=element_text(size=12),
            axis.title=element_text(size=12, face="bold"),
            axis.text=element_text(size=12))
ggsave('Fig3_Inversion_ABC.pdf', width=16.0, height=12.0, units="cm", dpi=500)
