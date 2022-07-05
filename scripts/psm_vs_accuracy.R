library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

construct_psm_acc_df <- function(rda_file) {
    load(rda_file)
    psm_vs_acc_df = result$bigpsm %>%
        data.frame() %>%
        setNames(1:30) %>%
        mutate(i=1:30,
               i_true_cluster=simulation$df$class) %>%
        pivot_longer(cols=-c(i, i_true_cluster), names_to="j", values_to="PSM") %>%
        mutate(j_true_cluster = rep(simulation$df$class, 30),
               rda_file = rda_file,
               same_true_cluster = (i_true_cluster == j_true_cluster) * 1) %>%
        filter(i < j)
    return(psm_vs_acc_df)
}

U_values = c(5, 7, 10, 13)
N_values = c(4, 9, 12, 16)
seeds = 0:99
runs = do.call(paste0, expand.grid("simplex_results/2_", N_values, "_", U_values, "/", seeds, "/results.rda"))

combined_psm_vs_acc_df = do.call(rbind, lapply(runs, construct_psm_acc_df))
combined_psm_vs_acc_df = combined_psm_vs_acc_df %>%
    mutate(psm_bin = cut(PSM, breaks=seq(0, 1, by=0.1), include.lowest=TRUE),
           psm_bin_mid = as.numeric(psm_bin) * 0.1 - 0.05,
           U = as.integer(str_extract(rda_file, "(?<=/2_)\\d+(?=_)")),
           N = as.integer(str_extract(rda_file, "(?<=/2_\\d{1,2}_)\\d+(?=/)")),
           U_N = str_extract(rda_file, "2_(\\d+)_(\\d+)"))

se <- function(x) sqrt(var(x) / length(x))

summary_df_U_N = combined_psm_vs_acc_df %>%
    group_by(psm_bin, U, N) %>%
    summarise(mean_acc = mean(same_true_cluster),
              count = n(),
              se_acc = se(same_true_cluster), .groups="keep")

#rmsq_df = summary_df_U_N %>%
#   mutate(psm_bin_mid = as.numeric(psm_bin) * 0.1 - 0.05,
#          rmsq = (mean_acc - psm_bin_mid)**2) %>%
#    group_by(U, N) %>%
#    summarise(rmsq_total = sum(rmsq))

g = ggplot(summary_df_U_N,
       aes(x=psm_bin, y=mean_acc)) +
    geom_point(aes(size=count)) +
    facet_grid(N ~ U, labeller = label_both) +
#    geom_errorbar(aes(ymin=mean_acc-1.96*se_acc, ymax=mean_acc+1.96*se_acc), width=.2) +
    geom_abline(slope=0.1, intercept=-0.05, colour="green") +
    theme(axis.text=element_text(size=8),
          axis.text.x= element_text(angle = 90, size=8, vjust=0.5)) +
    labs(x="Binned posterior similarity score",
         size="Number of pairs",
         y="Proportion of pairs in same true cluster") +
    scale_y_continuous(limits=c(0, 1))
ggsave("plots/psm_vs_acc_U_N_faceted.pdf", g, width=8, height=4, units="in")

g = ggplot(summary_df_U_N,
       aes(x=psm_bin, y=mean_acc, colour=U)) +
    geom_point() +
    facet_grid(N ~ .) +
#    geom_errorbar(aes(ymin=mean_acc-1.96*se_acc, ymax=mean_acc+1.96*se_acc), width=.2) +
    geom_abline(slope=0.1, intercept=-0.05, colour="green") +
    labs(x="Binned posterior similarity score",
         y="Proportion of pairs in same true cluster") +
    scale_y_continuous(limits=c(0, 1))
ggsave("plots/psm_vs_acc_U_N.pdf", g, width=10, height=6, units="in")

summary_df = combined_psm_vs_acc_df %>%
    group_by(psm_bin) %>%
    summarise(mean_acc = mean(same_true_cluster),
              se_acc = se(same_true_cluster), .groups="keep")

g = ggplot(summary_df,
       aes(x=psm_bin, y=mean_acc)) +
#    geom_errorbar(aes(ymin=mean_acc-1.96*se_acc, ymax=mean_acc+1.96*se_acc), width=.2) +
    geom_abline(slope=0.1, intercept=-0.05, colour="green4", size=2) +
    geom_point(size=4) +
    labs(x="Binned posterior similarity score",
         y="Proportion of pairs in same true cluster") +
    scale_y_continuous(limits=c(0, 1))
ggsave("plots/psm_vs_acc.pdf", g, width=6, height=6, units="in")

g = ggplot(combined_psm_vs_acc_df, aes(x=PSM)) +
    geom_histogram(bins=10) +
    labs(x = "Posterior similarity",
         y = "") +
    theme(axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank()) +
    facet_grid(N ~ U, labeller = label_both)
ggsave("plots/psm_hist.pdf", g, width=8, height=4, units="in")
