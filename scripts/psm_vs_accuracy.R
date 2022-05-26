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
               same_true_cluster = (i_true_cluster == j_true_cluster) * 1)
    return(psm_vs_acc_df)
}

U_values = c(5, 7, 10, 13)
N_values = c(4, 9, 12, 16)
seeds = 1:10
runs = do.call(paste0, expand.grid("simplex_results/2_", N_values, "_", U_values, "/", seeds, "/results.rda"))

combined_psm_vs_acc_df = do.call(rbind, lapply(runs, construct_psm_acc_df))

se <- function(x) sqrt(var(x) / length(x))

summary_df_U_N = combined_psm_vs_acc_df %>%
    mutate(psm_bin = cut(PSM, breaks=seq(0, 1, by=0.1)),
           U = as.integer(str_extract(rda_file, "(?<=/2_)\\d+(?=_)")),
           N = as.integer(str_extract(rda_file, "(?<=/2_\\d{1,2}_)\\d+(?=/)")),
           U_N = str_extract(rda_file, "2_(\\d+)_(\\d+)")) %>%
    drop_na(psm_bin) %>%
    group_by(psm_bin, U, N) %>%
    summarise(mean_acc = mean(same_true_cluster),
              se_acc = se(same_true_cluster), .groups="keep")

ggplot(summary_df_U_N,
       aes(x=psm_bin, y=mean_acc, colour=U)) +
    geom_point() +
    facet_grid(N ~ .) +
    geom_errorbar(aes(ymin=mean_acc-1.96*se_acc, ymax=mean_acc+1.96*se_acc), width=.2) +
    geom_abline(slope=0.1, intercept=-0.05, colour="green") +
    labs(x="Binned posteriror similarity score",
         y="Proportion of pairs in same true cluster") +
    scale_y_continuous(limits=c(0, 1))
ggsave("plots/psm_vs_acc_U_N.png", width=10, height=6, units="in")

summary_df = combined_psm_vs_acc_df %>%
    mutate(psm_bin = cut(PSM, breaks=seq(0, 1, by=0.1))) %>%
    drop_na(psm_bin) %>%
    group_by(psm_bin) %>%
    summarise(mean_acc = mean(same_true_cluster),
              se_acc = se(same_true_cluster), .groups="keep")

ggplot(summary_df,
       aes(x=psm_bin, y=mean_acc)) +
    geom_point() +
    geom_errorbar(aes(ymin=mean_acc-1.96*se_acc, ymax=mean_acc+1.96*se_acc), width=.2) +
    geom_abline(slope=0.1, intercept=-0.05, colour="green") +
    labs(x="Binned posteriror similarity score",
         y="Proportion of pairs in same true cluster") +
    scale_y_continuous(limits=c(0, 1))
ggsave("plots/psm_vs_acc.png", width=6, height=6, units="in")
