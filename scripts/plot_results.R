library(dplyr)
library(ggplot2)

std <- function(x) sd(x)/sqrt(length(x))

simplex_results = read.csv("simplex_results/full_results.csv") #snakemake@input[[1]])
#simplex_results = read.csv(snakemake@input[[1]])
cbbPalette = c("#0072B2", "#D55E00", "#009E73", "#E69F00", "#00FFFF")

summary_df = simplex_results %>%
  filter(d == 2, inferred_K) %>%
  group_by(method, data, var_latents, noise_factor) %>%
  summarise(nruns = n(),
            ari_mean = mean(ari),
            ari_se = std(ari),
            true_k_mean = mean(true_k),
            estimated_k_mean = mean(estimated_k),
            estimated_k_se = std(estimated_k), .groups="keep") %>%
  mutate(data = recode(data, x = 'observed', z = 'latent'),
         method = recode(method, DPMUnc_novar = "DPMZeroUnc"),
         U = noise_factor,
         N = var_latents)
print(table(summary_df$nruns))
g = ggplot(summary_df, aes(colour=method, y=ari_mean, x=U)) +
  geom_line() +
  scale_colour_manual(values=cbbPalette) +
  geom_errorbar(aes(ymin=ari_mean - ari_se,
                    ymax=ari_mean + ari_se),
                width=0.2) +
  labs(y="Mean accuracy of clustering (ARI)",
       colour="Method") +
  theme(text=element_text(size=15),
        legend.position = "bottom") +
  scale_x_continuous(breaks = unique(summary_df$U)) +
  facet_grid(" N ~ data", labeller = label_both)

ggsave(snakemake@output[[1]], g, width=6, height=8, units="in")
