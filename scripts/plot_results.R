library(dplyr)
library(ggplot2)

std <- function(x) sd(x)/sqrt(length(x))

simplex_results = read.csv(snakemake@input[[1]])

summary_df = simplex_results %>%
  filter(d == 2, inferred_K) %>%
  group_by(method, data, var_latents, noise_factor) %>%
  summarise(nruns = n(),
            ari_mean = mean(ari),
            ari_se = std(ari),
            true_k_mean = mean(true_k),
            estimated_k_mean = mean(estimated_k),
            estimated_k_se = std(estimated_k), .groups="keep") %>%
  mutate(data = recode(data, x = 'observed', z = 'latent'))
print(table(summary_df$nruns))
g = ggplot(summary_df, aes(colour=method, y=ari_mean, x=noise_factor)) +
  geom_line() +
  geom_errorbar(aes(ymin=ari_mean - ari_se,
                    ymax=ari_mean + ari_se),
                width=0.2) +
  facet_grid(" var_latents ~ data" )

ggsave(snakemake@output[[1]], g, width=10, height=12, units="in")
