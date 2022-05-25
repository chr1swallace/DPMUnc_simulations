library(ggplot2)

n_simplex <- function(n) {
  extra_point = (1 + sqrt(n+1))/n * rep(1, n)
  basis_elements = lapply(1:n, function(x) { y = rep(0, n); y[x] = 1; y})
  all_points = basis_elements
  all_points[[n+1]] = extra_point

  df = data.frame(do.call(rbind, all_points))
  colnames(df) = paste0("dim", 1:n)
  return (df)
}

generate_basic_uncertain_data <- function(n, d, k, var_latents, noise_factor, group_means) {
  classes = sample(1:k, n, replace=TRUE)
  simulation_df = data.frame(class = classes, sample_id = paste0("sample_", 1:n))

  obsVars_df = data.frame(matrix(rchisq(n*d,1), nrow=n, ncol=d)) * noise_factor
  colnames(obsVars_df) = paste0("sigmasq", 1:d)
  simulation_df = cbind(simulation_df, obsVars_df)

  colnames(group_means) = paste0("mu", 1:d)
  group_means$class = 1:k

  simulation_df = merge(simulation_df, group_means, by="class")

  latents = matrix(rnorm(n*d,
                         mean=as.vector(as.matrix(simulation_df[, paste0("mu", 1:d)])),
                         sd=sqrt(var_latents)),
                   nrow=n, ncol=d)
  latents = data.frame(latents)
  colnames(latents) = paste0("z", 1:d)

  simulation_df = cbind(simulation_df, latents)

  observed = matrix(rnorm(n*d,
                          mean=as.vector(as.matrix(latents)),
                          sd=sqrt(as.vector(as.matrix(obsVars_df)))),
                    nrow=n, ncol=d)
  observed = data.frame(observed)
  colnames(observed) = paste0("x", 1:d)

  simulation_df = cbind(simulation_df, observed)

  rownames(simulation_df) = simulation_df$sample_id
  head(simulation_df)

  obsData = as.matrix(simulation_df[, paste0("x", 1:d)])
  obsVars = as.matrix(simulation_df[, paste0("sigmasq", 1:d)])

  pca = prcomp(obsData)
  print(summary(pca))
  pca_df = data.frame(pca$x)
  colnames(pca_df) = paste0("xPC", 1:d)

  simulation_df = cbind(simulation_df, pca_df)

  head(simulation_df)

  return (list(df=simulation_df,
               pca=pca,
               obsVars=obsVars,
               obsData=obsData))
}

d = 2
group_means = data.frame(n_simplex(d))

n = d * 10
k = d + 1
set.seed(1314)
combinations = expand.grid(noise_factor=c(5, 7, 10, 13), var_latents=c(4, 9, 12, 16))
all_generated = do.call(rbind, apply(combinations, MARGIN=1, FUN=function(x) {
  df = generate_basic_uncertain_data(n=n, d=d, k=k, var_latents=x[['var_latents']],
                                     noise_factor=x[['noise_factor']],
                                     group_means=group_means * 10)$df
  df$N = x[['var_latents']]
  df$U = x[['noise_factor']]
  df
}))

cbbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
g = ggplot(all_generated, aes(x=z1, y=z2, colour=factor({class}))) +
  geom_point(size=2, shape=1) +
  geom_point(size=2, mapping=aes(x=x1, y=x2), alpha=0.5) +
  geom_segment(aes(xend=x1, yend=x2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  scale_colour_manual(values=cbbPalette) +
  geom_point(size=1, shape=2, mapping=aes(x=mu1, y=mu2)) +
#  stat_ellipse() +
  facet_grid("N ~ U", labeller = label_both) +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        text=element_text(size=15))
ggsave(snakemake@output[[1]], width=6, height=6, units="in", dpi=600)
