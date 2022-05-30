library(cluster)
library(dplyr)
library(ggplot2)
library(mclust)
library(mcclust)
library(stringr)
library(tidyr)

library(DPMUnc)

source("scripts/utils.R")

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

seed=13
d=2
var_latents=2
noise_factor=16
k = 1
n = 6

# These settings look good
set.seed(13)
group_means = data.frame(n_simplex(d))
simulation = generate_basic_uncertain_data(n=n, d=d, k=k, var_latents=var_latents, noise_factor=noise_factor,
                                           group_means=group_means[1, ])

cbbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
g = ggplot(simulation$df, aes(x=z1, y=z2, colour=factor({class}))) +
  geom_point(size=2, shape=1) +
  geom_point(mapping=aes(x=x1, y=x2, size=sigmasq2), alpha=0.5) +
  geom_segment(aes(xend=x1, yend=x2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  scale_colour_manual(values=cbbPalette) +
  geom_point(size=1, shape=2, mapping=aes(x=mu1, y=mu2)) +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        text=element_text(size=15))
ggsave("plots/mean_shift.png", width=3, height=3, units="in", dpi=600)

# Set hyperparameters for model
alpha0 = 2; kappa0 = 0.5; beta0 = 0.2 * mean(apply(simulation$obsData, 2, var))

gap_stat <- clusGap(simulation$obsData, FUN = kmeans, nstart = 25, K.max = k * 2, B = 50)
best_K = which.max(gap_stat$Tab[, 3])
kmeans_solution = kmeans(simulation$obsData, centers=best_K, nstart=25)

mclust_solution <- Mclust(simulation$obsData, x=mclustBIC(simulation$obsData))

outputdir = paste0("mean_shift/", str_replace(as.character(as.POSIXlt(Sys.time())), " ", "_"), "/")
dir.create(outputdir, recursive = TRUE)
print(dim(simulation$obsData))
DPMUnc(simulation$obsData, simulation$obsVars, saveFileDir = outputdir, seed=seed,
       kappa0=kappa0, alpha0=alpha0, beta0=beta0,
       saveLatentObs = TRUE,
       nIts=10000, scaleData=FALSE)
result = calc_psms(outputdir)
calls=maxpear(result$bigpsm, method="comp")

outputdir_novar = paste0("mean_shift/", str_replace(as.character(as.POSIXlt(Sys.time())), " ", "_"), "/novar/")
dir.create(outputdir_novar, recursive = TRUE)
DPMUnc(simulation$obsData, simulation$obsVars / 1000000000, saveFileDir = outputdir_novar, seed=seed,
       kappa0=kappa0, alpha0=alpha0, beta0=beta0,
       saveLatentObs = TRUE,
       nIts=10000, scaleData=FALSE)
result_novar = calc_psms(outputdir_novar)
calls_novar=maxpear(result_novar$bigpsm, method="comp")

psm_links_df = result$bigpsm %>%
    data.frame() %>%
    setNames(paste0("sample_", 1:6)) %>%
    mutate(i=paste0("sample_", 1:6)) %>%
    pivot_longer(cols=-c(i), names_to="j", values_to="PSM") %>%
    filter(i < j) %>%
    merge(simulation$df %>% select(sample_id, x1, x2), by.x="i", by.y="sample_id") %>%
    merge(simulation$df %>% select(sample_id, x1, x2), by.x="j", by.y="sample_id", suffixes=c("_from", "_to")) %>%
    mutate(res = "DPMUnc")

psm_links_df_novar = result_novar$bigpsm %>%
    data.frame() %>%
    setNames(paste0("sample_", 1:6)) %>%
    mutate(i=paste0("sample_", 1:6)) %>%
    pivot_longer(cols=-c(i), names_to="j", values_to="PSM") %>%
    filter(i < j) %>%
    merge(simulation$df %>% select(sample_id, x1, x2), by.x="i", by.y="sample_id") %>%
    merge(simulation$df %>% select(sample_id, x1, x2), by.x="j", by.y="sample_id", suffixes=c("_from", "_to")) %>%
    mutate(res = "DPMZeroUnc")

psm_links_df_both = rbind(psm_links_df, psm_links_df_novar)

ggplot(NULL) +
  geom_segment(data=psm_links_df_both, aes(x=x1_from, y=x2_from, xend=x1_to, yend=x2_to, size=PSM, colour=PSM), alpha=0.5) +
  geom_point(data = simulation$df, aes(x=z1, y=z2), size=2, shape=1) +
  geom_point(data = simulation$df, size=2, mapping=aes(x=x1, y=x2), alpha=0.5) +
  geom_segment(data = simulation$df, aes(x=z1, y=z2, xend=x1, yend=x2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = simulation$df, size=1, shape=2, mapping=aes(x=mu1, y=mu2)) +
  scale_size(limits=c(0, 1), guide="none", range=c(0, 1)) +
  scale_colour_distiller(limits=c(0, 1), palette="Reds", direction=1) +
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        text=element_text(size=15)) +
  facet_grid(. ~ res)
ggsave("plots/mean_shift_psm_graph.png", width=6, height=3, units="in", dpi=600)

latents = read_latent_obs(file = paste0(outputdir, "/latentObservations.csv"), d=d)
mean_latents = matrix(colMeans(latents), ncol=d, nrow=n) %>% data.frame() %>% setNames(c("zhat1", "zhat2")) %>% mutate(res = "DPMUnc")

latents_novar = read_latent_obs(file = paste0(outputdir_novar, "/latentObservations.csv"), d=d)
mean_latents_novar = matrix(colMeans(latents_novar), ncol=d, nrow=n) %>% data.frame() %>% setNames(c("zhat1", "zhat2")) %>% mutate(res = "DPMZeroUnc")

mean_latents_both = rbind(mean_latents, mean_latents_novar)

clusterMeans = readClusterParams(file = paste0(outputdir, "/clusterMeans.csv"), nDim=2, nlines_skip=500)
clusterMeans_novar = readClusterParams(file = paste0(outputdir_novar, "/clusterMeans.csv"), nDim=2, nlines_skip=500)

clusterMeans_both = rbind(clusterMeans %>% mutate(res = "DPMUnc"),
						  clusterMeans_novar %>% mutate(res = "DPMZeroUnc"))

empirical_x_means = colMeans(simulation$obsData) %>% t() %>% data.frame()
ggplot(NULL) +
  geom_point(data = simulation$df, aes(x=z1, y=z2), size=2, shape=1) +
  geom_point(data = simulation$df, mapping=aes(x=x1, y=x2, size=sigmasq2), alpha=0.5) +
  geom_segment(data = simulation$df, aes(x=z1, y=z2, xend=x1, yend=x2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = simulation$df, size=1, shape=2, mapping=aes(x=mu1, y=mu2)) +
  geom_point(data = mean_latents_both, mapping=aes(x=zhat1, y=zhat2), colour = "green", shape="square") + 
  geom_segment(data = cbind(simulation$df, mean_latents), aes(x=x1, y=x2, xend=zhat1, yend=zhat2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = empirical_x_means, mapping=aes(x=x1, y=x2), colour="deepskyblue", size=10, shape="x") +
  geom_density_2d(data = clusterMeans_both, mapping=aes(x=X1, y=X2), colour="goldenrod2") +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        text=element_text(size=15)) +
  facet_grid(. ~ res)
ggsave("plots/mean_shift_latents_clustermeans_both.png", width=6, height=3, units="in", dpi=600)

ggplot(NULL) +
  geom_point(data = simulation$df, aes(x=z1, y=z2), size=2, shape=1) +
  geom_point(data = simulation$df, mapping=aes(x=x1, y=x2, size=sigmasq2), alpha=0.5) +
  geom_segment(data = simulation$df, aes(x=z1, y=z2, xend=x1, yend=x2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = simulation$df, size=1, shape=2, mapping=aes(x=mu1, y=mu2)) +
  geom_point(data = mean_latents, mapping=aes(x=zhat1, y=zhat2), colour = "green", shape="square") + 
  geom_segment(data = cbind(simulation$df, mean_latents), aes(x=x1, y=x2, xend=zhat1, yend=zhat2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_density_2d(data = clusterMeans, mapping=aes(x=X1, y=X2), colour="goldenrod2") +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        text=element_text(size=15))
ggsave("plots/mean_shift_latents_clustermeans.png", width=3, height=3, units="in", dpi=600)

ggplot(NULL) +
  geom_point(data = simulation$df, aes(x=z1, y=z2), size=2, shape=1) +
  geom_point(data = simulation$df, mapping=aes(x=x1, y=x2, size=sigmasq2), alpha=0.5) +
  geom_segment(data = simulation$df, aes(x=z1, y=z2, xend=x1, yend=x2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = simulation$df, size=1, shape=2, mapping=aes(x=mu1, y=mu2)) +
  geom_point(data = mean_latents, mapping=aes(x=zhat1, y=zhat2), colour = "green", shape="square") + 
  geom_segment(data = cbind(simulation$df, mean_latents), aes(x=x1, y=x2, xend=zhat1, yend=zhat2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        text=element_text(size=15))
ggsave("plots/mean_shift_latents.png", width=3, height=3, units="in", dpi=600)

ggplot(NULL) +
  geom_point(data = simulation$df, aes(x=z1, y=z2), size=2, shape=1) +
  geom_point(data = simulation$df, mapping=aes(x=x1, y=x2, size=sigmasq2), alpha=0.5) +
  geom_segment(data = simulation$df, aes(x=z1, y=z2, xend=x1, yend=x2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = simulation$df, size=1, shape=2, mapping=aes(x=mu1, y=mu2)) +
  geom_point(data = mean_latents, mapping=aes(x=zhat1, y=zhat2), colour = "green", shape="square") + 
  geom_segment(data = cbind(simulation$df, mean_latents), aes(x=x1, y=x2, xend=zhat1, yend=zhat2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = empirical_x_means, mapping=aes(x=x1, y=x2), colour="deepskyblue", size=10, shape="x") +
  geom_density_2d(data = clusterMeans, mapping=aes(x=X1, y=X2), colour="goldenrod2") +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        text=element_text(size=15))
ggsave("plots/mean_shift_latents_clustermeans_empirical.png", width=3, height=3, units="in", dpi=600)

clusterMeans = readClusterParamsSingle(file = paste0(outputdir, "/clusterMeans.csv"), nDim=2, nlines_skip=500)
clusterMeans_novar = readClusterParamsSingle(file = paste0(outputdir_novar, "/clusterMeans.csv"), nDim=2, nlines_skip=500)

clusterMeans_both = rbind(clusterMeans %>% mutate(res = "DPMUnc"),
						  clusterMeans_novar %>% mutate(res = "DPMZeroUnc"))

empirical_x_means = colMeans(simulation$obsData) %>% t() %>% data.frame()
ggplot(NULL) +
  geom_point(data = simulation$df, aes(x=z1, y=z2), size=2, shape=1) +
  geom_point(data = simulation$df, mapping=aes(x=x1, y=x2, size=sigmasq2), alpha=0.5) +
  geom_segment(data = simulation$df, aes(x=z1, y=z2, xend=x1, yend=x2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = simulation$df, size=1, shape=2, mapping=aes(x=mu1, y=mu2)) +
  geom_point(data = mean_latents_both, mapping=aes(x=zhat1, y=zhat2), colour = "green", shape="square") + 
  geom_segment(data = cbind(simulation$df, mean_latents), aes(x=x1, y=x2, xend=zhat1, yend=zhat2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = empirical_x_means, mapping=aes(x=x1, y=x2), colour="deepskyblue", size=10, shape="x") +
  geom_density_2d(data = clusterMeans_both, mapping=aes(x=X1, y=X2), colour="goldenrod2") +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        text=element_text(size=15)) +
  facet_grid(. ~ res)
ggsave("plots/mean_shift_latents_clustermeans_both_mu1.png", width=6, height=3, units="in", dpi=600)

ggplot(NULL) +
  geom_point(data = simulation$df, aes(x=z1, y=z2), size=2, shape=1) +
  geom_point(data = simulation$df, mapping=aes(x=x1, y=x2, size=sigmasq2), alpha=0.5) +
  geom_segment(data = simulation$df, aes(x=z1, y=z2, xend=x1, yend=x2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = simulation$df, size=1, shape=2, mapping=aes(x=mu1, y=mu2)) +
  geom_point(data = mean_latents, mapping=aes(x=zhat1, y=zhat2), colour = "green", shape="square") + 
  geom_segment(data = cbind(simulation$df, mean_latents), aes(x=x1, y=x2, xend=zhat1, yend=zhat2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_density_2d(data = clusterMeans, mapping=aes(x=X1, y=X2), colour="goldenrod2") +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        text=element_text(size=15))
ggsave("plots/mean_shift_latents_clustermeans_mu1.png", width=3, height=3, units="in", dpi=600)

ggplot(NULL) +
  geom_point(data = simulation$df, aes(x=z1, y=z2), size=2, shape=1) +
  geom_point(data = simulation$df, mapping=aes(x=x1, y=x2, size=sigmasq2), alpha=0.5) +
  geom_segment(data = simulation$df, aes(x=z1, y=z2, xend=x1, yend=x2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = simulation$df, size=1, shape=2, mapping=aes(x=mu1, y=mu2)) +
  geom_point(data = mean_latents, mapping=aes(x=zhat1, y=zhat2), colour = "green", shape="square") + 
  geom_segment(data = cbind(simulation$df, mean_latents), aes(x=x1, y=x2, xend=zhat1, yend=zhat2), arrow=arrow(length = unit(0.01, "npc")), colour="grey") +
  geom_point(data = empirical_x_means, mapping=aes(x=x1, y=x2), colour="deepskyblue", size=10, shape="x") +
  geom_density_2d(data = clusterMeans, mapping=aes(x=X1, y=X2), colour="goldenrod2") +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        text=element_text(size=15))
ggsave("plots/mean_shift_latents_clustermeans_empirical_mu1.png", width=3, height=3, units="in", dpi=600)
