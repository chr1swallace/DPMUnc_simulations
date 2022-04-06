library(DPMUnc)
library(mclust)
library(mcclust)

source("scripts/utils.R")

generate_basic_uncertain_data <- function(n, d, k, var_latents, var_means, noise_factor) {
  classes = sample(1:k, n, replace=TRUE)
  simulation_df = data.frame(class = classes, sample_id = paste0("sample_", 1:n))

  obsVars_df = data.frame(matrix(rchisq(n*d,1), nrow=n, ncol=d)) * noise_factor
  colnames(obsVars_df) = paste0("sigmasq", 1:d)
  simulation_df = cbind(simulation_df, obsVars_df)

  group_means = data.frame(matrix(rnorm(d*k, mean=0, sd=sqrt(var_means)), nrow=k, ncol=d))
  colnames(group_means) = paste0("mu", 1:d)
  group_means$class = 1:k

  simulation_df = merge(simulation_df, group_means, by="class")

  true_noise = data.frame(matrix(rnorm(n*d, sd=sqrt(var_latents)), nrow=n, ncol=d))
  colnames(true_noise) = paste0("z", 1:d, "_minus_mu", 1:d)
  latents = simulation_df[, paste0("mu", 1:d)] + true_noise
  colnames(latents) = paste0("z", 1:d)
  simulation_df = cbind(simulation_df, true_noise)
  simulation_df = cbind(simulation_df, latents)
  observation_noise = data.frame(matrix(rnorm(n*d, mean=0, sd=sqrt(as.vector(as.matrix(obsVars_df)))),
                                        ncol=d, nrow=n))
  colnames(observation_noise) = paste0("x", 1:d, "_minus_z", 1:d)
  observed = observation_noise + latents
  colnames(observed) = paste0("x", 1:d)

  simulation_df = cbind(simulation_df, observation_noise)
  simulation_df = cbind(simulation_df, observed)

  rownames(simulation_df) = simulation_df$sample_id
  head(simulation_df)

  obsData = as.matrix(simulation_df[, paste0("x", 1:d)])
  obsVars = as.matrix(simulation_df[, paste0("sigmasq", 1:d)])

  scaled = DPMUnc::scale_data(obsData, obsVars)

  scale_to_var1 <- function(mat) {
    scale(mat, attr(scaled$data, "scaled:center"), attr(scaled$data, "scaled:scale"))
  }

  obsData = as.matrix(scaled$data)
  obsVars = as.matrix(scaled$vars)

  pca = prcomp(obsData)
  print(summary(pca))
  pca_df = data.frame(pca$x)
  colnames(pca_df) = paste0("xPC", 1:d)

  rescaled_latents = data.frame(scale(latents, center = pca$center, scale=pca$scale) %*% pca$rotation)
  colnames(rescaled_latents) = paste0("zPC", 1:d)

  rescaled_centers = data.frame(scale(simulation_df[, paste0("mu", 1:d)],
                                      center = pca$center, scale=pca$scale) %*% pca$rotation)
  colnames(rescaled_centers) = paste0("muPC", 1:d)

  simulation_df = cbind(simulation_df, pca_df, rescaled_latents, rescaled_centers)

  head(simulation_df)

  return (list(df=simulation_df,
               pca=pca,
               obsVars=obsVars,
               obsData=obsData,
               scale_to_var1=scale_to_var1))
}

n=as.numeric(snakemake@wildcards[['n']])
seed=as.numeric(snakemake@wildcards[['seed']])
d=as.numeric(snakemake@wildcards[['d']])
k=as.numeric(snakemake@wildcards[['k']])
var_latents=as.numeric(snakemake@wildcards[['var_latents']])
var_means=as.numeric(snakemake@wildcards[['var_means']])
noise_factor=as.numeric(snakemake@wildcards[['noise_factor']])

set.seed(seed)
simulation = generate_basic_uncertain_data(n=n, d=d, k=k, var_latents=var_latents, var_means=var_means, noise_factor=noise_factor)

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#c4a4a8")

scaled_latents = simulation$scale_to_var1(simulation$df[, paste0("z", 1:d)])

kmeans_solution = kmeans(simulation$obsData, centers=k)
kmeans_solution_latents = kmeans(scaled_latents, centers=k)

mclust_solution <- Mclust(simulation$obsData, x=mclustBIC(simulation$obsData))
mclust_solution_latents <- Mclust(scaled_latents, x=mclustBIC(scaled_latents))

outputdir = dirname(snakemake@output[["clusterAllocations"]])
DPMUnc(simulation$obsData, simulation$obsVars, saveFileDir = outputdir, seed=seed, nIts=10000)
result = calc_psms(outputdir)
calls=maxpear(result$bigpsm, method="comp")
outputdir = dirname(snakemake@output[["clusterAllocationsNovar"]])
DPMUnc(simulation$obsData, simulation$obsVars / 100000, saveFileDir = outputdir, seed=seed, nIts=10000)
result_novar = calc_psms(outputdir)
calls_novar=maxpear(result_novar$bigpsm, method="comp")

outputdir = dirname(snakemake@output[["clusterAllocationsLatents"]])
DPMUnc(scaled_latents, simulation$obsVars / 100000, saveFileDir = outputdir, seed=seed, nIts=10000)
result_latents = calc_psms(outputdir)
calls_latents=maxpear(result_latents$bigpsm, method="comp")

results = data.frame(method = rep(c("kmeans", "mclust", "DPMUnc", "DPMUnc_novar"),
                                  c(2,2,1,2)),
                     data = rep(c("x", "z"), 4)[1:7],
                     ari = c(adjustedRandIndex(kmeans_solution$cluster, simulation$df$class),
                             adjustedRandIndex(kmeans_solution_latents$cluster, simulation$df$class),
                             adjustedRandIndex(mclust_solution$classification, simulation$df$class),
                             adjustedRandIndex(mclust_solution_latents$classification, simulation$df$class),
                             adjustedRandIndex(calls$cl, simulation$df$class),
                             adjustedRandIndex(calls_latents$cl, simulation$df$class),
                             adjustedRandIndex(calls_novar$cl, simulation$df$class)),
                     seed=seed,
                     n=n,
                     d=d,
                     k=k,
                     noise_factor=noise_factor,
                     var_means=var_means,
                     var_latents=var_latents)

results

save(list=ls(), file=snakemake@output[["rda"]])
write.csv(results, snakemake@output[["summary"]])
