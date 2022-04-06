library(cluster)
library(DPMUnc)
library(mclust)
library(mcclust)

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

seed=as.numeric(snakemake@wildcards[['seed']])
d=as.numeric(snakemake@wildcards[['d']])
var_latents=as.numeric(snakemake@wildcards[['var_latents']])
noise_factor=as.numeric(snakemake@wildcards[['noise_factor']])
k = d + 1
n = k * 10

set.seed(seed)
group_means = data.frame(n_simplex(d))
simulation = generate_basic_uncertain_data(n=n, d=d, k=k, var_latents=var_latents, noise_factor=noise_factor,
                                           group_means=group_means * 10)
true_k = length(unique(simulation$df$class))

produce_summary <- function(method, x_or_z, clusters) {
    list(method=method,
         data=x_or_z,
         estimated_k=length(unique(clusters)),
         ari=adjustedRandIndex(clusters, simulation$df$class),
         seed=seed,
         n=n,
         d=d,
         k=k,
         true_k=true_k,
         noise_factor=noise_factor,
         var_latents=var_latents)
}

scaled_latents = simulation$scale_to_var1(simulation$df[, paste0("z", 1:d)])

gap_stat <- clusGap(simulation$obsData,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = k * 2,
                    B = 50)
print(gap_stat)
best_K = which.max(gap_stat$Tab[, 3])
kmeans_solution = kmeans(simulation$obsData, centers=best_K, nstart=25)
summary_kmeans = produce_summary("kmeans", "x", kmeans_solution$cluster)

kmeans_solution_latents = kmeans(scaled_latents, centers=best_K, nstart=25)
summary_kmeans_latents = produce_summary("kmeans", "z", kmeans_solution_latents$cluster)

kmeans_solution_true = kmeans(simulation$obsData, centers=true_k, nstart=25)
summary_kmeans_true = produce_summary("kmeans_true", "x", kmeans_solution_true$cluster)

kmeans_solution_true_latents = kmeans(scaled_latents, centers=true_k, nstart=25)
summary_kmeans_true_latents = produce_summary("kmeans_true", "z", kmeans_solution_true_latents$cluster)

mclust_solution <- Mclust(simulation$obsData, x=mclustBIC(simulation$obsData))
summary_mclust = produce_summary("mclust", "x", mclust_solution$classification)

mclust_solution_latents <- Mclust(scaled_latents, x=mclustBIC(scaled_latents))
summary_mclust_latents = produce_summary("mclust", "z", mclust_solution_latents$classification)

outputdir = dirname(snakemake@output[["clusterAllocations"]])
DPMUnc(simulation$obsData, simulation$obsVars, saveFileDir = outputdir, seed=seed, nIts=10000)
result = calc_psms(outputdir)
calls=maxpear(result$bigpsm, method="comp")
summary_dpmunc = produce_summary("DPMUnc", "x", calls$cl)

outputdir = dirname(snakemake@output[["clusterAllocationsNovar"]])
DPMUnc(simulation$obsData, simulation$obsVars / 100000, saveFileDir = outputdir, seed=seed, nIts=10000)
result_novar = calc_psms(outputdir)
calls_novar=maxpear(result_novar$bigpsm, method="comp")
summary_dpmuncnovar = produce_summary("DPMUnc_novar", "x", calls_novar$cl)

outputdir = dirname(snakemake@output[["clusterAllocationsLatents"]])
DPMUnc(scaled_latents, simulation$obsVars / 100000, saveFileDir = outputdir, seed=seed, nIts=10000)
result_latents = calc_psms(outputdir)
calls_latents=maxpear(result_latents$bigpsm, method="comp")
summary_dpmuncnovar_latents = produce_summary("DPMUnc_novar", "z", calls_latents$cl)

results = as.data.frame(do.call(rbind,
                                list(summary_kmeans,
                                     summary_kmeans_latents,
                                     summary_kmeans_true,
                                     summary_kmeans_true_latents,
                                     summary_mclust,
                                     summary_mclust_latents,
                                     summary_dpmunc,
                                     summary_dpmuncnovar,
                                     summary_dpmuncnovar_latents)))
results
print(class(results))

write.csv(apply(results, MARGIN=2, FUN=as.character), snakemake@output[["summary"]])
save(list=ls(), file=snakemake@output[["rda"]])
