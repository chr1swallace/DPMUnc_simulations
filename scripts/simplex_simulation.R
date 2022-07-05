library(cluster)
library(DPMUnc)
library(mclust)
library(mcclust)

source("scripts/utils.R")

# Cut a DPMUnc result at a specific K
cut_at_k <- function(dpmunc_result, k) {
    psm = dpmunc_result$bigpsm
    hclust.comp = hclust(as.dist(1 - psm), method="complete")
    calls_at_true_k = cutree(hclust.comp, k=k)
}

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
# Set hyperparameters for model
alpha0 = 2; kappa0 = 0.5; beta0 = 0.2 * mean(apply(simulation$obsData, 2, var))
true_k = length(unique(simulation$df$class))

produce_summary <- function(method, x_or_z, clusters, inferred_K=TRUE) {
    print(paste0("Summarising results from method ", method))
    print(class(clusters))
    print(clusters)
    print(length(clusters))
    print(length(simulation$df$class))
    print(dim(simulation$df))
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
         var_latents=var_latents,
         inferred_K=inferred_K)
}

latents = as.matrix(simulation$df[, paste0("z", 1:d)])

gap_stat <- clusGap(simulation$obsData, FUN = kmeans, nstart = 25, K.max = k * 2, B = 50)
best_K = which.max(gap_stat$Tab[, 3])
kmeans_solution = kmeans(simulation$obsData, centers=best_K, nstart=25)
summary_kmeans = produce_summary("kmeans", "x", kmeans_solution$cluster)

gap_stat <- clusGap(latents, FUN = kmeans, nstart = 25, K.max = k * 2, B = 50)
best_K = which.max(gap_stat$Tab[, 3])
kmeans_solution_latents = kmeans(latents, centers=best_K, nstart=25)
summary_kmeans_latents = produce_summary("kmeans", "z", kmeans_solution_latents$cluster)

kmeans_solution_true = kmeans(simulation$obsData, centers=true_k, nstart=25)
summary_kmeans_true = produce_summary("kmeans", "x", kmeans_solution_true$cluster, FALSE)

kmeans_solution_true_latents = kmeans(latents, centers=true_k, nstart=25)
summary_kmeans_true_latents = produce_summary("kmeans", "z", kmeans_solution_true_latents$cluster, FALSE)

mclust_solution <- Mclust(simulation$obsData, modelNames=c("VVI"))
summary_mclust = produce_summary("mclust_VVI", "x", mclust_solution$classification)

mclust_solution_latents <- Mclust(latents, modelNames=c("VVI"))
summary_mclust_latents = produce_summary("mclust_VVI", "z", mclust_solution_latents$classification)

summary_mclust_true = produce_summary("mclust_EII", "x", Mclust(simulation$obsData, modelNames=c("EII"))$classification)
summary_mclust_latents_true = produce_summary("mclust_EII", "z", Mclust(latents, modelNames=c("EII"))$classification)

outputdir = dirname(snakemake@output[["clusterAllocations"]])
print(dim(simulation$obsData))
DPMUnc(simulation$obsData, simulation$obsVars, saveFileDir = outputdir, seed=seed,
       kappa0=kappa0, alpha0=alpha0, beta0=beta0,
       nIts=10000, scaleData=FALSE)
result = calc_psms(outputdir)
calls=maxpear(result$bigpsm, method="comp")
print(dim(result$bigpsm))
summary_dpmunc = produce_summary("DPMUnc", "x", calls$cl)

outputdir = dirname(snakemake@output[["clusterAllocationsNovar"]])
DPMUnc(simulation$obsData, simulation$obsVars * 1e-10, saveFileDir = outputdir, seed=seed,
       kappa0=kappa0, alpha0=alpha0, beta0=beta0,
       nIts=10000, scaleData=FALSE)
result_novar = calc_psms(outputdir)
calls_novar=maxpear(result_novar$bigpsm, method="comp")
summary_dpmuncnovar = produce_summary("DPMUnc_novar", "x", calls_novar$cl)

outputdir = dirname(snakemake@output[["clusterAllocationsLatents"]])
DPMUnc(latents, simulation$obsVars * 1e-10 , saveFileDir = outputdir, seed=seed,
       kappa0=kappa0, alpha0=alpha0, beta0=beta0,
       nIts=10000, scaleData=FALSE)
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
                                     summary_mclust_true,
                                     summary_mclust_latents_true,
                                     summary_dpmunc,
                                     produce_summary("DPMUnc", "x", cut_at_k(result, true_k), FALSE),
                                     summary_dpmuncnovar,
                                     summary_dpmuncnovar_latents,
                                     produce_summary("DPMUnc_novar", "x", cut_at_k(result_novar, true_k), FALSE),
                                     produce_summary("DPMUnc_novar", "z", cut_at_k(result_latents, true_k), FALSE))))
results
print(class(results))

write.csv(apply(results, MARGIN=2, FUN=as.character), snakemake@output[["summary"]])
save(list=ls(), file=snakemake@output[["rda"]])
