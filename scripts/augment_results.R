library(cluster)
library(mclust)
library(mcclust)

#args = commandArgs(trailingOnly=TRUE)

#rda_file = paste0("simplex_results/", args[[1]], "/results.rda")
new_snakemake = snakemake
rda_file = snakemake@input[["rda"]]

load(rda_file)
snakemake = new_snakemake

cut_at_k <- function(dpmunc_result, k) {
    psm = dpmunc_result$bigpsm
    hclust.comp = hclust(as.dist(1 - psm), method="complete")
    calls_at_true_k = cutree(hclust.comp, k=k)
}

augmented_results = rbind(results,
      produce_summary("mclust", "x", Mclust(simulation$obsData, x=mclustBIC(simulation$obsData), G=true_k)$classification),
      produce_summary("mclust", "z", Mclust(scaled_latents, x=mclustBIC(scaled_latents), G=true_k)$classification),
      produce_summary("DPMUnc", "x", cut_at_k(result, true_k)),
      produce_summary("DPMUnc_novar", "z", cut_at_k(result_latents, true_k)),
      produce_summary("DPMUnc_novar", "x", cut_at_k(result_novar, true_k)))

augmented_results$inferred_k = TRUE
augmented_results[augmented_results$method == "kmeans_true", ]$inferred_k = FALSE
augmented_results[10:14, ]$inferred_k = FALSE

augmented_results[augmented_results$method == "kmeans_true", ]$method = "kmeans"

print(augmented_results)
write.csv(apply(augmented_results, MARGIN=2, FUN=as.character), file=snakemake@output[["summary"]], row.names=FALSE)
print("Saved to")
print(snakemake@output[["summary"]])
